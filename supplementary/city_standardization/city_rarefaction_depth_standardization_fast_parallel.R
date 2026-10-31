required_pkgs <- c(
  "data.table",
  "readxl",
  "iNEXT",
  "ggplot2",
  "MASS",
  "emmeans"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(data.table)
library(readxl)
library(iNEXT)
library(ggplot2)
library(MASS)
library(emmeans)

inext_version <- as.character(utils::packageVersion("iNEXT"))


get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)

  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
  }

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, winslash = "/", mustWork = TRUE)))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) >= 1) {
  normalizePath(args[1], winslash = "/", mustWork = TRUE)
} else {
  get_script_dir()
}

rpkm_file <- file.path(work_dir, "RPKM_allvirus_matrix_clean.tsv")
metadata_file <- file.path(work_dir, "metadata.xlsx")
depth_file <- file.path(work_dir, "mapping_rate.tsv")
species_cluster_file <- file.path(work_dir, "species_clusters_clean.tsv")
genus_cluster_file <- file.path(work_dir, "genus_clusters_clean.xlsx")
family_cluster_file <- file.path(work_dir, "family_clusters_clean.xlsx")

input_files <- c(
  rpkm_file,
  metadata_file,
  depth_file,
  species_cluster_file,
  genus_cluster_file,
  family_cluster_file
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

out_dir <- file.path(
  work_dir,
  "city_rarefaction_standardization_fast_parallel"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


seed_main <- 20260819L
estimate_boot <- 1000L
curve_boot <- 1000L
coverage_tolerance <- 5e-4

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1) detected_cores <- 1L
requested_cores <- suppressWarnings(as.integer(Sys.getenv("NCORES", "0")))
if (!is.finite(requested_cores) || requested_cores < 1L) {
  requested_cores <- max(1L, detected_cores - 2L)
}
inext_workers <- min(12L, requested_cores)
model_workers <- min(6L, requested_cores)

save_pdf <- identical(Sys.getenv("SAVE_PDF", "0"), "1")


host_order_override <- NA_character_
host_min_samples_per_city <- 5L
host_min_cities <- 5L

depth_support_quantiles <- c(0.05, 0.95)

all_color <- "#4E79A7"
ge2_color <- "#D3724F"
curve_color <- "#4E79A7"
ci_fill <- "#DCE6F2"

normalize_id <- function(x) {
  tolower(trimws(as.character(x)))
}

safe_quantile <- function(x, p) {
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

write_csv <- function(x, name) {
  write.csv(
    as.data.frame(x),
    file.path(out_dir, name),
    row.names = FALSE,
    na = ""
  )
}

save_plot <- function(plot_object, filename, width, height) {
  ggsave(
    file.path(out_dir, paste0(filename, ".png")),
    plot_object,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE
  )
  if (save_pdf) {
    ggsave(
      file.path(out_dir, paste0(filename, ".pdf")),
      plot_object,
      width = width,
      height = height,
      limitsize = FALSE
    )
  }
}

get_set_labels <- function(level) {
  if (level == "Species") {
    return(c(All = "All vOTUs", GE2 = "vOTUs represented by >=2 sequences"))
  }
  if (level == "Genus") {
    return(c(All = "All genera", GE2 = "Genera represented by >=2 vOTUs"))
  }
  c(All = "All families", GE2 = "Families represented by >=2 vOTUs")
}

extract_first_column <- function(x, candidates, new_name) {
  old_name <- intersect(candidates, names(x))
  if (length(old_name) == 0) {
    stop(
      "Could not find any of these columns: ",
      paste(candidates, collapse = ", ")
    )
  }
  x[[new_name]] <- as.numeric(x[[old_name[1]]])
  x
}


read_species_clusters <- function(file) {
  x <- data.table::fread(
    file,
    sep = "\t",
    header = FALSE,
    data.table = FALSE,
    fill = TRUE,
    select = 1:2
  )

  names(x) <- c("Representative", "Members")
  x$Representative <- trimws(as.character(x$Representative))
  x$Members <- as.character(x$Members)
  x <- x[!is.na(x$Representative) & x$Representative != "", , drop = FALSE]

  if (anyDuplicated(x$Representative)) {
    stop("Duplicate species representative IDs were found.")
  }

  x$Cluster_size <- ifelse(
    is.na(x$Members) | trimws(x$Members) == "",
    0L,
    lengths(strsplit(x$Members, ",", fixed = TRUE))
  )
  x
}

clean_character_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "character"
  x <- apply(x, 2, function(z) trimws(as.character(z)))
  x[x == ""] <- NA_character_
  x
}

read_higher_clusters <- function(file, level_name) {
  raw <- readxl::read_xlsx(
    file,
    col_names = FALSE,
    .name_repair = "minimal"
  )
  mat <- clean_character_matrix(raw)

  row_members <- lapply(seq_len(nrow(mat)), function(i) {
    unique(mat[i, !is.na(mat[i, ]) & mat[i, ] != ""])
  })

  representative <- vapply(
    row_members,
    function(z) if (length(z) == 0) NA_character_ else z[1],
    character(1)
  )

  keep <- !is.na(representative) & representative != ""
  representative <- representative[keep]
  row_members <- row_members[keep]

  split_rows <- split(seq_along(representative), representative)
  merged_members <- lapply(split_rows, function(ii) {
    unique(unlist(row_members[ii], use.names = FALSE))
  })

  duplicate_rows_merged <- length(representative) - length(split_rows)
  if (duplicate_rows_merged > 0) {
  }

  list(
    representative = names(split_rows),
    members = unname(merged_members),
    cluster_size = lengths(merged_members),
    duplicate_rows_merged = duplicate_rows_merged
  )
}

aggregate_species_pa_to_level <- function(pa_species, cluster_object, level_name) {
  matched_members <- lapply(cluster_object$members, function(z) {
    intersect(z, rownames(pa_species))
  })

  keep <- lengths(matched_members) > 0
  if (!any(keep)) {
    stop("No ", level_name, " clusters matched species-level vOTUs.")
  }

  matched_members <- matched_members[keep]
  representatives <- cluster_object$representative[keep]
  source_cluster_size <- cluster_object$cluster_size[keep]

  out <- matrix(
    0L,
    nrow = length(representatives),
    ncol = ncol(pa_species),
    dimnames = list(representatives, colnames(pa_species))
  )

  for (i in seq_along(representatives)) {
    out[i, ] <- as.integer(
      colSums(pa_species[matched_members[[i]], , drop = FALSE]) > 0
    )
  }

  list(
    pa = out,
    source_cluster_size = source_cluster_size,
    matched_member_count = lengths(matched_members),
    n_unmatched_clusters = sum(!keep)
  )
}


meta <- as.data.frame(readxl::read_excel(metadata_file))
required_meta <- c("Library", "City", "Host_order")
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0) {
  stop("metadata.xlsx is missing: ", paste(missing_meta, collapse = ", "))
}

meta$sample_key <- normalize_id(meta$Library)
meta$City <- trimws(as.character(meta$City))
meta$Host_order <- trimws(as.character(meta$Host_order))
meta$Host_order[is.na(meta$Host_order) | meta$Host_order == ""] <- "Unknown"

if (anyDuplicated(meta$sample_key)) {
  stop("Duplicate Library names were found after normalization.")
}
if (any(is.na(meta$City) | meta$City == "")) {
  stop("Missing City values were found in metadata.xlsx.")
}

depth <- data.table::fread(depth_file, data.table = FALSE)
required_depth <- c("sample", "total_reads")
missing_depth <- setdiff(required_depth, names(depth))
if (length(missing_depth) > 0) {
  stop("mapping_rate.tsv is missing: ", paste(missing_depth, collapse = ", "))
}

depth$sample_key <- normalize_id(depth$sample)
depth$sequencing_depth <- as.numeric(depth$total_reads)
if (anyDuplicated(depth$sample_key)) {
  stop("Duplicate sample names were found in mapping_rate.tsv.")
}
if (any(!is.finite(depth$sequencing_depth) | depth$sequencing_depth <= 0)) {
  stop("All total_reads values must be finite and positive.")
}

rpkm <- data.table::fread(
  rpkm_file,
  data.table = TRUE,
  check.names = FALSE
)
if (ncol(rpkm) < 2) {
  stop("RPKM matrix must contain one vOTU column and sample columns.")
}
data.table::setnames(rpkm, names(rpkm)[1], "vOTU")
rpkm$vOTU <- trimws(as.character(rpkm$vOTU))
if (anyDuplicated(rpkm$vOTU)) {
  stop("Duplicate vOTU rows were found in the RPKM matrix.")
}

sample_cols <- names(rpkm)[-1]
sample_keys <- normalize_id(sample_cols)
if (anyDuplicated(sample_keys)) {
  stop("Duplicate RPKM sample columns were found after normalization.")
}

rpkm_name_by_key <- setNames(sample_cols, sample_keys)
species_pa <- as.matrix(rpkm[, ..sample_cols]) > 0
storage.mode(species_pa) <- "integer"
rownames(species_pa) <- rpkm$vOTU
colnames(species_pa) <- sample_cols
rm(rpkm)
invisible(gc())


species_clusters <- read_species_clusters(species_cluster_file)
species_ge2_ids <- intersect(
  species_clusters$Representative[species_clusters$Cluster_size >= 2],
  rownames(species_pa)
)
if (length(species_ge2_ids) == 0) {
  stop("No Species cluster_size >=2 representatives matched the RPKM matrix.")
}

genus_clusters <- read_higher_clusters(genus_cluster_file, "Genus")
genus_agg <- aggregate_species_pa_to_level(species_pa, genus_clusters, "Genus")

family_clusters <- read_higher_clusters(family_cluster_file, "Family")
family_agg <- aggregate_species_pa_to_level(species_pa, family_clusters, "Family")

level_matrices <- list(
  Species = list(
    All = species_pa,
    GE2 = species_pa[species_ge2_ids, , drop = FALSE]
  ),
  Genus = list(
    All = genus_agg$pa,
    GE2 = genus_agg$pa[
      genus_agg$source_cluster_size >= 2,
      ,
      drop = FALSE
    ]
  ),
  Family = list(
    All = family_agg$pa,
    GE2 = family_agg$pa[
      family_agg$source_cluster_size >= 2,
      ,
      drop = FALSE
    ]
  )
)

if (any(vapply(level_matrices, function(z) {
  any(vapply(z, nrow, integer(1)) == 0)
}, logical(1)))) {
  stop("At least one taxonomic level/set matrix is empty.")
}


common_meta_rpkm <- intersect(meta$sample_key, sample_keys)
common_all3 <- Reduce(
  intersect,
  list(meta$sample_key, sample_keys, depth$sample_key)
)

if (length(common_meta_rpkm) == 0) {
  stop("No samples matched metadata and RPKM.")
}
if (length(common_all3) != length(common_meta_rpkm)) {
  stop(
    "Sequencing depth is missing for ",
    length(common_meta_rpkm) - length(common_all3),
    " metadata/RPKM-matched sample(s)."
  )
}

meta_acc <- meta[meta$sample_key %in% common_all3, , drop = FALSE]
meta_acc$rpkm_sample <- unname(rpkm_name_by_key[meta_acc$sample_key])

city_counts <- as.data.frame(table(meta_acc$City), stringsAsFactors = FALSE)
names(city_counts) <- c("City", "N_samples")
city_counts <- city_counts[order(city_counts$N_samples, decreasing = TRUE), ]

host_city_table <- table(meta_acc$Host_order, meta_acc$City)
eligible_city_count <- rowSums(host_city_table >= host_min_samples_per_city)
eligible_sample_count <- vapply(seq_len(nrow(host_city_table)), function(i) {
  sum(host_city_table[i, host_city_table[i, ] >= host_min_samples_per_city])
}, numeric(1))

host_candidates <- data.frame(
  Host_order = rownames(host_city_table),
  Eligible_cities = as.integer(eligible_city_count),
  Eligible_samples = as.integer(eligible_sample_count),
  stringsAsFactors = FALSE
)
host_candidates <- host_candidates[
  host_candidates$Host_order != "Unknown" &
    host_candidates$Eligible_cities >= host_min_cities,
  ,
  drop = FALSE
]

if (nrow(host_candidates) == 0) {
  stop(
    "No Host_order has at least ",
    host_min_samples_per_city,
    " samples in at least ",
    host_min_cities,
    " cities."
  )
}

if (!is.na(host_order_override)) {
  if (!host_order_override %in% host_candidates$Host_order) {
    stop("host_order_override does not meet the common-support thresholds.")
  }
  selected_host <- host_order_override
} else {
  host_candidates <- host_candidates[
    order(
      host_candidates$Eligible_cities,
      host_candidates$Eligible_samples,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  selected_host <- host_candidates$Host_order[1]
}

selected_host_counts <- host_city_table[selected_host, ]
selected_host_cities <- names(selected_host_counts)[
  selected_host_counts >= host_min_samples_per_city
]

meta_host <- meta_acc[
  meta_acc$Host_order == selected_host &
    meta_acc$City %in% selected_host_cities,
  ,
  drop = FALSE
]

host_selection_summary <- data.frame(
  Selected_Host_order = selected_host,
  Minimum_samples_per_city = host_min_samples_per_city,
  Cities_included = length(selected_host_cities),
  Samples_included = nrow(meta_host),
  Cities_excluded = paste(
    setdiff(sort(unique(meta_acc$City)), sort(selected_host_cities)),
    collapse = "; "
  ),
  Season_in_model = "No",
  stringsAsFactors = FALSE
)

write_csv(city_counts, "city_sample_counts.csv")
write_csv(host_candidates, "host_order_support_candidates.csv")
write_csv(host_selection_summary, "host_order_sensitivity_selection.csv")

matching_summary <- data.frame(
  Item = c(
    "Samples in metadata",
    "Samples in RPKM matrix",
    "Samples in mapping_rate.tsv",
    "Matched across all three files",
    "Cities in full analysis",
    "Cities in selected Host_order sensitivity analysis",
    "Species All",
    "Species >=2 sequences",
    "Genus All after duplicate merge and unmatched-cluster removal",
    "Genus >=2 vOTUs",
    "Family All after duplicate merge and unmatched-cluster removal",
    "Family >=2 vOTUs",
    "Duplicated genus representative rows merged",
    "Duplicated family representative rows merged"
  ),
  N = c(
    nrow(meta),
    length(sample_keys),
    nrow(depth),
    length(common_all3),
    length(unique(meta_acc$City)),
    length(selected_host_cities),
    nrow(level_matrices$Species$All),
    nrow(level_matrices$Species$GE2),
    nrow(level_matrices$Genus$All),
    nrow(level_matrices$Genus$GE2),
    nrow(level_matrices$Family$All),
    nrow(level_matrices$Family$GE2),
    genus_clusters$duplicate_rows_merged,
    family_clusters$duplicate_rows_merged
  ),
  stringsAsFactors = FALSE
)
write_csv(matching_summary, "matching_summary.csv")


run_inext_analysis <- function(
    pa_matrix,
    meta_use,
    level,
    set_name,
    scope,
    scope_stub,
    seed_offset
) {
  labels <- get_set_labels(level)
  group_label <- unname(labels[set_name])
  cities <- sort(unique(meta_use$City))
  file_base <- paste(
    tolower(level),
    tolower(set_name),
    scope_stub,
    sep = "_"
  )

  city_n <- as.data.frame(table(meta_use$City), stringsAsFactors = FALSE)
  names(city_n) <- c("City", "Reference_samples")
  n_min <- min(city_n$Reference_samples)

  if (n_min < 2) {
    stop(scope, " contains a city with fewer than two samples.")
  }

  incidence_list <- setNames(vector("list", length(cities)), cities)
  for (city in cities) {
    city_samples <- meta_use$rpkm_sample[meta_use$City == city]
    city_samples <- intersect(city_samples, colnames(pa_matrix))
    if (length(city_samples) < n_min) {
      stop("Unexpected sample loss for ", city, " at ", level, "/", set_name)
    }
    incidence_list[[city]] <- iNEXT::as.incfreq(
      pa_matrix[, city_samples, drop = FALSE]
    )
  }

  set.seed(seed_main + seed_offset)
  size_est <- iNEXT::estimateD(
    incidence_list,
    q = 0,
    datatype = "incidence_freq",
    base = "size",
    level = n_min,
    nboot = estimate_boot,
    conf = 0.95
  )
  size_est <- extract_first_column(
    size_est,
    c("m", "t"),
    "Returned_sample_size"
  )
  size_est$Target_sample_size <- n_min

  if (any(abs(size_est$Returned_sample_size - n_min) > 1e-6)) {
    stop("iNEXT did not return the requested common sample size.")
  }

  coverage_limit_list <- lapply(cities, function(city) {
    reference_n <- city_n$Reference_samples[
      match(city, city_n$City)
    ]

    at_bounds <- iNEXT::estimateD(
      setNames(list(incidence_list[[city]]), city),
      q = 0,
      datatype = "incidence_freq",
      base = "size",
      level = c(1, 2 * reference_n),
      nboot = 0,
      conf = 0.95
    )
    at_bounds <- extract_first_column(
      at_bounds,
      c("m", "t"),
      "Returned_sample_size"
    )

    one_row <- which(abs(at_bounds$Returned_sample_size - 1) <= 1e-6)
    double_row <- which(
      abs(at_bounds$Returned_sample_size - 2 * reference_n) <= 1e-6
    )
    if (length(one_row) != 1L || length(double_row) != 1L) {
      stop("Unexpected SC(1)/SC(2n) output for city ", city, ".")
    }

    bound_sc <- as.numeric(at_bounds$SC[c(one_row, double_row)])
    if (any(!is.finite(bound_sc)) || any(bound_sc <= 0) || any(bound_sc > 1)) {
      stop("Invalid SC(1)/SC(2n) returned for city ", city, ".")
    }
    if (bound_sc[1] > bound_sc[2] + coverage_tolerance) {
      stop("Non-monotonic coverage bounds returned for city ", city, ".")
    }

    data.frame(
      Assemblage = city,
      Reference_samples = as.numeric(reference_n),
      Double_reference_sample_size = 2 * as.numeric(reference_n),
      Coverage_at_one_sample = bound_sc[1],
      Coverage_at_double_reference = bound_sc[2],
      stringsAsFactors = FALSE
    )
  })

  coverage_limits <- data.table::rbindlist(
    coverage_limit_list,
    use.names = TRUE,
    fill = TRUE
  )

  write_csv(
    coverage_limits,
    paste0(file_base, "_coverage_attainable_bounds.csv")
  )

  common_coverage_lower_bound <- max(
    coverage_limits$Coverage_at_one_sample
  )
  common_coverage_upper_bound <- min(
    coverage_limits$Coverage_at_double_reference
  )

  lower_city <- coverage_limits$Assemblage[
    which.max(coverage_limits$Coverage_at_one_sample)
  ]
  upper_city <- coverage_limits$Assemblage[
    which.min(coverage_limits$Coverage_at_double_reference)
  ]

  coverage_standardization_available <-
    common_coverage_lower_bound <=
      common_coverage_upper_bound + coverage_tolerance

  coverage_est <- data.frame(
    Assemblage = character(0),
    SC = numeric(0),
    Equivalent_sample_size = numeric(0),
    qD = numeric(0),
    qD.LCL = numeric(0),
    qD.UCL = numeric(0),
    Level = character(0),
    Set = character(0),
    Group = character(0),
    Scope = character(0),
    Selected_Host_order = character(0),
    stringsAsFactors = FALSE
  )
  target_coverage <- NA_real_
  coverage_spread <- NA_real_
  coverage_target_max_abs_difference <- NA_real_
  maximum_extrapolation_multiple <- NA_real_
  coverage_validation_pass <- NA
  two_times_limit_pass <- NA
  coverage_note <- ""

  if (coverage_standardization_available) {
    target_coverage <- common_coverage_upper_bound


    set.seed(seed_main + seed_offset + 1L)
    coverage_est <- iNEXT::estimateD(
      incidence_list,
      q = 0,
      datatype = "incidence_freq",
      base = "coverage",
      level = target_coverage,
      nboot = estimate_boot,
      conf = 0.95
    )
    coverage_est <- extract_first_column(
      coverage_est,
      c("t", "m"),
      "Equivalent_sample_size"
    )

    coverage_spread <- max(coverage_est$SC) - min(coverage_est$SC)
    coverage_target_max_abs_difference <- max(
      abs(coverage_est$SC - target_coverage)
    )

    coverage_est <- merge(
      coverage_est,
      coverage_limits,
      by = "Assemblage",
      all.x = TRUE,
      sort = FALSE
    )
    coverage_est$Common_target_coverage <- target_coverage
    coverage_est$Coverage_difference <- coverage_est$SC - target_coverage
    coverage_est$Extrapolation_multiple <-
      coverage_est$Equivalent_sample_size / coverage_est$Reference_samples

    maximum_extrapolation_multiple <- max(
      coverage_est$Extrapolation_multiple,
      na.rm = TRUE
    )
    coverage_validation_pass <-
      coverage_target_max_abs_difference <= coverage_tolerance
    two_times_limit_pass <- all(
      coverage_est$Extrapolation_multiple <= 2 + coverage_tolerance
    )

    if (!coverage_validation_pass) {
      stop(
        "Coverage standardization was expected to be attainable but iNEXT ",
        "did not return the requested common coverage for ",
        level, "/", set_name, "/", scope,
        ". Target = ", signif(target_coverage, 8),
        "; maximum absolute SC difference = ",
        signif(coverage_target_max_abs_difference, 6), "."
      )
    }
    if (!two_times_limit_pass) {
      stop(
        "Coverage standardization exceeded the two-times-reference limit for ",
        level, "/", set_name, "/", scope, "."
      )
    }

    coverage_note <- "Exact common coverage available and validated."
  } else {
    coverage_note <- paste0(
      "No exact common coverage within 1 to 2n for all assemblages. ",
      "Highest SC(1) = ",
      signif(common_coverage_lower_bound, 8),
      " (", lower_city, "); lowest SC(2n) = ",
      signif(common_coverage_upper_bound, 8),
      " (", upper_city, ")."
    )

  }

  set.seed(seed_main + seed_offset + 2L)
  inext_object <- iNEXT::iNEXT(
    incidence_list,
    q = 0,
    datatype = "incidence_freq",
    endpoint = NULL,
    se = TRUE,
    conf = 0.95,
    nboot = curve_boot
  )

  size_curve <- as.data.frame(inext_object$iNextEst$size_based)
  coverage_curve <- as.data.frame(inext_object$iNextEst$coverage_based)
  size_curve <- extract_first_column(size_curve, c("m", "t"), "Sample_size")

  for (x in list(size_est, size_curve, coverage_curve)) {
    if (!all(c("Assemblage", "qD", "qD.LCL", "qD.UCL") %in% names(x))) {
      stop("Unexpected iNEXT output columns; check the installed iNEXT version.")
    }
  }
  if (
    coverage_standardization_available &&
      !all(c("Assemblage", "qD", "qD.LCL", "qD.UCL") %in% names(coverage_est))
  ) {
    stop("Unexpected coverage-standardization output columns from iNEXT.")
  }

  add_metadata <- function(x) {
    x$Level <- level
    x$Set <- set_name
    x$Group <- group_label
    x$Scope <- scope
    x$Selected_Host_order <- if (scope_stub == "all_cities") {
      NA_character_
    } else {
      selected_host
    }
    x
  }

  size_est <- add_metadata(size_est)
  if (coverage_standardization_available) {
    coverage_est <- add_metadata(coverage_est)
  }
  size_curve <- add_metadata(size_curve)
  coverage_curve <- add_metadata(coverage_curve)


  validation <- data.frame(
    Level = level,
    Set = set_name,
    Group = group_label,
    Scope = scope,
    Cities = length(cities),
    Common_sample_size = n_min,
    Coverage_standardization_available =
      coverage_standardization_available,
    Common_coverage_lower_bound = common_coverage_lower_bound,
    Common_coverage_upper_bound = common_coverage_upper_bound,
    Common_target_coverage = target_coverage,
    Lower_bound_city = lower_city,
    Upper_bound_city = upper_city,
    Coverage_spread = coverage_spread,
    Coverage_target_max_abs_difference =
      coverage_target_max_abs_difference,
    Maximum_extrapolation_multiple = maximum_extrapolation_multiple,
    Coverage_validation_pass = coverage_validation_pass,
    Two_times_limit_pass = two_times_limit_pass,
    Coverage_note = coverage_note,
    stringsAsFactors = FALSE
  )

  list(
    size_est = size_est,
    coverage_est = coverage_est,
    size_curve = size_curve,
    coverage_curve = coverage_curve,
    validation = validation
  )
}


analysis_jobs <- list()
job_index <- 1L
for (level in c("Species", "Genus", "Family")) {
  for (set_name in c("All", "GE2")) {
    analysis_jobs[[job_index]] <- list(
      level = level,
      set_name = set_name,
      scope_type = "all",
      seed_offset = job_index * 20L
    )
    job_index <- job_index + 1L
    analysis_jobs[[job_index]] <- list(
      level = level,
      set_name = set_name,
      scope_type = "host",
      seed_offset = job_index * 20L
    )
    job_index <- job_index + 1L
  }
}

run_analysis_job <- function(job) {
  level <- job$level
  set_name <- job$set_name

  if (job$scope_type == "all") {
    meta_use <- meta_acc
    scope <- "All matched samples; observed host composition"
    scope_stub <- "all_cities"
  } else {
    meta_use <- meta_host
    scope <- paste0(
      "Host-order-standardized sensitivity restricted to ",
      selected_host,
      " (cities with >=",
      host_min_samples_per_city,
      " samples)"
    )
    scope_stub <- paste0(
      "host_",
      gsub("[^A-Za-z0-9]+", "_", tolower(selected_host))
    )
  }


  run_inext_analysis(
    level_matrices[[level]][[set_name]],
    meta_use,
    level,
    set_name,
    scope,
    scope_stub,
    seed_offset = job$seed_offset
  )
}

if (.Platform$OS.type == "unix" && inext_workers > 1L) {
  all_inext_results <- parallel::mclapply(
    analysis_jobs,
    run_analysis_job,
    mc.cores = inext_workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  all_inext_results <- lapply(analysis_jobs, run_analysis_job)
}

failed_jobs <- vapply(all_inext_results, inherits, logical(1), what = "try-error")
if (any(failed_jobs)) {
  stop("At least one parallel iNEXT job failed.")
}

combine_result <- function(component) {
  data.table::rbindlist(
    lapply(all_inext_results, `[[`, component),
    use.names = TRUE,
    fill = TRUE
  )
}

size_est_df <- combine_result("size_est")
coverage_est_df <- combine_result("coverage_est")
size_curve_df <- combine_result("size_curve")
coverage_curve_df <- combine_result("coverage_curve")
coverage_validation_df <- combine_result("validation")

write_csv(size_est_df, "exact_sample_size_standardized_richness.csv")
write_csv(coverage_est_df, "coverage_standardized_richness.csv")
write_csv(size_curve_df, "inext_size_curve_data.csv")
write_csv(coverage_curve_df, "inext_coverage_curve_data.csv")
write_csv(coverage_validation_df, "coverage_target_validation.csv")


method_values <- c(
  Rarefaction = "solid",
  Observed = "solid",
  Extrapolation = "dashed"
)

scope_filter <- function(x, host_scope) {
  if (host_scope) {
    x[!is.na(x$Selected_Host_order), , drop = FALSE]
  } else {
    x[is.na(x$Selected_Host_order), , drop = FALSE]
  }
}

for (level in c("Species", "Genus", "Family")) {
  for (host_scope in c(FALSE, TRUE)) {
    scope_stub <- if (host_scope) {
      paste0("host_", gsub("[^A-Za-z0-9]+", "_", tolower(selected_host)))
    } else {
      "all_cities"
    }
    scope_title <- if (host_scope) {
      paste0("Host-order sensitivity: ", selected_host)
    } else {
      "All matched samples; observed host composition"
    }

    size_dat <- scope_filter(
      size_curve_df[size_curve_df$Level == level, , drop = FALSE],
      host_scope
    )
    size_dat$Set <- factor(size_dat$Set, levels = c("All", "GE2"))

    p_size <- ggplot(
      size_dat,
      aes(
        x = Sample_size,
        y = qD,
        color = Set,
        fill = Set,
        group = interaction(Set, Method)
      )
    ) +
      geom_ribbon(
        aes(ymin = qD.LCL, ymax = qD.UCL),
        alpha = 0.12,
        color = NA
      ) +
      geom_line(
        aes(linetype = Method),
        linewidth = 0.7
      ) +
      geom_point(
        data = size_dat[size_dat$Method == "Observed", , drop = FALSE],
        size = 1.4
      ) +
      facet_wrap(~Assemblage, ncol = 4, scales = "free_y") +
      scale_color_manual(
        values = c(All = all_color, GE2 = ge2_color),
        labels = get_set_labels(level)
      ) +
      scale_fill_manual(
        values = c(All = all_color, GE2 = ge2_color),
        labels = get_set_labels(level)
      ) +
      scale_linetype_manual(values = method_values, drop = FALSE) +
      labs(
        x = "Number of samples",
        y = "Estimated richness (q = 0)",
        color = NULL,
        fill = NULL,
        linetype = NULL,
        title = paste0(level, "-level sample-size-based accumulation curves"),
        subtitle = paste0(scope_title, "; All and >=2 are shown together.")
      ) +
      theme_bw(base_size = 10) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title.position = "plot"
      )

    save_plot(
      p_size,
      paste0(tolower(level), "_all_ge2_", scope_stub, "_size_curves"),
      12,
      15
    )

    coverage_dat <- scope_filter(
      coverage_curve_df[coverage_curve_df$Level == level, , drop = FALSE],
      host_scope
    )
    coverage_dat$Set <- factor(coverage_dat$Set, levels = c("All", "GE2"))

    p_cov <- ggplot(
      coverage_dat,
      aes(
        x = SC,
        y = qD,
        color = Set,
        fill = Set,
        group = interaction(Set, Method)
      )
    ) +
      geom_ribbon(
        aes(ymin = qD.LCL, ymax = qD.UCL),
        alpha = 0.12,
        color = NA
      ) +
      geom_line(
        aes(linetype = Method),
        linewidth = 0.7
      ) +
      geom_point(
        data = coverage_dat[coverage_dat$Method == "Observed", , drop = FALSE],
        size = 1.4
      ) +
      facet_wrap(~Assemblage, ncol = 4, scales = "free_y") +
      scale_color_manual(
        values = c(All = all_color, GE2 = ge2_color),
        labels = get_set_labels(level)
      ) +
      scale_fill_manual(
        values = c(All = all_color, GE2 = ge2_color),
        labels = get_set_labels(level)
      ) +
      scale_linetype_manual(values = method_values, drop = FALSE) +
      scale_x_continuous(limits = c(0, 1)) +
      labs(
        x = "Sample coverage",
        y = "Estimated richness (q = 0)",
        color = NULL,
        fill = NULL,
        linetype = NULL,
        title = paste0(level, "-level coverage-based accumulation curves"),
        subtitle = paste0(scope_title, "; All and >=2 are shown together.")
      ) +
      theme_bw(base_size = 10) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title.position = "plot"
      )

    save_plot(
      p_cov,
      paste0(tolower(level), "_all_ge2_", scope_stub, "_coverage_curves"),
      12,
      15
    )

    point_dat <- scope_filter(
      size_est_df[size_est_df$Level == level, , drop = FALSE],
      host_scope
    )
    point_dat$Set <- factor(point_dat$Set, levels = c("All", "GE2"))
    all_rows <- point_dat[point_dat$Set == "All", , drop = FALSE]
    city_order <- all_rows$Assemblage[order(all_rows$qD)]
    point_dat$Assemblage <- factor(point_dat$Assemblage, levels = city_order)

    p_point <- ggplot(
      point_dat,
      aes(x = Assemblage, y = qD, color = Set)
    ) +
      geom_errorbar(
        aes(ymin = qD.LCL, ymax = qD.UCL),
        width = 0.18,
        position = position_dodge(width = 0.55)
      ) +
      geom_point(
        size = 2,
        position = position_dodge(width = 0.55)
      ) +
      coord_flip() +
      scale_color_manual(
        values = c(All = all_color, GE2 = ge2_color),
        labels = get_set_labels(level)
      ) +
      labs(
        x = NULL,
        y = "Equal-sample-size standardized richness",
        color = NULL,
        title = paste0(level, "-level standardized richness: All vs >=2"),
        subtitle = scope_title
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title.position = "plot"
      )

    save_plot(
      p_point,
      paste0(tolower(level), "_all_ge2_", scope_stub, "_size_standardized"),
      9,
      9
    )
  }
}


sample_richness_list <- list()
sample_index <- 1L

for (level in c("Species", "Genus", "Family")) {
  labels <- get_set_labels(level)
  for (set_name in c("All", "GE2")) {
    pa_matrix <- level_matrices[[level]][[set_name]]
    richness <- colSums(pa_matrix > 0)
    sample_richness_list[[sample_index]] <- data.frame(
      sample_key = normalize_id(names(richness)),
      RPKM_sample = names(richness),
      Level = level,
      Set = set_name,
      Group = unname(labels[set_name]),
      Richness = as.numeric(richness),
      stringsAsFactors = FALSE
    )
    sample_index <- sample_index + 1L
  }
}

sample_richness_df <- data.table::rbindlist(
  sample_richness_list,
  use.names = TRUE,
  fill = TRUE
)

model_dat <- merge(
  meta_host[, c("Library", "City", "Host_order", "sample_key")],
  sample_richness_df,
  by = "sample_key",
  all = FALSE
)
model_dat <- merge(
  model_dat,
  depth[, c("sample_key", "sequencing_depth")],
  by = "sample_key",
  all = FALSE
)
model_dat$log10_depth <- log10(model_dat$sequencing_depth)
model_dat$City <- factor(model_dat$City)

write_csv(
  model_dat,
  "depth_model_data_selected_host_order.csv"
)

depth_sample_dat <- unique(
  model_dat[, c("sample_key", "City", "sequencing_depth", "log10_depth")]
)

depth_interval_list <- lapply(
  split(depth_sample_dat$log10_depth, depth_sample_dat$City),
  function(z) {
    c(
      Lower = safe_quantile(z, depth_support_quantiles[1]),
      Upper = safe_quantile(z, depth_support_quantiles[2]),
      Median = median(z),
      N = length(z)
    )
  }
)

depth_interval <- data.frame(
  City = names(depth_interval_list),
  do.call(rbind, depth_interval_list),
  row.names = NULL,
  stringsAsFactors = FALSE
)

candidate_depths <- sort(unique(c(depth_interval$Lower, depth_interval$Upper)))
pooled_median_depth <- median(depth_sample_dat$log10_depth)
candidate_support <- vapply(candidate_depths, function(z) {
  sum(z >= depth_interval$Lower & z <= depth_interval$Upper)
}, integer(1))

maximum_supported_cities <- max(candidate_support)
best_candidates <- candidate_depths[
  candidate_support == maximum_supported_cities
]
common_log_depth <- best_candidates[
  which.min(abs(best_candidates - pooled_median_depth))
]

depth_model_cities <- depth_interval$City[
  common_log_depth >= depth_interval$Lower &
    common_log_depth <= depth_interval$Upper
]

if (length(depth_model_cities) < host_min_cities) {
  stop(
    "Common sequencing-depth support remains in only ",
    length(depth_model_cities),
    " cities; review depth_support_quantiles or the input data."
  )
}

depth_interval$Selected_for_depth_model <-
  depth_interval$City %in% depth_model_cities
depth_interval$Common_log10_depth <- common_log_depth
depth_interval$Common_sequencing_depth <- 10^common_log_depth
depth_interval$Common_depth_in_central_interval <-
  common_log_depth >= depth_interval$Lower &
    common_log_depth <= depth_interval$Upper
depth_interval$Lower_quantile <- depth_support_quantiles[1]
depth_interval$Upper_quantile <- depth_support_quantiles[2]

write_csv(depth_interval, "sequencing_depth_common_support.csv")

model_dat <- model_dat[
  as.character(model_dat$City) %in% depth_model_cities,
  ,
  drop = FALSE
]
model_dat$City <- droplevels(model_dat$City)

model_jobs <- list()
model_index <- 1L
for (level in c("Species", "Genus", "Family")) {
  for (set_name in c("All", "GE2")) {
    model_jobs[[model_index]] <- list(
      level = level,
      set_name = set_name,
      model_index = model_index
    )
    model_index <- model_index + 1L
  }
}

run_depth_model_job <- function(job) {
  level <- job$level
  set_name <- job$set_name
  labels <- get_set_labels(level)


  dat <- model_dat[
    model_dat$Level == level & model_dat$Set == set_name,
    ,
    drop = FALSE
  ]
  dat$City <- droplevels(dat$City)

  full_model <- MASS::glm.nb(
    Richness ~ City + splines::ns(log10_depth, df = 3),
    data = dat,
    control = glm.control(maxit = 100)
  )
  reduced_model <- MASS::glm.nb(
    Richness ~ splines::ns(log10_depth, df = 3),
    data = dat,
    control = glm.control(maxit = 100)
  )

  if (!isTRUE(full_model$converged) || !isTRUE(reduced_model$converged)) {
    stop("Negative-binomial model did not converge for ", level, "/", set_name)
  }

  lrt <- as.data.frame(anova(reduced_model, full_model, test = "Chisq"))
  lrt$Model <- rownames(lrt)
  rownames(lrt) <- NULL
  lrt$Level <- level
  lrt$Set <- set_name
  lrt$Group <- unname(labels[set_name])
  lrt$Selected_Host_order <- selected_host

  capture.output(
    summary(full_model),
    file = file.path(
      out_dir,
      paste0(
        tolower(level),
        "_",
        tolower(set_name),
        "_depth_model_",
        gsub("[^A-Za-z0-9]+", "_", tolower(selected_host)),
        ".txt"
      )
    )
  )

  emm <- emmeans::emmeans(
    full_model,
    ~ City,
    at = list(log10_depth = common_log_depth),
    type = "response"
  )
  emm_df <- as.data.frame(summary(emm, infer = c(TRUE, TRUE), type = "response"))

  estimate_col <- intersect(c("response", "rate", "emmean"), names(emm_df))
  lower_col <- intersect(c("asymp.LCL", "lower.CL"), names(emm_df))
  upper_col <- intersect(c("asymp.UCL", "upper.CL"), names(emm_df))
  if (length(estimate_col) == 0 || length(lower_col) == 0 || length(upper_col) == 0) {
    stop("Unexpected emmeans output columns for ", level, "/", set_name)
  }

  adjusted <- data.frame(
    Level = level,
    Set = set_name,
    Group = unname(labels[set_name]),
    City = as.character(emm_df$City),
    Selected_Host_order = selected_host,
    Common_sequencing_depth = 10^common_log_depth,
    Adjusted_alpha_richness = as.numeric(emm_df[[estimate_col[1]]]),
    SE = as.numeric(emm_df$SE),
    CI95_lower = as.numeric(emm_df[[lower_col[1]]]),
    CI95_upper = as.numeric(emm_df[[upper_col[1]]]),
    stringsAsFactors = FALSE
  )

  model_support <- depth_interval[
    depth_interval$Selected_for_depth_model,
    ,
    drop = FALSE
  ]
  model_support$Level <- level
  model_support$Set <- set_name
  model_support$Selected_Host_order <- selected_host

  pairwise <- as.data.frame(
    emmeans::contrast(
      emm,
      method = "pairwise",
      adjust = "tukey",
      type = "response"
    )
  )
  pairwise$Level <- level
  pairwise$Set <- set_name
  pairwise$Group <- unname(labels[set_name])
  pairwise$Selected_Host_order <- selected_host
  pairwise$Common_sequencing_depth <- 10^common_log_depth

  list(
    adjusted = adjusted,
    support = model_support,
    global_test = lrt,
    pairwise = pairwise
  )
}

if (.Platform$OS.type == "unix" && model_workers > 1L) {
  model_results_parallel <- parallel::mclapply(
    model_jobs,
    run_depth_model_job,
    mc.cores = model_workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  model_results_parallel <- lapply(model_jobs, run_depth_model_job)
}

depth_adjusted_list <- lapply(model_results_parallel, `[[`, "adjusted")
depth_support_list <- lapply(model_results_parallel, `[[`, "support")
global_test_list <- lapply(model_results_parallel, `[[`, "global_test")
pairwise_list <- lapply(model_results_parallel, `[[`, "pairwise")

depth_adjusted_df <- data.table::rbindlist(
  depth_adjusted_list,
  use.names = TRUE,
  fill = TRUE
)
depth_support_df <- data.table::rbindlist(
  depth_support_list,
  use.names = TRUE,
  fill = TRUE
)
global_test_df <- data.table::rbindlist(
  global_test_list,
  use.names = TRUE,
  fill = TRUE
)
pairwise_df <- data.table::rbindlist(
  pairwise_list,
  use.names = TRUE,
  fill = TRUE
)

p_col <- intersect(c("Pr(Chi)", "Pr(>Chi)", "Pr(>Chisq)"), names(global_test_df))
global_test_df$Global_city_P_BH <- NA_real_
if (length(p_col) > 0) {
  comparison_rows <- !is.na(global_test_df[[p_col[1]]])
  global_test_df$Global_city_P_BH[comparison_rows] <- p.adjust(
    global_test_df[[p_col[1]]][comparison_rows],
    method = "BH"
  )
}

write_csv(depth_adjusted_df, "depth_adjusted_alpha_richness_selected_host_order.csv")
write_csv(depth_support_df, "sequencing_depth_common_support_by_model.csv")
write_csv(global_test_df, "depth_model_global_city_tests.csv")
write_csv(pairwise_df, "depth_model_pairwise_city_contrasts_tukey.csv")

for (level in c("Species", "Genus", "Family")) {
  plot_dat <- depth_adjusted_df[depth_adjusted_df$Level == level, , drop = FALSE]
  city_order <- plot_dat$City[
    plot_dat$Set == "All"
  ][order(plot_dat$Adjusted_alpha_richness[plot_dat$Set == "All"])]
  plot_dat$City <- factor(plot_dat$City, levels = city_order)

  p <- ggplot(
    plot_dat,
    aes(x = City, y = Adjusted_alpha_richness, color = Set)
  ) +
    geom_errorbar(
      aes(ymin = CI95_lower, ymax = CI95_upper),
      width = 0.18,
      position = position_dodge(width = 0.5)
    ) +
    geom_point(
      size = 2,
      position = position_dodge(width = 0.5)
    ) +
    coord_flip() +
    scale_color_manual(
      values = c(All = all_color, GE2 = ge2_color),
      labels = get_set_labels(level)
    ) +
    labs(
      x = NULL,
      y = "Depth-adjusted mean sample-level richness",
      color = NULL,
      title = paste0(level, "-level depth-adjusted alpha richness"),
      subtitle = paste0(
        "Restricted to ",
        selected_host,
        "; common depth chosen inside the 5%-95% depth interval of every included city."
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title.position = "plot"
    )

  save_plot(
    p,
    paste0(tolower(level), "_depth_adjusted_alpha_selected_host_order"),
    9,
    9
  )
}


if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(
    list(
      Matching_summary = matching_summary,
      City_sample_counts = city_counts,
      Host_support = host_candidates,
      Host_selection = host_selection_summary,
      Sample_size_standardized = as.data.frame(size_est_df),
      Coverage_standardized = as.data.frame(coverage_est_df),
      Coverage_validation = as.data.frame(coverage_validation_df),
      Depth_adjusted_alpha = as.data.frame(depth_adjusted_df),
      Depth_support = depth_interval,
      Depth_support_by_model = as.data.frame(depth_support_df),
      Global_city_tests = as.data.frame(global_test_df),
      Pairwise_city_contrasts = as.data.frame(pairwise_df)
    ),
    file = file.path(
      out_dir,
      "city_rarefaction_standardization_fast_parallel.xlsx"
    ),
    overwrite = TRUE
  )
}

attempted_coverage <- coverage_validation_df$Coverage_standardization_available

if (
  any(
    !coverage_validation_df$Coverage_validation_pass[attempted_coverage],
    na.rm = TRUE
  )
) {
  stop("At least one attempted coverage standardization failed validation.")
}
if (
  any(
    !coverage_validation_df$Two_times_limit_pass[attempted_coverage],
    na.rm = TRUE
  )
) {
  stop("At least one attempted coverage estimate exceeded the two-times limit.")
}
if (!all(depth_support_df$Common_depth_in_central_interval)) {
  stop("At least one depth-adjusted prediction is outside city depth support.")
}
