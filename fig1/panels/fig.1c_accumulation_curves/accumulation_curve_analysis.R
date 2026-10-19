# Accumulation curve analysis and stability assessment

library(readr)
library(stringr)
library(readxl)
library(vegan)
library(openxlsx)
library(parallel)

seeds <- c(123, 456, 789, 2025, 2026)
n_randomizations <- 1000L
n_cores <- if (.Platform$OS.type == "windows") 1L else 2L
plot_seed <- 123

output_xlsx <- "accumulation_curve_stability_ge2.xlsx"

run_accumulation <- function(comm, seeds, permutations = 1000L, n_cores = 2L) {
  runs <- mclapply(
    seeds,
    function(s) {
      set.seed(s)
      specaccum(
        comm,
        method = "random",
        permutations = permutations
      )
    },
    mc.cores = min(n_cores, length(seeds)),
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )

  failed <- vapply(runs, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("Accumulation analysis failed for seed(s): ",
         paste(seeds[failed], collapse = ", "))
  }

  richness <- sapply(runs, function(x) x$richness)
  if (is.null(dim(richness))) {
    richness <- matrix(richness, ncol = length(seeds))
  }

  colnames(richness) <- paste0("Seed_", seeds)
  n_samples <- runs[[1]]$sites

  mean_richness <- rowMeans(richness)
  sd_richness <- apply(richness, 1, sd)
  minimum <- apply(richness, 1, min)
  maximum <- apply(richness, 1, max)
  range_value <- maximum - minimum
  cv <- sd_richness / pmax(mean_richness, 1) * 100
  relative_range <- range_value / pmax(mean_richness, 1) * 100

  detail <- data.frame(
    Number_of_samples = n_samples,
    richness,
    Mean_richness = mean_richness,
    SD = sd_richness,
    CV_percent = cv,
    Minimum = minimum,
    Maximum = maximum,
    Range = range_value,
    Relative_range_percent = relative_range,
    check.names = FALSE
  )

  idx <- round(c(0.25, 0.50, 0.75, 1.00) * length(n_samples))
  idx <- pmax(idx, 1)

  check <- data.frame(
    Sampling_fraction = c("25%", "50%", "75%", "100%"),
    Number_of_samples = n_samples[idx],
    richness[idx, , drop = FALSE],
    Mean_richness = mean_richness[idx],
    SD = sd_richness[idx],
    CV_percent = cv[idx],
    Minimum = minimum[idx],
    Maximum = maximum[idx],
    Range = range_value[idx],
    Relative_range_percent = relative_range[idx],
    check.names = FALSE
  )

  summary <- data.frame(
    Metric = c(
      "Number of independent random seeds",
      "Randomizations per seed",
      "Maximum absolute range among runs",
      "Maximum relative range among runs (%)",
      "Maximum CV among runs (%)",
      "Mean CV among runs (%)"
    ),
    Value = c(
      length(seeds),
      permutations,
      max(range_value),
      max(relative_range),
      max(cv),
      mean(cv)
    )
  )

  rm(runs, richness)
  invisible(gc())

  list(
    detail = detail,
    check = check,
    summary = summary
  )
}

add_results <- function(wb, prefix, result) {
  addWorksheet(wb, paste0(prefix, "_Detail"))
  writeData(wb, paste0(prefix, "_Detail"), result$detail)

  addWorksheet(wb, paste0(prefix, "_Check"))
  writeData(wb, paste0(prefix, "_Check"), result$check)

  addWorksheet(wb, paste0(prefix, "_Summary"))
  writeData(wb, paste0(prefix, "_Summary"), result$summary)
}

plot_accumulation <- function(all_result, ge2_result, file, main, ylab,
                              all_label, ge2_label) {
  seed_col <- paste0("Seed_", plot_seed)

  pdf(file, width = 6, height = 5)

  plot(
    all_result$detail$Number_of_samples,
    all_result$detail[[seed_col]],
    type = "l",
    col = "#7792c8",
    lwd = 12,
    xlab = "Number of samples",
    ylab = ylab,
    main = main
  )

  lines(
    ge2_result$detail$Number_of_samples,
    ge2_result$detail[[seed_col]],
    col = "#c16e54",
    lwd = 12
  )

  legend(
    "bottomright",
    legend = c(all_label, ge2_label),
    col = c("#7792c8", "#c16e54"),
    lwd = 3,
    bty = "n"
  )

  dev.off()
}

build_cluster_matrix <- function(cluster_file, species_pa) {
  x <- read_xlsx(
    cluster_file,
    col_names = FALSE,
    .name_repair = "minimal"
  )

  x <- as.matrix(x)
  x <- apply(x, 2, as.character)

  representatives <- x[, 1]
  cluster_size <- rowSums(!is.na(x) & x != "")

  cluster_pa <- matrix(
    0L,
    nrow = length(representatives),
    ncol = ncol(species_pa),
    dimnames = list(representatives, colnames(species_pa))
  )

  for (i in seq_len(nrow(x))) {
    members <- x[i, ]
    members <- members[!is.na(members) & members != ""]
    members <- intersect(members, rownames(species_pa))

    if (length(members) == 0L) {
      next
    }

    cluster_pa[i, colSums(species_pa[members, , drop = FALSE]) > 0] <- 1L
  }

  list(
    all = cluster_pa,
    ge2 = cluster_pa[cluster_size >= 2L, , drop = FALSE]
  )
}


abu <- read_tsv(
  "RPKM_allvirus_matrix_clean.tsv",
  show_col_types = FALSE,
  progress = FALSE
)

contigs <- abu$Contig
abu$Contig <- NULL

abu_mat <- as.matrix(abu)
rownames(abu_mat) <- contigs

pa_species <- (abu_mat > 0) * 1L

rm(abu, abu_mat, contigs)
invisible(gc())

wb <- createWorkbook()


# Species level


species_clusters <- read_tsv(
  "species_clusters_clean.tsv",
  col_names = c("vOTU", "members"),
  col_types = cols(
    vOTU = col_character(),
    members = col_character()
  ),
  show_col_types = FALSE,
  progress = FALSE
)

species_cluster_size <- ifelse(
  is.na(species_clusters$members) | species_clusters$members == "",
  0L,
  str_count(species_clusters$members, fixed(",")) + 1L
)

species_ge2_ids <- species_clusters$vOTU[species_cluster_size >= 2L]
species_ge2_ids <- intersect(species_ge2_ids, rownames(pa_species))
pa_species_ge2 <- pa_species[species_ge2_ids, , drop = FALSE]

rm(species_clusters, species_cluster_size, species_ge2_ids)
invisible(gc())

species_all <- run_accumulation(
  t(pa_species),
  seeds,
  n_randomizations,
  n_cores
)

species_ge2 <- run_accumulation(
  t(pa_species_ge2),
  seeds,
  n_randomizations,
  n_cores
)

add_results(wb, "Species_All", species_all)
add_results(wb, "Species_ge2", species_ge2)

plot_accumulation(
  species_all,
  species_ge2,
  "species_accumulation_curves_ge2_clean.pdf",
  "Species-level accumulation curves",
  "Accumulated vOTUs",
  "All vOTUs",
  "vOTUs >=2 sequences"
)

saveWorkbook(wb, output_xlsx, overwrite = TRUE)

rm(pa_species_ge2, species_all, species_ge2)
invisible(gc())


# Genus level


genus <- build_cluster_matrix(
  "genus_clusters_clean.xlsx",
  pa_species
)

genus_all <- run_accumulation(
  t(genus$all),
  seeds,
  n_randomizations,
  n_cores
)

genus_ge2 <- run_accumulation(
  t(genus$ge2),
  seeds,
  n_randomizations,
  n_cores
)

add_results(wb, "Genus_All", genus_all)
add_results(wb, "Genus_ge2", genus_ge2)

plot_accumulation(
  genus_all,
  genus_ge2,
  "genus_accumulation_curves_ge2_clean.pdf",
  "Genus-level accumulation curves",
  "Accumulated genus vOTUs",
  "All genus vOTUs",
  "Genus vOTUs >=2 sequences"
)

saveWorkbook(wb, output_xlsx, overwrite = TRUE)

rm(genus, genus_all, genus_ge2)
invisible(gc())


# Family level


family <- build_cluster_matrix(
  "family_clusters_clean.xlsx",
  pa_species
)

rm(pa_species)
invisible(gc())

family_all <- run_accumulation(
  t(family$all),
  seeds,
  n_randomizations,
  n_cores
)

family_ge2 <- run_accumulation(
  t(family$ge2),
  seeds,
  n_randomizations,
  n_cores
)

add_results(wb, "Family_All", family_all)
add_results(wb, "Family_ge2", family_ge2)

plot_accumulation(
  family_all,
  family_ge2,
  "family_accumulation_curves_ge2_clean.pdf",
  "Family-level accumulation curves",
  "Accumulated family vOTUs",
  "All family vOTUs",
  "Family vOTUs >=2 sequences"
)

saveWorkbook(wb, output_xlsx, overwrite = TRUE)

rm(family, family_all, family_ge2)
invisible(gc())
