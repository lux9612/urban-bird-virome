# Spatial vOTU similarity network

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(vegan)
library(openxlsx)

# Input files
meta_file <- "metadata2.xlsx"
abundance_file <- "RPKM_veter_percent_matrix_clean.tsv"

# Fig.3a basemap: Tianditu via GeoJSON.CN; source-reported GS(2024)0650
basemap_file <- "map_only/data/china_tianditu.geojson"
if (!file.exists(basemap_file)) {
  basemap_file <- "https://geojson.cn/api/tiandi/100000.json"
}

# Network settings
top_k <- 4
min_shared <- 30
min_similarity <- 0.06

# Permutation settings
previous_permutations <- 499
n_permutations <- 9999
main_seed <- 1
stability_seeds <- c(1, 123, 456, 2025, 2026)

# Output
output_xlsx <- "spatial_network_permutation_stability_9999.xlsx"
output_png <- "geographical_network_9999perm_seed1.png"
output_pdf <- "geographical_network_9999perm_seed1.pdf"


# Map data
crs_target <- 3857
basemap <- st_read(basemap_file, quiet = TRUE)
stopifnot(!is.na(st_crs(basemap)))
basemap <- st_transform(basemap, crs_target)
gtype <- as.character(st_geometry_type(basemap))
province_3857 <- basemap[gtype %in% c("POLYGON", "MULTIPOLYGON"), ]
official_lines_3857 <- basemap[gtype %in% c("LINESTRING", "MULTILINESTRING"), ]
stopifnot(nrow(province_3857) == 34, nrow(official_lines_3857) > 0)
map_bbox <- st_bbox(basemap)


# Metadata and abundance matrix
meta <- read_excel(meta_file) %>%
  rename(SampleID = Library) %>%
  mutate(
    SampleID = trimws(as.character(SampleID)),
    SampleID = gsub("Birdfe", "BirdFe", SampleID),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  filter(
    !is.na(SampleID), SampleID != "",
    !is.na(lon), !is.na(lat)
  ) %>%
  mutate(
    PointID = paste0(
      sprintf("%.5f", lon), "_", sprintf("%.5f", lat)
    )
  )

abundance <- read.table(
  abundance_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

colnames(abundance)[1] <- "vOTU"
colnames(abundance) <- gsub("Birdfe", "BirdFe", colnames(abundance))

shared_samples <- intersect(colnames(abundance)[-1], meta$SampleID)
if (length(shared_samples) < 3) {
  stop("Too few matched samples. Check SampleID names.")
}

abundance <- abundance[, c("vOTU", shared_samples), drop = FALSE]
meta_used <- meta %>% filter(SampleID %in% shared_samples)


# Site-level presence/absence matrix
site_votu <- abundance %>%
  pivot_longer(
    cols = -vOTU,
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  inner_join(
    meta_used %>% select(SampleID, PointID),
    by = "SampleID"
  ) %>%
  mutate(PA = as.integer(Abundance > 0)) %>%
  group_by(PointID, vOTU) %>%
  summarise(PA = max(PA), .groups = "drop") %>%
  pivot_wider(
    names_from = vOTU,
    values_from = PA,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(site_votu) <- site_votu$PointID
site_votu$PointID <- NULL

X_pa <- as.matrix(site_votu)
storage.mode(X_pa) <- "numeric"


# Observed network
similarity_matrix <- 1 - as.matrix(
  vegdist(X_pa, method = "jaccard", binary = TRUE)
)

similarity_table <- as.data.frame(
  as.table(similarity_matrix),
  stringsAsFactors = FALSE
) %>%
  transmute(
    p1 = as.character(Var1),
    p2 = as.character(Var2),
    similarity = as.numeric(Freq)
  ) %>%
  filter(p1 != p2)

candidate_edges <- similarity_table %>%
  group_by(p1) %>%
  arrange(desc(similarity), p2) %>%
  slice_head(n = top_k) %>%
  ungroup() %>%
  mutate(
    u = pmin(p1, p2),
    v = pmax(p1, p2)
  ) %>%
  distinct(u, v, .keep_all = TRUE) %>%
  transmute(p1 = u, p2 = v, similarity = similarity)

count_shared <- function(a, b, X) {
  sum(X[a, ] == 1 & X[b, ] == 1)
}

candidate_edges$shared_vOTUs <- mapply(
  count_shared,
  candidate_edges$p1,
  candidate_edges$p2,
  MoreArgs = list(X = X_pa)
)

candidate_edges <- candidate_edges %>%
  mutate(
    edge_id = paste(pmin(p1, p2), pmax(p1, p2), sep = "__")
  )


# Null model
upper_index <- which(upper.tri(similarity_matrix))

run_permutation <- function(seed, n_perm) {
  set.seed(seed)

  null_values <- matrix(
    NA_real_,
    nrow = n_perm,
    ncol = length(upper_index)
  )

  for (i in seq_len(n_perm)) {
    X_perm <- apply(X_pa, 2, sample)

    if (is.null(dim(X_perm))) {
      X_perm <- matrix(X_perm, nrow = nrow(X_pa))
    }

    rownames(X_perm) <- rownames(X_pa)
    colnames(X_perm) <- colnames(X_pa)

    perm_similarity <- 1 - as.matrix(
      vegdist(X_perm, method = "jaccard", binary = TRUE)
    )

    null_values[i, ] <- perm_similarity[upper_index]

    if (i %% 1000 == 0 || i == n_perm) {
    }
  }

  threshold_95 <- as.numeric(
    quantile(null_values, 0.95, na.rm = TRUE)
  )
  threshold_99 <- as.numeric(
    quantile(null_values, 0.99, na.rm = TRUE)
  )

  retained_edges <- candidate_edges %>%
    filter(
      similarity >= threshold_99,
      shared_vOTUs >= min_shared,
      similarity >= min_similarity
    )

  list(
    seed = seed,
    permutations = n_perm,
    threshold_95 = threshold_95,
    threshold_99 = threshold_99,
    retained_edges = retained_edges
  )
}

previous_result <- run_permutation(
  main_seed,
  previous_permutations
)

stability_results <- lapply(
  stability_seeds,
  function(seed) run_permutation(seed, n_permutations)
)
names(stability_results) <- paste0("seed_", stability_seeds)

main_result <- stability_results[[paste0("seed_", main_seed)]]
main_edge_ids <- sort(main_result$retained_edges$edge_id)

stability_summary <- do.call(
  rbind,
  lapply(stability_results, function(x) {
    ids <- sort(x$retained_edges$edge_id)
    overlap_n <- length(intersect(ids, main_edge_ids))
    union_n <- length(union(ids, main_edge_ids))

    data.frame(
      seed = x$seed,
      permutations = x$permutations,
      threshold_95 = x$threshold_95,
      threshold_99 = x$threshold_99,
      retained_edges = nrow(x$retained_edges),
      overlap_with_seed1 = overlap_n,
      edge_set_jaccard_vs_seed1 =
        ifelse(union_n == 0, NA_real_, overlap_n / union_n),
      exact_same_as_seed1 = identical(ids, main_edge_ids)
    )
  })
)

previous_ids <- sort(previous_result$retained_edges$edge_id)
current_ids <- sort(main_result$retained_edges$edge_id)

shared_edge_n <- length(intersect(previous_ids, current_ids))
union_edge_n <- length(union(previous_ids, current_ids))

previous_vs_current <- data.frame(
  comparison = "Seed 1: 499 vs 9999 permutations",
  previous_permutations = previous_permutations,
  current_permutations = n_permutations,
  previous_threshold_95 = previous_result$threshold_95,
  current_threshold_95 = main_result$threshold_95,
  previous_threshold_99 = previous_result$threshold_99,
  current_threshold_99 = main_result$threshold_99,
  previous_retained_edges = nrow(previous_result$retained_edges),
  current_retained_edges = nrow(main_result$retained_edges),
  shared_edges = shared_edge_n,
  union_edges = union_edge_n,
  edge_set_jaccard =
    ifelse(union_edge_n == 0, NA_real_, shared_edge_n / union_edge_n),
  exact_same_edge_set = identical(previous_ids, current_ids)
)

all_retained_ids <- sort(
  unique(
    unlist(
      lapply(
        stability_results,
        function(x) x$retained_edges$edge_id
      )
    )
  )
)

edge_stability <- candidate_edges %>%
  filter(edge_id %in% all_retained_ids) %>%
  select(edge_id, p1, p2, similarity, shared_vOTUs)

for (seed in stability_seeds) {
  x <- stability_results[[paste0("seed_", seed)]]
  edge_stability[[paste0("retained_seed_", seed)]] <-
    edge_stability$edge_id %in% x$retained_edges$edge_id
}

retained_columns <- paste0("retained_seed_", stability_seeds)

edge_stability$n_runs_retained <- rowSums(
  edge_stability[, retained_columns, drop = FALSE]
)
edge_stability$run_frequency <-
  edge_stability$n_runs_retained / length(stability_seeds)

edge_stability <- edge_stability %>%
  arrange(
    desc(n_runs_retained),
    desc(similarity),
    desc(shared_vOTUs)
  )

threshold_statistics <- data.frame(
  metric = c(
    "Mean 99th-percentile threshold",
    "SD of 99th-percentile threshold",
    "CV of 99th-percentile threshold (%)",
    "Minimum 99th-percentile threshold",
    "Maximum 99th-percentile threshold",
    "Range of 99th-percentile threshold",
    "Relative range of 99th-percentile threshold (%)",
    "Minimum retained edges",
    "Maximum retained edges"
  ),
  value = c(
    mean(stability_summary$threshold_99),
    sd(stability_summary$threshold_99),
    100 * sd(stability_summary$threshold_99) /
      mean(stability_summary$threshold_99),
    min(stability_summary$threshold_99),
    max(stability_summary$threshold_99),
    diff(range(stability_summary$threshold_99)),
    100 * diff(range(stability_summary$threshold_99)) /
      mean(stability_summary$threshold_99),
    min(stability_summary$retained_edges),
    max(stability_summary$retained_edges)
  )
)

parameters <- data.frame(
  item = c(
    "Similarity metric",
    "Candidate-edge rule",
    "Minimum shared vOTUs",
    "Minimum Jaccard similarity",
    "Previous permutations",
    "Current permutations per run",
    "Independent seeds",
    "Main figure seed",
    "Permutation structure",
    "Preserved property",
    "Randomized property",
    "Null-model interpretation"
  ),
  value = c(
    "Binary Jaccard similarity",
    paste0("Top-", top_k, " similarity edges per site before deduplication"),
    min_shared,
    min_similarity,
    previous_permutations,
    n_permutations,
    paste(stability_seeds, collapse = ", "),
    main_seed,
    "Each vOTU presence/absence column was independently shuffled across sampling sites",
    "Overall site occupancy frequency of each vOTU",
    "Identity of sampling sites at which each vOTU occurred",
    "Observed between-site similarity compared with a null distribution retaining vOTU occupancy while randomizing geographic allocation"
  )
)


# Export results
wb <- createWorkbook()

write_sheet <- function(sheet, x) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, x)
  setColWidths(wb, sheet, cols = seq_len(ncol(x)), widths = "auto")
}

write_sheet("Parameters", parameters)
write_sheet("Previous_vs_9999", previous_vs_current)
write_sheet("Stability_Summary", stability_summary)
write_sheet("Threshold_Stats", threshold_statistics)
write_sheet("Edge_Stability", edge_stability)
write_sheet("Main_Edges_Seed1", main_result$retained_edges)
write_sheet("Previous_Edges_499", previous_result$retained_edges)
write_sheet("Candidate_Edges", candidate_edges)

saveWorkbook(wb, output_xlsx, overwrite = TRUE)


# Network map
significant_edges <- main_result$retained_edges
if (nrow(significant_edges) == 0) {
  stop("No edges retained under the 9999-permutation threshold.")
}

similarity_breaks <- unique(
  as.numeric(
    quantile(
      significant_edges$similarity,
      probs = c(0, 0.33, 0.66, 1),
      na.rm = TRUE
    )
  )
)

if (length(similarity_breaks) < 4) {
  similarity_breaks <- unique(
    seq(
      min(significant_edges$similarity),
      max(significant_edges$similarity),
      length.out = 4
    )
  )
}

if (length(similarity_breaks) < 2) {
  epsilon <- 1e-8
  similarity_breaks <- c(
    min(significant_edges$similarity) - epsilon,
    max(significant_edges$similarity) + epsilon
  )
}

significant_edges <- significant_edges %>%
  mutate(
    similarity_group = cut(
      similarity,
      breaks = similarity_breaks,
      include.lowest = TRUE,
      dig.lab = 4
    )
  )

base_colors <- c("#F6E8C3", "#DFC27D", "#BF812D", "#8C510A")
n_groups <- nlevels(significant_edges$similarity_group)

network_colors <- base_colors[
  seq_len(min(n_groups, length(base_colors)))
]
names(network_colors) <- levels(
  significant_edges$similarity_group
)[seq_along(network_colors)]

points <- meta_used %>%
  group_by(PointID) %>%
  summarise(
    lon = mean(lon),
    lat = mean(lat),
    n_libraries = n(),
    .groups = "drop"
  )

points_sf <- st_as_sf(
  points,
  coords = c("lon", "lat"),
  crs = 4326
) %>%
  st_transform(crs_target)

point_coordinates <- st_coordinates(points_sf)
point_index <- setNames(seq_len(nrow(points)), points$PointID)

bezier_line <- function(p1, p2, strength = 0.18, n = 80) {
  v <- p2 - p1
  d <- sqrt(sum(v^2))

  if (d == 0) {
    return(st_linestring(rbind(p1, p2)))
  }

  perpendicular <- c(-v[2], v[1]) / d
  midpoint <- (p1 + p2) / 2
  control <- midpoint + perpendicular * (d * strength)
  t <- seq(0, 1, length.out = n)

  xy <- cbind(
    (1 - t)^2 * p1[1] +
      2 * (1 - t) * t * control[1] +
      t^2 * p2[1],
    (1 - t)^2 * p1[2] +
      2 * (1 - t) * t * control[2] +
      t^2 * p2[2]
  )

  st_linestring(xy)
}

make_edge <- function(a, b) {
  p1 <- point_coordinates[point_index[a], ]
  p2 <- point_coordinates[point_index[b], ]
  bezier_line(p1, p2)
}

edge_geometries <- mapply(
  make_edge,
  significant_edges$p1,
  significant_edges$p2,
  SIMPLIFY = FALSE
)

edges_sf <- st_sf(
  similarity = significant_edges$similarity,
  shared_vOTUs = significant_edges$shared_vOTUs,
  similarity_group = significant_edges$similarity_group,
  geometry = st_sfc(edge_geometries, crs = crs_target)
)

network_plot <- ggplot() +
  geom_sf(data = province_3857, fill = NA, color = "black", linewidth = 0.35) +
  geom_sf(data = official_lines_3857, color = "black", linewidth = 0.4) +
  geom_sf(
    data = edges_sf,
    aes(
      color = similarity_group,
      linewidth = shared_vOTUs
    ),
    alpha = 0.9,
    lineend = "round"
  ) +
  geom_sf(
    data = points_sf,
    aes(size = n_libraries),
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  scale_color_manual(
    values = network_colors,
    name = "Jaccard similarity (vOTU)"
  ) +
  scale_linewidth_continuous(
    name = "Shared vOTUs",
    range = c(0.6, 2.2)
  ) +
  scale_size_continuous(
    name = "Libraries per site",
    range = c(2.8, 6.2)
  ) +
  coord_sf(
    xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
    ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
    expand = TRUE
  ) +
  labs(caption = "Basemap: Tianditu via GeoJSON.CN; source-reported GS(2024)0650") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.caption = element_text(size = 7, hjust = 0)
  ) +
  ggtitle(
    "Geographical network of shared avian viromes (vOTU-based)"
  )

print(network_plot)

ggsave(
  output_png,
  plot = network_plot,
  width = 10,
  height = 8,
  dpi = 600
)

ggsave(
  output_pdf,
  plot = network_plot,
  width = 10,
  height = 8
)
