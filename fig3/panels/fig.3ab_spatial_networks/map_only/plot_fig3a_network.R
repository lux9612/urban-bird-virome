# Saved seed-1 network; no permutation analysis.
# options(fig3a.data_dir = "path/to/data", fig3a.output_dir = "path/to/output")
required <- c("sf", "ggplot2", "dplyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
library(sf)
library(ggplot2)
library(dplyr)

script_dirs <- character()
for (fr in sys.frames()) {
  if (!is.null(fr$ofile)) {
    script_dirs <- c(script_dirs, dirname(normalizePath(fr$ofile, mustWork = TRUE)))
  }
}
farg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_dirs) && length(farg)) {
  script_dirs <- dirname(normalizePath(sub("^--file=", "", farg[1]), mustWork = TRUE))
}
input_names <- c("fig3a_metadata_matched.csv", "fig3a_edges.csv", "china_tianditu.geojson")
data_dir <- getOption("fig3a.data_dir")
if (is.null(data_dir)) {
  candidates <- c(file.path(rev(script_dirs), "data"), "data",
                  "fig3/panels/fig.3ab_spatial_networks/map_only/data", "map_only/data")
  found <- candidates[vapply(candidates, function(p) all(file.exists(file.path(p, input_names))), logical(1))]
  if (!length(found)) stop("Set options(fig3a.data_dir = ...) to the directory containing: ", paste(input_names, collapse = ", "))
  data_dir <- found[1]
}
absent <- input_names[!file.exists(file.path(data_dir, input_names))]
if (length(absent)) stop("Missing input files: ", paste(absent, collapse = ", "))
output_dir <- getOption("fig3a.output_dir", if (length(script_dirs)) tail(script_dirs, 1) else getwd())
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
meta_used <- read.csv(file.path(data_dir, "fig3a_metadata_matched.csv"), stringsAsFactors = FALSE)
significant_edges <- read.csv(file.path(data_dir, "fig3a_edges.csv"), stringsAsFactors = FALSE)
check_columns <- function(x, columns, label) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
  if (!nrow(x) || anyNA(x[, columns, drop = FALSE])) stop(label, " is empty or contains missing values")
}
check_columns(meta_used, c("SampleID", "PointID", "lon", "lat"), "Metadata")
check_columns(significant_edges, c("p1", "p2", "similarity", "shared_vOTUs"), "Edges")
stopifnot(!anyDuplicated(meta_used$SampleID), all(nzchar(meta_used$SampleID)),
          all(nzchar(meta_used$PointID)),
          is.numeric(meta_used$lon), is.numeric(meta_used$lat),
          all(is.finite(meta_used$lon)), all(is.finite(meta_used$lat)),
          all(abs(meta_used$lon) <= 180), all(abs(meta_used$lat) < 90),
          is.numeric(significant_edges$similarity), is.numeric(significant_edges$shared_vOTUs),
          all(is.finite(significant_edges$similarity)),
          all(significant_edges$similarity >= 0 & significant_edges$similarity <= 1),
          all(is.finite(significant_edges$shared_vOTUs)),
          all(significant_edges$shared_vOTUs > 0),
          all(significant_edges$shared_vOTUs == floor(significant_edges$shared_vOTUs)),
          all(significant_edges$p1 != significant_edges$p2))
edge_keys <- paste(pmin(significant_edges$p1, significant_edges$p2),
                   pmax(significant_edges$p1, significant_edges$p2), sep = "__")
stopifnot(!anyDuplicated(edge_keys))
coordinates_per_site <- meta_used %>% distinct(PointID, lon, lat) %>% count(PointID)
stopifnot(all(coordinates_per_site$n == 1))
crs_target <- 3857
basemap <- st_read(file.path(data_dir, "china_tianditu.geojson"), quiet = TRUE)
stopifnot(!is.na(st_crs(basemap)))
basemap <- st_transform(basemap, crs_target)
gtype <- as.character(st_geometry_type(basemap))
province_3857 <- basemap[gtype %in% c("POLYGON", "MULTIPOLYGON"), ]
boundary_lines_3857 <- basemap[gtype %in% c("LINESTRING", "MULTILINESTRING"), ]
stopifnot(nrow(province_3857) == 34, nrow(boundary_lines_3857) > 0)
map_bbox <- st_bbox(basemap)
output_png <- file.path(output_dir, "Fig3a_network.png")
output_pdf <- file.path(output_dir, "Fig3a_network.pdf")
# Network map
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

stopifnot(all(significant_edges$p1 %in% points$PointID),
          all(significant_edges$p2 %in% points$PointID))
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
  geom_sf(data = boundary_lines_3857, color = "black", linewidth = 0.4) +
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
