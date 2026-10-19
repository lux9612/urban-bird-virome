library(tidyverse)
library(readxl)
library(vegan)

# Parameters

N_PERM <- 9999
K_NEIGHBOR <- 4
MIN_SHARED <- 30
MIN_SIM <- 0.06

TOP_N_PARTNER <- 12
TOP_K_PARTNER <- 6

TOP_N_FAMILY <- 15
TOP_FAM_GLOBAL <- 18
UNIQUE_GLOBAL <- TRUE

set.seed(1)


meta <- read_excel("metadata2.xlsx") %>%
  rename(SampleID = Library) %>%
  mutate(
    SampleID = trimws(as.character(SampleID)),
    SampleID = gsub("Birdfe", "BirdFe", SampleID),
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    PointID = paste0(sprintf("%.5f", lon), "_", sprintf("%.5f", lat))
  ) %>%
  filter(
    !is.na(SampleID),
    SampleID != "",
    !is.na(lon),
    !is.na(lat)
  )

abu <- read.table(
  "RPKM_veter_percent_matrix_clean.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

colnames(abu)[1] <- "vOTU"
colnames(abu) <- gsub("Birdfe", "BirdFe", colnames(abu))

shared_samples <- intersect(colnames(abu)[-1], meta$SampleID)

if (length(shared_samples) < 3) {
  stop("Too few matched samples. Check SampleID in metadata2.xlsx and the abundance matrix.")
}

abu2 <- abu[, c("vOTU", shared_samples), drop = FALSE]
meta2 <- meta %>%
  filter(SampleID %in% shared_samples)


abu_long <- abu2 %>%
  pivot_longer(
    cols = -vOTU,
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  inner_join(
    meta2 %>% select(SampleID, PointID),
    by = "SampleID"
  ) %>%
  mutate(PA = as.integer(Abundance > 0))

point_votu <- abu_long %>%
  group_by(PointID, vOTU) %>%
  summarise(PA = max(PA), .groups = "drop") %>%
  pivot_wider(
    names_from = vOTU,
    values_from = PA,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(point_votu) <- point_votu$PointID
point_votu$PointID <- NULL

X_pa <- as.matrix(point_votu)


dist_obs <- vegdist(X_pa, method = "jaccard", binary = TRUE)
sim_mat <- 1 - as.matrix(dist_obs)

sim_df <- as.data.frame(
  as.table(sim_mat),
  stringsAsFactors = FALSE
) %>%
  transmute(
    p1 = as.character(Var1),
    p2 = as.character(Var2),
    sim = as.numeric(Freq)
  ) %>%
  filter(p1 != p2)

edges0 <- sim_df %>%
  group_by(p1) %>%
  arrange(desc(sim), p2) %>%
  slice_head(n = K_NEIGHBOR) %>%
  ungroup() %>%
  mutate(
    u = pmin(p1, p2),
    v = pmax(p1, p2)
  ) %>%
  distinct(u, v, .keep_all = TRUE) %>%
  transmute(
    p1 = u,
    p2 = v,
    sim = sim
  )


upper_idx <- which(upper.tri(sim_mat))
null_vals <- matrix(
  NA_real_,
  nrow = N_PERM,
  ncol = length(upper_idx)
)

for (b in seq_len(N_PERM)) {
  X_perm <- apply(X_pa, 2, sample)

  sim_perm <- 1 - as.matrix(
    vegdist(X_perm, method = "jaccard", binary = TRUE)
  )

  null_vals[b, ] <- sim_perm[upper_idx]

  if (b %% 1000 == 0) {
  }
}

thr95 <- as.numeric(
  quantile(null_vals, 0.95, na.rm = TRUE)
)

thr99 <- as.numeric(
  quantile(null_vals, 0.99, na.rm = TRUE)
)


get_shared <- function(a, b, X) {
  sum(X[a, ] == 1 & X[b, ] == 1)
}

edges0$shared <- mapply(
  get_shared,
  edges0$p1,
  edges0$p2,
  MoreArgs = list(X = X_pa)
)

edges_sig <- edges0 %>%
  filter(
    sim >= thr99,
    shared >= MIN_SHARED,
    sim >= MIN_SIM
  ) %>%
  arrange(desc(shared), desc(sim))

if (nrow(edges_sig) == 0) {
  stop("No edges passed the thresholds.")
}

write.table(
  edges_sig,
  "Significant_network_edges_9999.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    permutations = N_PERM,
    threshold_95 = thr95,
    threshold_99 = thr99,
    min_shared = MIN_SHARED,
    min_similarity = MIN_SIM,
    top_k_neighbor = K_NEIGHBOR
  ),
  "Network_null_thresholds_9999.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


edge_long <- bind_rows(
  edges_sig %>%
    transmute(
      site = as.character(p1),
      partner = as.character(p2),
      shared = as.numeric(shared)
    ),
  edges_sig %>%
    transmute(
      site = as.character(p2),
      partner = as.character(p1),
      shared = as.numeric(shared)
    )
)

top_sites_tbl <- edge_long %>%
  group_by(site) %>%
  summarise(
    total_shared = sum(shared, na.rm = TRUE),
    degree = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(total_shared)) %>%
  slice_head(n = TOP_N_PARTNER)

top_sites_partner <- top_sites_tbl$site

partner_plot_df <- edge_long %>%
  filter(site %in% top_sites_partner) %>%
  group_by(site) %>%
  arrange(desc(shared), partner) %>%
  mutate(
    partner_plot = if_else(
      row_number() <= TOP_K_PARTNER,
      partner,
      "Others"
    )
  ) %>%
  ungroup() %>%
  group_by(site, partner_plot) %>%
  summarise(
    shared = sum(shared, na.rm = TRUE),
    .groups = "drop"
  )

site_order_partner <- partner_plot_df %>%
  group_by(site) %>%
  summarise(
    total_shared = sum(shared),
    .groups = "drop"
  ) %>%
  arrange(total_shared) %>%
  pull(site)

partner_plot_df <- partner_plot_df %>%
  mutate(
    site = factor(site, levels = site_order_partner)
  )

p_partner <- ggplot(
  partner_plot_df,
  aes(x = site, y = shared, fill = partner_plot)
) +
  geom_col(width = 0.78) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Shared vOTUs",
    fill = "Partner site"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

print(p_partner)

ggsave(
  "Top_sites_shared_vOTUs_by_partner.png",
  p_partner,
  width = 10,
  height = 6.8,
  dpi = 600
)

write.table(
  top_sites_tbl,
  "Top_sites_total_shared.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  partner_plot_df,
  "Top_sites_shared_vOTUs_by_partner_plotdata.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
