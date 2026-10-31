library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(forcats)
library(vegan)


meta <- read_excel("metadata2.xlsx") %>%
  rename(SampleID = Library) %>%
  mutate(
    SampleID = trimws(as.character(SampleID)),
    SampleID = gsub("Birdfe", "BirdFe", SampleID),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  filter(
    !is.na(SampleID),
    SampleID != "",
    !is.na(lon),
    !is.na(lat)
  ) %>%
  mutate(
    PointID = paste0(
      sprintf("%.5f", lon),
      "_",
      sprintf("%.5f", lat)
    )
  )

abu <- read.table(
  "RPKM_veter_percent_matrix_clean.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

colnames(abu)[1] <- "vOTU"
colnames(abu) <- gsub("Birdfe", "BirdFe", colnames(abu))

shared_samples <- intersect(
  colnames(abu)[-1],
  meta$SampleID
)

if (length(shared_samples) < 3) {
  stop("Too few matched samples between metadata and the RPKM matrix.")
}

abu2 <- abu[
  ,
  c("vOTU", shared_samples),
  drop = FALSE
]

meta2 <- meta %>%
  filter(SampleID %in% shared_samples)

abu_long <- abu2 %>%
  pivot_longer(
    cols = -vOTU,
    names_to = "SampleID",
    values_to = "Abu"
  ) %>%
  inner_join(
    meta2 %>% select(SampleID, PointID),
    by = "SampleID"
  ) %>%
  mutate(
    PA = as.integer(Abu > 0)
  )

point_votu <- abu_long %>%
  group_by(PointID, vOTU) %>%
  summarise(
    PA = max(PA),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = vOTU,
    values_from = PA,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(point_votu) <- point_votu$PointID
point_votu$PointID <- NULL

X_pa <- as.matrix(point_votu)


dist_obs <- vegdist(
  X_pa,
  method = "jaccard",
  binary = TRUE
)

sim_mat <- 1 - as.matrix(dist_obs)

K <- 4

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
  slice_head(n = K) %>%
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


set.seed(1)

B <- 9999
upper_idx <- which(upper.tri(sim_mat))

null_vals <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(upper_idx)
)

for (b in seq_len(B)) {

  Xb <- apply(
    X_pa,
    2,
    sample
  )

  sim_b <- 1 - as.matrix(
    vegdist(
      Xb,
      method = "jaccard",
      binary = TRUE
    )
  )

  null_vals[b, ] <- sim_b[upper_idx]

  if (b %% 1000 == 0) {
  }
}

thr95 <- as.numeric(
  quantile(
    null_vals,
    0.95,
    na.rm = TRUE
  )
)

thr99 <- as.numeric(
  quantile(
    null_vals,
    0.99,
    na.rm = TRUE
  )
)


get_shared <- function(a, b, X) {
  sum(
    X[a, ] == 1 &
      X[b, ] == 1
  )
}

edges0$shared <- mapply(
  get_shared,
  edges0$p1,
  edges0$p2,
  MoreArgs = list(X = X_pa)
)

min_shared <- 30
min_sim <- 0.06

edges_sig <- edges0 %>%
  filter(
    sim >= thr99,
    shared >= min_shared,
    sim >= min_sim
  )


if (nrow(edges_sig) == 0) {
  stop("No edges passed the specified thresholds.")
}

write.table(
  edges_sig,
  "significant_edges_9999.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    permutations = B,
    threshold_95 = thr95,
    threshold_99 = thr99,
    min_shared = min_shared,
    min_similarity = min_sim,
    top_k = K
  ),
  "network_thresholds_9999.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

rm(
  dist_obs,
  sim_df,
  null_vals
)

invisible(gc())


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

TOP_N_PARTNER <- 12
TOP_K_PARTNER <- 6

top_sites_tbl <- edge_long %>%
  group_by(site) %>%
  summarise(
    total_shared = sum(shared, na.rm = TRUE),
    degree = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(total_shared)) %>%
  slice_head(n = TOP_N_PARTNER)

top_sites <- top_sites_tbl$site

write.table(
  top_sites_tbl,
  "Top_sites_total_shared_clean.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_df <- edge_long %>%
  filter(site %in% top_sites) %>%
  group_by(site) %>%
  arrange(desc(shared), partner) %>%
  mutate(
    partner2 = ifelse(
      row_number() <= TOP_K_PARTNER,
      partner,
      "Others"
    )
  ) %>%
  ungroup() %>%
  group_by(site, partner2) %>%
  summarise(
    shared = sum(shared, na.rm = TRUE),
    .groups = "drop"
  )

site_order <- plot_df %>%
  group_by(site) %>%
  summarise(
    total = sum(shared),
    .groups = "drop"
  ) %>%
  arrange(total) %>%
  pull(site)

plot_df <- plot_df %>%
  mutate(
    site = factor(
      site,
      levels = site_order
    )
  )

p_stack <- ggplot(
  plot_df,
  aes(
    x = site,
    y = shared,
    fill = partner2
  )
) +
  geom_col(width = 0.78) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Shared vOTUs (counts)",
    fill = "Shared with (partner site)",
    title = paste0(
      "Shared vOTUs composition across top ",
      TOP_N_PARTNER,
      " sites"
    ),
    subtitle = paste0(
      "Per site: top ",
      TOP_K_PARTNER,
      " partner sites shown; remaining partners merged as Others"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_stack)

ggsave(
  "Top_sites_shared_vOTUs_stacked_by_partner_clean.pdf",
  p_stack,
  width = 10,
  height = 6.8
)

write.table(
  plot_df,
  "Top_sites_shared_vOTUs_stacked_plotdata_clean.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(forcats)

# Input checks
stopifnot(exists("edges_sig"))
stopifnot(all(c("p1","p2","shared") %in% colnames(edges_sig)))
stopifnot(exists("X_pa"))
stopifnot(!is.null(rownames(X_pa)))
stopifnot(!is.null(colnames(X_pa)))

TOP_N <- 15      # number of sites to display, ranked by total_shared
TOP_K <- 99999   # use all significant neighbors by default; keep this value large to avoid limiting neighbors

# Definition of unique vOTUs:
# TRUE  = globally unique (observed at only one site across all sites)
# FALSE = relatively unique (present minus shared; not shared with significant neighbors)
UNIQUE_GLOBAL <- TRUE

# To avoid overcrowding, retain the most abundant families globally and merge the remainder as Others
TOP_FAM_GLOBAL <- 18

tax <- read_excel("veter_votus_clean.xlsx") %>%
  rename(vOTU = votus) %>%
  mutate(
    vOTU = trimws(as.character(vOTU)),
    Family = trimws(as.character(Family))
  ) %>%
  filter(!is.na(vOTU), vOTU != "", !is.na(Family), Family != "")

tax_map <- setNames(tax$Family, tax$vOTU)
votu_names <- colnames(X_pa)

edge_long <- bind_rows(
  edges_sig %>% transmute(site = as.character(p1), partner = as.character(p2), shared = as.numeric(shared)),
  edges_sig %>% transmute(site = as.character(p2), partner = as.character(p1), shared = as.numeric(shared))
)

top_sites <- edge_long %>%
  group_by(site) %>%
  summarise(total_shared = sum(shared, na.rm = TRUE), degree = n(), .groups = "drop") %>%
  arrange(desc(total_shared)) %>%
  slice_head(n = TOP_N) %>%
  pull(site)

present_idx <- lapply(rownames(X_pa), function(s) which(X_pa[s, ] == 1))
names(present_idx) <- rownames(X_pa)

get_shared_votus_vec <- function(a, b) {
  ia <- present_idx[[a]]
  ib <- present_idx[[b]]
  if (length(ia) == 0 || length(ib) == 0) return(character(0))
  votu_names[intersect(ia, ib)]
}

site_shared <- vector("list", length(top_sites))
names(site_shared) <- top_sites

for (s in top_sites) {
  nb_tbl <- edge_long %>%
    filter(site == s) %>%
    arrange(desc(shared)) %>%
    slice_head(n = TOP_K)
  
  nb <- nb_tbl$partner
  if (length(nb) == 0) {
    site_shared[[s]] <- character(0)
  } else {
    site_shared[[s]] <- unique(unlist(lapply(nb, function(x) get_shared_votus_vec(s, x)), use.names = FALSE))
  }
}

site_unique <- vector("list", length(top_sites))
names(site_unique) <- top_sites

if (UNIQUE_GLOBAL) {
  # Globally unique: the vOTU occurs at only one site across all sites
  col_occ <- colSums(X_pa == 1, na.rm = TRUE)  # number of sites in which each vOTU occurs
  for (s in top_sites) {
    idx_s <- present_idx[[s]]
    if (length(idx_s) == 0) {
      site_unique[[s]] <- character(0)
    } else {
      v_s <- votu_names[idx_s]
      site_unique[[s]] <- v_s[col_occ[v_s] == 1]
    }
  }
} else {
  # Relatively unique: present minus shared (not shared with significant neighbors)
  for (s in top_sites) {
    idx_s <- present_idx[[s]]
    v_s <- if (length(idx_s) == 0) character(0) else votu_names[idx_s]
    site_unique[[s]] <- setdiff(v_s, site_shared[[s]])
  }
}

count_family <- function(vv) {
  if (length(vv) == 0) return(tibble(Family = "None", n_vOTU = 0))
  fam <- unname(tax_map[vv])
  fam <- fam[!is.na(fam)]
  if (length(fam) == 0) return(tibble(Family = "Unclassified", n_vOTU = length(vv)))
  tibble(Family = fam) %>% count(Family, name = "n_vOTU")
}

df_shared <- bind_rows(lapply(names(site_shared), function(s) {
  count_family(site_shared[[s]]) %>% mutate(site = s, Type = "Shared")
}))

df_unique <- bind_rows(lapply(names(site_unique), function(s) {
  count_family(site_unique[[s]]) %>% mutate(site = s, Type = "Unique")
}))

df_fam <- bind_rows(df_shared, df_unique)

write.table(df_fam, "Site_shared_unique_vOTU_family_counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

top_fams <- df_fam %>%
  filter(Family != "None") %>%
  group_by(Family) %>%
  summarise(total = sum(n_vOTU), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = TOP_FAM_GLOBAL) %>%
  pull(Family)

plot_df <- df_fam %>%
  mutate(Family2 = ifelse(Family %in% top_fams, Family, "Others")) %>%
  group_by(site, Type, Family2) %>%
  summarise(n_vOTU = sum(n_vOTU), .groups = "drop")

# Order families within each panel by n_vOTU
plot_df <- plot_df %>%
  group_by(site, Type) %>%
  mutate(Family2 = fct_reorder(Family2, n_vOTU, .desc = TRUE)) %>%
  ungroup()

# Order sites by total shared vOTUs
site_order <- plot_df %>%
  group_by(site, Type) %>%
  summarise(total = sum(n_vOTU), .groups = "drop") %>%
  filter(Type == "Shared") %>%
  arrange(desc(total)) %>%
  pull(site)

plot_df$site <- factor(plot_df$site, levels = rev(site_order))  # reverse order for facet display

p <- ggplot(plot_df, aes(x = Family2, y = n_vOTU)) +
  geom_col(width = 0.82) +
  facet_wrap(vars(site, Type), ncol = 5, scales = "free_y") +
  labs(
    x = NULL,
    y = "Number of vOTUs",
    title = "Family composition of shared vs unique vOTUs per site",
    subtitle = paste0(
      "Shared = union of vOTUs shared with significant neighbors (edges_sig); ",
      ifelse(UNIQUE_GLOBAL,
             "Unique = vOTUs observed only at this site across all sites",
             "Unique = present vOTUs minus shared vOTUs (relative to significant neighbors)"),
      "; Top ", TOP_FAM_GLOBAL, " families shown, others merged as 'Others'."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    strip.text = element_text(face = "bold", size = 9)
  )

print(p)
ggsave("Site_shared_vs_unique_vOTU_family_facets.pdf", p, width = 14, height = 8)
