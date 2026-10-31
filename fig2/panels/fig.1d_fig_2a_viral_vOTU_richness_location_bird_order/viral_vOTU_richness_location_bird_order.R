# Viral vOTU richness across locations and bird orders

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Input files

abundance_file <- "RPKM_allvirus_percent_matrix_clean.tsv"
location_file <- "location.xlsx"
eukaryotic_votu_file <- "veter_otus_clean.xlsx"
bird_file <- "bird_species.xlsx"


abu <- read.table(
  abundance_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

loc <- read_xlsx(location_file)
euk_ids <- read_xlsx(eukaryotic_votu_file)
bird <- read_xlsx(bird_file)

# Standardize column names
names(loc)[names(loc) == "ID"] <- "Library"
names(euk_ids)[1] <- "Contig"

# Assign missing bird orders to "Unknown"
bird$Order[is.na(bird$Order)] <- "Unknown"


abu_long <- abu %>%
  pivot_longer(
    cols = -Contig,
    names_to = "Library",
    values_to = "value"
  )


abu_meta <- abu_long %>%
  left_join(loc, by = "Library") %>%
  left_join(bird, by = "Library") %>%
  filter(!is.na(location))


abu_pa <- abu_meta %>%
  mutate(
    present = value > 0
  )


combo_votu <- abu_pa %>%
  group_by(location, Order, Contig) %>%
  summarise(
    any_present = any(present),
    .groups = "drop"
  )


votu_rich_combo <- combo_votu %>%
  group_by(location, Order) %>%
  summarise(
    votu_richness = sum(any_present),
    .groups = "drop"
  )


sample_n_combo <- abu_meta %>%
  group_by(location, Order) %>%
  summarise(
    n_samples = n_distinct(Library),
    .groups = "drop"
  )

# Combine vOTU richness and sample size
div_combo <- votu_rich_combo %>%
  left_join(
    sample_n_combo,
    by = c("location", "Order")
  )


p_bird_order <- ggplot(
  div_combo,
  aes(
    x = location,
    y = Order
  )
) +
  geom_point(
    aes(
      size = votu_richness,
      fill = n_samples
    ),
    shape = 21,
    color = "black",
    alpha = 0.8
  ) +
  scale_size_continuous(
    name = "vOTU richness"
  ) +
  scale_fill_gradient(
    name = "Sample size",
    low = "lightblue",
    high = "darkblue"
  ) +
  coord_flip() +
  labs(
    x = "Location",
    y = "Bird order",
    title = "Viral vOTU richness by location and bird order"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_bird_order
