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


abu_loc <- abu_long %>%
  left_join(loc, by = "Library") %>%
  filter(!is.na(location))


abu_pa_type <- abu_loc %>%
  mutate(present = value > 0) %>%
  left_join(
    euk_ids %>%
      mutate(Type = "Eukaryotic vOTUs"),
    by = "Contig"
  ) %>%
  mutate(
    Type = ifelse(
      is.na(Type),
      "Phage vOTUs",
      Type
    )
  )


votu_rich_type <- abu_pa_type %>%
  group_by(location, Type, Contig) %>%
  summarise(
    any_present = any(present),
    .groups = "drop"
  ) %>%
  group_by(location, Type) %>%
  summarise(
    votu_richness = sum(any_present),
    .groups = "drop"
  )

# Order locations according to total vOTU richness
votu_rich_type <- votu_rich_type %>%
  group_by(location) %>%
  mutate(
    total_votus = sum(votu_richness)
  ) %>%
  ungroup() %>%
  mutate(
    location = reorder(location, total_votus)
  )


p_location <- ggplot(
  votu_rich_type,
  aes(
    x = location,
    y = votu_richness,
    fill = Type
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Eukaryotic vOTUs" = "#afd6ec",
      "Phage vOTUs" = "#94c5aa"
    )
  ) +
  labs(
    x = "Location",
    y = "Number of detected viral vOTUs",
    fill = "vOTU type",
    title = "Viral vOTU richness across locations"
  ) +
  theme_classic(base_size = 14)

p_location
