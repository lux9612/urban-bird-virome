# Vertebrate-associated viral family heatmap

library(readxl)
library(dplyr)
library(pheatmap)

# Input files
abundance_file <- "RPKM_veter_percent_matrix_clean.tsv"
taxonomy_file <- "taxonomy.xlsx"
location_file <- "location.xlsx"


vmat <- read.table(
  abundance_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

tax <- read_xlsx(taxonomy_file)
loc <- read_xlsx(location_file)

# Standardize the location column name
names(loc)[names(loc) == "Location"] <- "location"


vmat2 <- vmat %>%
  left_join(tax, by = c("Contig" = "ID")) %>%
  filter(!is.na(Family))


family_mat <- vmat2 %>%
  group_by(Family) %>%
  summarise(
    across(where(is.numeric), sum),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(family_mat) <- family_mat$Family

family_mat <- family_mat[
  , names(family_mat) != "Family",
  drop = FALSE
]


heatmat <- log10(family_mat + 1)


sample_ids <- colnames(heatmat)

ann_col <- data.frame(
  ID = sample_ids,
  stringsAsFactors = FALSE
) %>%
  left_join(loc, by = "ID")

annotation_col <- data.frame(
  Location = ann_col$location
)

rownames(annotation_col) <- ann_col$ID


location_order <- c(
  "Aral, Xinjiang",
  "Beijing",
  "Chengdu, Sichuan",
  "Chifeng, Inner Mongolia",
  "Chongqing",
  "Guangzhou, Guangdong",
  "Haikou, Hainan",
  "Harbin, Heilongjiang",
  "Yining, Xinjiang",
  "Lanzhou, Gansu",
  "Nanjing, Jiangsu",
  "Nyingchi, Tibet",
  "Xining, Qinghai",
  "Shanghai",
  "Tai'an, Shandong",
  "Ürümqi, Xinjiang",
  "Wuhan, Hubei",
  "Xi'an, Shaanxi",
  "Xianyang, Shaanxi",
  "Yinchuan, Ningxia",
  "Zhengzhou, Henan",
  "Zhenjiang, Jiangsu"
)

annotation_col$Location <- factor(
  annotation_col$Location,
  levels = location_order
)


ord <- order(
  annotation_col$Location,
  rownames(annotation_col)
)

heatmat_ord <- heatmat[
  , ord,
  drop = FALSE
]

annotation_col_ord <- annotation_col[
  ord,
  ,
  drop = FALSE
]


pheatmap(
  heatmat_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(
    c("white", "deepskyblue", "lightsalmon")
  )(100),
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 6,
  annotation_col = annotation_col_ord,
  main = "Vertebrate-associated viral families (log10 RPKM + 1)"
)
