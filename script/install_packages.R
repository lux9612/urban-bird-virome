cran <- c("RColorBrewer", "circlize", "cowplot", "dplyr", "forcats", "ggplot2", "ggraph", "ggrepel", "igraph", "knitr", "openxlsx", "parallel", "patchwork", "pheatmap", "purrr", "readr", "readxl", "rmarkdown", "scales", "sf", "stringr", "tibble", "tidygraph", "tidyr", "tidyverse", "vegan")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install("ComplexHeatmap", ask=FALSE, update=FALSE)
}
