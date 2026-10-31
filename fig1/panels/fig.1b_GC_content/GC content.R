# GC content and genome length distribution of viral sequences

library(ggplot2)
library(patchwork)

# Load input data
df <- read.table(
  "species_votus_length_gc_clean.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Colors for CheckV quality categories
quality_colors <- c(
  "Complete"       = "#ce9255",
  "High-quality"   = "#efd98d",
  "Medium-quality" = "#cfebc6",
  "Low-quality"    = "#afd6ec",
  "Not-determined" = "#7f7f7f"
)

# Set the order of CheckV quality categories
df$checkv_quality <- factor(
  df$checkv_quality,
  levels = c(
    "Not-determined",
    "Low-quality",
    "Medium-quality",
    "High-quality",
    "Complete"
  )
)

# Main scatter plot
p_main <- ggplot(
  df,
  aes(
    x = Length_bp_log10,
    y = GC_percent,
    color = checkv_quality
  )
) +
  geom_point(
    size = 1.5,
    alpha = 1,
    shape = 16
  ) +
  scale_color_manual(values = quality_colors) +
  labs(
    x = "Genome length (log10 bp)",
    y = "GC content (%)",
    color = "CheckV quality"
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  )

# Histogram of genome length
p_top <- ggplot(
  df,
  aes(
    x = Length_bp_log10,
    fill = checkv_quality
  )
) +
  geom_histogram(
    binwidth = 0.01,
    color = NA,
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(values = quality_colors) +
  labs(y = "Count") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 5, 0, 5)
  )

# Histogram of GC content
p_right <- ggplot(
  df,
  aes(
    x = GC_percent,
    fill = checkv_quality
  )
) +
  geom_histogram(
    binwidth = 0.1,
    color = NA,
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(values = quality_colors) +
  coord_flip() +
  labs(x = "Count") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 0)
  )

# Combine the plots
left <- p_top / p_main +
  plot_layout(heights = c(1, 4))

right <- plot_spacer() / p_right +
  plot_layout(heights = c(1, 4))

final_plot <- (left | right) +
  plot_layout(widths = c(4, 1))

final_plot
