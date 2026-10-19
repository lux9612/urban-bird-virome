required_packages <- c("readxl", "dplyr", "igraph", "tidygraph", "ggraph")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(readxl)
library(dplyr)
library(igraph)
library(tidygraph)
library(ggraph)

input_file <- "virus_host.xlsx"
output_dir <- "figure6_results"
if (!file.exists(input_file)) stop("Input file not found: ", input_file, call. = FALSE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
df <- read_excel(input_file)
colnames(df) <- trimws(colnames(df))
required_columns <- c(
  "qseqid", "pident", "qcovhsp", "Virus family", "Virus genus",
  "Host class", "Host genus", "Host Species"
)
missing_columns <- setdiff(required_columns, colnames(df))
if (length(missing_columns) > 0) {
  stop("Missing required columns in ", input_file, ": ",
       paste(missing_columns, collapse = ", "), call. = FALSE)
}

prepare_association_data <- function(data) {


  selected_columns <- c(
    "qseqid",
    "pident",
    "qcovhsp",
    "Virus genus",
    "Virus family",
    "Host Species",
    "Host genus",
    "Host class"
  )


  out <- data %>%
    select(all_of(selected_columns)) %>%
    rename(
      virus = qseqid,
      genus = `Virus genus`,
      family = `Virus family`,
      species = `Host Species`,
      host_genus = `Host genus`,
      class = `Host class`
    )


  out <- out %>%
    filter(
      pident >= 99,
      qcovhsp >= 80
    )


  out %>%
    filter(
      !is.na(species) |
        !is.na(host_genus)
    ) %>%
    mutate(
      host = case_when(
        !is.na(species) & species != "" ~ species,
        TRUE ~ host_genus
      )
    ) %>%
    filter(
      !is.na(host),
      host != ""
    ) %>%
    distinct(
      virus,
      host,
      .keep_all = TRUE
    )
}

plot_association_network <- function(net_unique, output_prefix) {

  virus_host_count <- net_unique %>%
    count(
      virus,
      name = "host_n"
    )

  multi_host_virus <- virus_host_count %>%
    filter(host_n >= 2) %>%
    pull(virus)

  edges <- net_unique %>%
    transmute(
      from = paste0("V_", virus),
      to = paste0("H_", host)
    )

  virus_nodes <- net_unique %>%
    distinct(
      virus,
      genus,
      family
    ) %>%
    mutate(
      tax = case_when(
        !is.na(genus) & genus != "" ~ genus,
        !is.na(family) & family != "" ~ family,
        TRUE ~ "Unclassified"
      ),
      label = paste0(
        virus,
        " (",
        tax,
        ")"
      ),
      name = paste0("V_", virus),
      type = "Virus",
      circulation = ifelse(
        virus %in% multi_host_virus,
        "High circulation",
        "Single circulation"
      )
    ) %>%
    select(
      name,
      type,
      label,
      circulation
    ) %>%
    distinct(
      name,
      .keep_all = TRUE
    )

  host_nodes <- net_unique %>%
    distinct(
      host,
      class
    ) %>%
    mutate(
      label = ifelse(
        class == "Mammalia",
        paste0("* ", host),
        host
      ),
      name = paste0("H_", host),
      type = "Host",
      circulation = "Host"
    ) %>%
    select(
      name,
      type,
      label,
      circulation
    ) %>%
    distinct(
      name,
      .keep_all = TRUE
    )

  nodes <- bind_rows(
    virus_nodes,
    host_nodes
  ) %>%
    distinct(
      name,
      .keep_all = TRUE
    )

  graph_object <- graph_from_data_frame(
    edges,
    vertices = nodes,
    directed = FALSE
  )

  tbl_graph_object <- as_tbl_graph(graph_object)

  p <- ggraph(
    tbl_graph_object,
    layout = "fr"
  ) +
    geom_edge_link(
      colour = "grey80",
      alpha = 0.6
    ) +
    geom_node_point(
      aes(
        color = circulation,
        shape = type
      ),
      size = 5
    ) +
    geom_node_text(
      aes(label = label),
      size = 3
    ) +
    scale_color_manual(
      values = c(
        "High circulation" = "#D55E00",
        "Single circulation" = "grey70",
        "Host" = "#0072B2"
      )
    ) +
    scale_shape_manual(
      values = c(
        Virus = 15,
        Host = 16
      )
    ) +
    theme_void()

  print(p)

  ggsave(
    file.path(
      output_dir,
      paste0(output_prefix, ".pdf")
    ),
    p,
    width = 20,
    height = 18
  )

  ggsave(
    file.path(
      output_dir,
      paste0(output_prefix, ".png")
    ),
    p,
    width = 20,
    height = 18,
    dpi = 600
  )

  write.table(
    edges,
    file.path(
      output_dir,
      paste0(output_prefix, "_edges.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    nodes,
    file.path(
      output_dir,
      paste0(output_prefix, "_nodes.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    virus_host_count,
    file.path(
      output_dir,
      paste0(output_prefix, "_virus_host_counts.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(
    list(
      plot = p,
      edges = edges,
      nodes = nodes,
      virus_host_count = virus_host_count
    )
  )
}

association_all <- prepare_association_data(df)


association_all_result <- plot_association_network(
  association_all,
  "virus_circulation_network_all"
)
