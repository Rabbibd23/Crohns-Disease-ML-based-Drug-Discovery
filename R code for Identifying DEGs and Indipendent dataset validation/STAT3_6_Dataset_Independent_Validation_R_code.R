# ============================================================
# STAT3 Independent Validation Across 6 GEO Datasets
# Final corrected version
#
# Expected panel order:
# Top row:    GSE66407 | GSE179285 | GSE95095
# Bottom row: GSE75214 | GSE102133 | GSE57945
#
# NOTE:
# The input file may be named GSE102134.csv, but the mRNA expression
# subseries with 12 controls + 65 CD samples is GSE102133. Therefore,
# the figure/manuscript label is shown as GSE102133.
# ============================================================

rm(list = ls())

# ============================================================
# 1. USER SETTINGS
# ============================================================

target_gene <- "STAT3"

output_folder_name <- paste0(
  target_gene,
  "_6_Dataset_Independent_Validation_Output"
)

plot_title <- paste0(
  target_gene,
  " expression validation across independent GEO datasets"
)

y_axis_title <- paste0(target_gene, " expression")

dpi_value <- 600

# 3 columns x 2 rows; large enough for manuscript use
width_in  <- 12.5
height_in <- 9.0

# ============================================================
# 2. INSTALL AND LOAD PACKAGES
# ============================================================

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "readr",
  "tools",
  "grid"
)

invisible(lapply(packages, install_if_missing))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tools)
  library(grid)
})

# ============================================================
# 3. SELECT THE 6 FINAL EXPRESSION CSV FILES
# ============================================================

cat("\nPlease select the SIX final expression CSV files.\n")
cat("Do not select metadata, sensitivity, summary, or STAT3-only files.\n\n")

select_csv_files <- function() {
  if (.Platform$OS.type == "windows") {
    files <- choose.files(
      caption = "Select the 6 final GEO expression CSV files",
      filters = matrix(
        c(
          "CSV files", "*.csv",
          "All files", "*.*"
        ),
        ncol = 2,
        byrow = TRUE
      ),
      multi = TRUE
    )
  } else {
    if (requireNamespace("tcltk", quietly = TRUE)) {
      files <- tcltk::tk_choose.files(
        caption = "Select the 6 final GEO expression CSV files",
        multi = TRUE,
        filters = matrix(
          c(
            "CSV files", "*.csv",
            "All files", "*.*"
          ),
          ncol = 2,
          byrow = TRUE
        )
      )
    } else {
      files <- file.choose()
    }
  }

  return(files)
}

input_files <- select_csv_files()

if (length(input_files) == 0) {
  stop("No CSV file selected.")
}

if (length(input_files) != 6) {
  stop(
    paste0(
      "Exactly 6 expression CSV files are required. You selected ",
      length(input_files),
      "."
    )
  )
}

cat("\nSelected files:\n")
print(input_files)

# ============================================================
# 4. CREATE OUTPUT FOLDER
# ============================================================

base_dir <- dirname(input_files[1])
output_dir <- file.path(base_dir, output_folder_name)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\nOutput folder:\n")
cat(output_dir, "\n")

# ============================================================
# 5. DATASET LABEL CLEANING
# ============================================================

get_dataset_label <- function(file_path) {
  file_name <- basename(file_path)

  dataset_id <- stringr::str_extract(
    file_name,
    "GSE[0-9]+"
  )

  if (is.na(dataset_id)) {
    dataset_id <- tools::file_path_sans_ext(file_name)
  }

  # GSE102134 is the parent SuperSeries.
  # The mRNA expression subseries used here is GSE102133.
  if (identical(dataset_id, "GSE102134")) {
    dataset_id <- "GSE102133"
  }

  return(dataset_id)
}

# Fixed figure order: old three on top, new three below
preferred_order <- c(
  "GSE66407",
  "GSE179285",
  "GSE95095",
  "GSE75214",
  "GSE102133",
  "GSE57945"
)

# ============================================================
# 6. PROCESS EACH DATASET
# ============================================================

process_dataset <- function(file_path, target_gene, output_dir) {

  dataset_id <- get_dataset_label(file_path)

  cat("\n============================================================\n")
  cat("Processing dataset:", dataset_id, "\n")
  cat("Input file:", basename(file_path), "\n")
  cat("============================================================\n")

  df <- read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = TRUE
  )

  # ----------------------------------------------------------
  # Detect gene-symbol column
  # ----------------------------------------------------------

  possible_gene_cols <- c(
    "gene_symbol",
    "Gene_symbol",
    "GENE_SYMBOL",
    "official_gene_symbol",
    "Gene",
    "gene",
    "Symbol",
    "symbol"
  )

  gene_col <- possible_gene_cols[
    possible_gene_cols %in% colnames(df)
  ][1]

  if (length(gene_col) == 0 || is.na(gene_col)) {
    gene_col <- colnames(df)[1]
    message(
      "No standard gene_symbol column found. First column used: ",
      gene_col
    )
  }

  df[[gene_col]] <- toupper(
    trimws(
      as.character(df[[gene_col]])
    )
  )

  # ----------------------------------------------------------
  # Detect Control and Disease columns
  # ----------------------------------------------------------

  control_cols <- grep(
    "^control(\\.|$)",
    colnames(df),
    value = TRUE,
    ignore.case = TRUE
  )

  disease_cols <- grep(
    "^disease(\\.|$)",
    colnames(df),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(control_cols) == 0 || length(disease_cols) == 0) {
    stop(
      paste0(
        "Could not detect control/disease columns in ",
        dataset_id,
        ". Sample columns must start with control or disease."
      )
    )
  }

  sample_cols <- c(control_cols, disease_cols)

  # ----------------------------------------------------------
  # Convert expression columns to numeric
  # ----------------------------------------------------------

  df[sample_cols] <- lapply(
    df[sample_cols],
    function(x) {
      suppressWarnings(
        as.numeric(as.character(x))
      )
    }
  )

  cat("Rows before duplicate-gene averaging:", nrow(df), "\n")
  cat("Duplicate gene rows:", sum(duplicated(df[[gene_col]])), "\n")
  cat("Control samples:", length(control_cols), "\n")
  cat("Disease samples:", length(disease_cols), "\n")

  # ----------------------------------------------------------
  # Remove empty gene symbols
  # ----------------------------------------------------------

  df <- df %>%
    filter(
      !is.na(.data[[gene_col]]),
      .data[[gene_col]] != ""
    )

  # ----------------------------------------------------------
  # Average duplicated genes per sample
  # ----------------------------------------------------------

  unique_df <- df %>%
    group_by(gene_symbol = .data[[gene_col]]) %>%
    summarise(
      across(
        all_of(sample_cols),
        ~ {
          if (all(is.na(.x))) {
            NA_real_
          } else {
            mean(.x, na.rm = TRUE)
          }
        }
      ),
      .groups = "drop"
    ) %>%
    mutate(
      across(
        all_of(sample_cols),
        ~ ifelse(is.nan(.x), NA_real_, .x)
      )
    )

  cat("Unique genes after averaging:", nrow(unique_df), "\n")

  # Save duplicate-averaged matrix for reproducibility
  unique_output_file <- file.path(
    output_dir,
    paste0(
      dataset_id,
      "_unique_gene_averaged_expression.csv"
    )
  )

  write.csv(
    unique_df,
    unique_output_file,
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # Extract target gene
  # ----------------------------------------------------------

  target_df <- unique_df %>%
    filter(
      gene_symbol == toupper(target_gene)
    )

  if (nrow(target_df) == 0) {
    stop(
      paste0(
        target_gene,
        " not found in dataset ",
        dataset_id,
        "."
      )
    )
  }

  if (nrow(target_df) != 1) {
    stop(
      paste0(
        "Expected one ",
        target_gene,
        " row after duplicate averaging in ",
        dataset_id,
        ", but found ",
        nrow(target_df),
        "."
      )
    )
  }

  # ----------------------------------------------------------
  # Convert target-gene row to long format
  # ----------------------------------------------------------

  long_df <- target_df %>%
    select(
      gene_symbol,
      all_of(sample_cols)
    ) %>%
    pivot_longer(
      cols = -gene_symbol,
      names_to = "Sample",
      values_to = "Expression"
    ) %>%
    mutate(
      Dataset = dataset_id,
      Group = case_when(
        grepl("^control", Sample, ignore.case = TRUE) ~ "Control",
        grepl("^disease", Sample, ignore.case = TRUE) ~ "Disease",
        TRUE ~ NA_character_
      ),
      Gene = gene_symbol,
      Expression = as.numeric(Expression)
    ) %>%
    filter(
      !is.na(Group),
      !is.na(Expression),
      is.finite(Expression)
    )

  if (sum(long_df$Group == "Control") == 0 ||
      sum(long_df$Group == "Disease") == 0) {
    stop(
      paste0(
        "No usable Control/Disease expression values remained for ",
        dataset_id,
        "."
      )
    )
  }

  return(long_df)
}

# ============================================================
# 7. PROCESS ALL SIX DATASETS
# ============================================================

all_long_list <- lapply(
  input_files,
  process_dataset,
  target_gene = target_gene,
  output_dir = output_dir
)

all_long_df <- bind_rows(all_long_list)

if (nrow(all_long_df) == 0) {
  stop("No STAT3 expression data found in selected files.")
}

# Group order
all_long_df$Group <- factor(
  all_long_df$Group,
  levels = c("Control", "Disease")
)

# Dataset order for 3 x 2 figure
present_datasets <- unique(as.character(all_long_df$Dataset))

dataset_order <- c(
  preferred_order[preferred_order %in% present_datasets],
  setdiff(present_datasets, preferred_order)
)

all_long_df$Dataset <- factor(
  all_long_df$Dataset,
  levels = dataset_order
)

# Save combined target-gene long-format data
combined_expression_file <- file.path(
  output_dir,
  paste0(
    target_gene,
    "_6_dataset_combined_expression_long_format.csv"
  )
)

write.csv(
  all_long_df,
  combined_expression_file,
  row.names = FALSE
)

# ============================================================
# 8. WILCOXON RANK-SUM TEST + BH CORRECTION ACROSS 6 DATASETS
# ============================================================

stat_table <- all_long_df %>%
  group_by(Dataset) %>%
  summarise(
    n_control = sum(Group == "Control"),
    n_disease = sum(Group == "Disease"),

    median_control = median(
      Expression[Group == "Control"],
      na.rm = TRUE
    ),

    median_disease = median(
      Expression[Group == "Disease"],
      na.rm = TRUE
    ),

    mean_control = mean(
      Expression[Group == "Control"],
      na.rm = TRUE
    ),

    mean_disease = mean(
      Expression[Group == "Disease"],
      na.rm = TRUE
    ),

    median_difference = median_disease - median_control,

    direction = case_when(
      median_disease > median_control ~ "Upregulated in disease",
      median_disease < median_control ~ "Downregulated in disease",
      TRUE ~ "No median change"
    ),

    p_value = tryCatch(
      wilcox.test(
        Expression ~ Group,
        alternative = "two.sided",
        exact = FALSE,
        correct = TRUE
      )$p.value,
      error = function(e) NA_real_
    ),

    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),

    significance = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    ),

    fdr_label = case_when(
      is.na(p_adj_BH) ~ "FDR = NA",
      p_adj_BH < 0.001 ~ "FDR < 0.001",
      TRUE ~ paste0(
        "FDR = ",
        signif(p_adj_BH, 3)
      )
    ),

    final_label = paste0(
      significance,
      "\n",
      fdr_label
    )
  ) %>%
  arrange(
    match(
      as.character(Dataset),
      dataset_order
    )
  )

statistics_file <- file.path(
  output_dir,
  paste0(
    target_gene,
    "_6_dataset_Wilcoxon_BH_statistics.csv"
  )
)

write.csv(
  stat_table,
  statistics_file,
  row.names = FALSE
)

cat("\nFinal statistical summary:\n")
print(stat_table)

# ============================================================
# 9. ANNOTATION POSITION FOR EACH FACET
# ============================================================

anno_df <- all_long_df %>%
  group_by(Dataset) %>%
  summarise(
    y_max = max(Expression, na.rm = TRUE),
    y_min = min(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_range = ifelse(
      y_range == 0,
      abs(y_max) + 1,
      y_range
    ),
    y_pos = y_max + 0.18 * y_range,
    x_pos = 1.5
  ) %>%
  left_join(
    stat_table %>%
      select(Dataset, final_label),
    by = "Dataset"
  )

# ============================================================
# 10. PUBLICATION-QUALITY 3 x 2 BOXPLOT
#     Original colors/style retained
# ============================================================

p <- ggplot(
  all_long_df,
  aes(
    x = Group,
    y = Expression,
    fill = Group
  )
) +

  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.92,
    linewidth = 0.65,
    color = "black"
  ) +

  geom_jitter(
    aes(color = Group),
    width = 0.12,
    size = 1.65,
    alpha = 0.65,
    show.legend = FALSE
  ) +

  facet_wrap(
    ~ Dataset,
    ncol = 3,
    nrow = 2,
    scales = "free_y"
  ) +

  geom_label(
    data = anno_df,
    aes(
      x = x_pos,
      y = y_pos,
      label = final_label
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.7,
    lineheight = 0.88,
    label.size = 0,
    label.padding = unit(0.18, "lines"),
    fill = "white",
    color = "black"
  ) +

  # ORIGINAL COLORS RETAINED
  scale_fill_manual(
    values = c(
      "Control" = "#4DBBD5",
      "Disease" = "#E64B35"
    )
  ) +

  scale_color_manual(
    values = c(
      "Control" = "#4DBBD5",
      "Disease" = "#E64B35"
    )
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(0.06, 0.35)
    )
  ) +

  labs(
    title = plot_title,
    x = "",
    y = y_axis_title
  ) +

  theme_classic(base_size = 14) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5,
      color = "black"
    ),

    axis.title.y = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),

    axis.text.x = element_text(
      face = "bold",
      size = 13,
      color = "black"
    ),

    axis.text.y = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),

    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.8
    ),

    strip.text = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),

    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),

    # Legend stays at the top, but the color boxes are larger
    # and the text is larger/bold for manuscript readability.
    legend.position = "top",
    legend.title = element_blank(),

    legend.text = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),

    legend.key.width = unit(1.05, "cm"),
    legend.key.height = unit(0.70, "cm"),
    legend.spacing.x = unit(0.20, "cm"),

    plot.margin = margin(
      10,
      10,
      10,
      10
    )
  ) +

  guides(
    fill = guide_legend(
      override.aes = list(
        alpha = 1
      )
    )
  )

print(p)

# ============================================================
# 11. SAVE FIGURES
# ============================================================

tiff_file <- file.path(
  output_dir,
  paste0(
    target_gene,
    "_6_dataset_independent_validation_3x2.tiff"
  )
)

png_file <- file.path(
  output_dir,
  paste0(
    target_gene,
    "_6_dataset_independent_validation_3x2.png"
  )
)

pdf_file <- file.path(
  output_dir,
  paste0(
    target_gene,
    "_6_dataset_independent_validation_3x2.pdf"
  )
)

ggsave(
  filename = tiff_file,
  plot = p,
  device = "tiff",
  width = width_in,
  height = height_in,
  units = "in",
  dpi = dpi_value,
  compression = "lzw"
)

ggsave(
  filename = png_file,
  plot = p,
  device = "png",
  width = width_in,
  height = height_in,
  units = "in",
  dpi = dpi_value
)

ggsave(
  filename = pdf_file,
  plot = p,
  device = "pdf",
  width = width_in,
  height = height_in,
  units = "in"
)

# ============================================================
# 12. FINAL MESSAGE
# ============================================================

cat("\n============================================================\n")
cat("STAT3 SIX-DATASET VALIDATION COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nOutput folder:\n", output_dir, "\n")

cat("\nSaved files:\n")
cat("1. Combined STAT3 long-format data:\n", combined_expression_file, "\n")
cat("2. Wilcoxon + BH statistics:\n", statistics_file, "\n")
cat("3. TIFF figure:\n", tiff_file, "\n")
cat("4. PNG figure:\n", png_file, "\n")
cat("5. PDF figure:\n", pdf_file, "\n")

cat("\nFigure layout:\n")
cat("Top row: GSE66407 | GSE179285 | GSE95095\n")
cat("Bottom row: GSE75214 | GSE102133 | GSE57945\n")
