# ============================================================
# STAT3 Independent Dataset Validation
# One KG: STAT3
# Duplicate genes averaged + clear significance labels
# ============================================================

rm(list = ls())

# ============================================================
# 1. USER SETTINGS
# ============================================================

target_gene <- "STAT3"

output_folder_name <- paste0(target_gene, "_Independent_Dataset_Validation_Output")

plot_title <- paste0(target_gene, " expression validation across independent GEO datasets")
y_axis_title <- paste0(target_gene, " expression")

dpi_value <- 600

width_in  <- 13
height_in <- 5.5

# ============================================================
# 2. INSTALL AND LOAD PACKAGES
# ============================================================

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

packages <- c("ggplot2", "dplyr", "tidyr", "stringr", "readr", "tools", "grid")

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
# 3. SELECT CSV FILES
# ============================================================

cat("\nPlease select your CSV files.\n")
cat("Select all 4 GSE CSV files together if possible.\n\n")

select_csv_files <- function() {
  
  if (.Platform$OS.type == "windows") {
    files <- choose.files(
      caption = "Select CSV files",
      filters = matrix(
        c("CSV files", "*.csv",
          "All files", "*.*"),
        ncol = 2,
        byrow = TRUE
      ),
      multi = TRUE
    )
  } else {
    if (requireNamespace("tcltk", quietly = TRUE)) {
      files <- tcltk::tk_choose.files(
        caption = "Select CSV files",
        multi = TRUE,
        filters = matrix(
          c("CSV files", "*.csv",
            "All files", "*.*"),
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
# 5. PROCESS EACH DATASET
# ============================================================

process_dataset <- function(file_path, target_gene, output_dir) {
  
  dataset_id <- file_path_sans_ext(basename(file_path))
  
  cat("\n============================================================\n")
  cat("Processing dataset:", dataset_id, "\n")
  cat("============================================================\n")
  
  df <- read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = TRUE
  )
  
  # Detect gene column
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
  
  gene_col <- possible_gene_cols[possible_gene_cols %in% colnames(df)][1]
  
  if (is.na(gene_col)) {
    gene_col <- colnames(df)[1]
    message("No standard gene_symbol column found. First column used as gene column: ", gene_col)
  }
  
  df[[gene_col]] <- trimws(as.character(df[[gene_col]]))
  df[[gene_col]] <- toupper(df[[gene_col]])
  
  # Detect control and disease columns
  control_cols <- grep("^control(\\.|$)", colnames(df), value = TRUE, ignore.case = TRUE)
  disease_cols <- grep("^disease(\\.|$)", colnames(df), value = TRUE, ignore.case = TRUE)
  
  if (length(control_cols) == 0 || length(disease_cols) == 0) {
    stop(
      "Could not detect control and disease columns in ", dataset_id, ".\n",
      "Your sample columns should start with control and disease."
    )
  }
  
  sample_cols <- c(control_cols, disease_cols)
  
  # Convert sample columns to numeric
  df[sample_cols] <- lapply(df[sample_cols], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  
  cat("Rows before averaging:", nrow(df), "\n")
  cat("Duplicate gene rows:", sum(duplicated(df[[gene_col]])), "\n")
  cat("Control samples:", length(control_cols), "\n")
  cat("Disease samples:", length(disease_cols), "\n")
  
  # Average duplicate genes
  unique_df <- df %>%
    group_by(gene_symbol = .data[[gene_col]]) %>%
    summarise(
      across(
        all_of(sample_cols),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  # Replace NaN with NA
  unique_df <- unique_df %>%
    mutate(across(all_of(sample_cols), ~ ifelse(is.nan(.x), NA, .x)))
  
  cat("Unique genes after averaging:", nrow(unique_df), "\n")
  
  # Save unique averaged matrix
  unique_output_file <- file.path(
    output_dir,
    paste0(dataset_id, "_unique_gene_averaged_expression.csv")
  )
  
  write.csv(unique_df, unique_output_file, row.names = FALSE)
  
  # Extract target gene
  target_df <- unique_df %>%
    filter(gene_symbol == toupper(target_gene))
  
  if (nrow(target_df) == 0) {
    warning(target_gene, " not found in dataset: ", dataset_id)
    return(NULL)
  }
  
  # Convert to long format
  long_df <- target_df %>%
    select(gene_symbol, all_of(sample_cols)) %>%
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
    filter(!is.na(Group), !is.na(Expression))
  
  return(long_df)
}

# ============================================================
# 6. PROCESS ALL FILES
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

all_long_df$Group <- factor(all_long_df$Group, levels = c("Control", "Disease"))

dataset_order <- file_path_sans_ext(basename(input_files))
dataset_order <- dataset_order[dataset_order %in% unique(all_long_df$Dataset)]
all_long_df$Dataset <- factor(all_long_df$Dataset, levels = dataset_order)

combined_expression_file <- file.path(
  output_dir,
  paste0(target_gene, "_combined_expression_long_format.csv")
)

write.csv(all_long_df, combined_expression_file, row.names = FALSE)

# ============================================================
# 7. WILCOXON TEST
# ============================================================

stat_table <- all_long_df %>%
  group_by(Dataset) %>%
  summarise(
    n_control = sum(Group == "Control"),
    n_disease = sum(Group == "Disease"),
    median_control = median(Expression[Group == "Control"], na.rm = TRUE),
    median_disease = median(Expression[Group == "Disease"], na.rm = TRUE),
    mean_control = mean(Expression[Group == "Control"], na.rm = TRUE),
    mean_disease = mean(Expression[Group == "Disease"], na.rm = TRUE),
    direction = case_when(
      median_disease > median_control ~ "Upregulated in disease",
      median_disease < median_control ~ "Downregulated in disease",
      TRUE ~ "No median change"
    ),
    p_value = tryCatch(
      wilcox.test(Expression ~ Group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
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
      TRUE ~ paste0("FDR = ", signif(p_adj_BH, 3))
    ),
    final_label = paste0(significance, "\n", fdr_label)
  )

statistics_file <- file.path(
  output_dir,
  paste0(target_gene, "_Wilcoxon_BH_statistics.csv")
)

write.csv(stat_table, statistics_file, row.names = FALSE)

cat("\nStatistical summary:\n")
print(stat_table)

# ============================================================
# 8. ANNOTATION POSITION
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
    y_range = ifelse(y_range == 0, abs(y_max) + 1, y_range),
    
    # Label will stay clearly inside each panel
    y_pos = y_max + 0.18 * y_range,
    x_pos = 1.5
  ) %>%
  left_join(
    stat_table %>% select(Dataset, final_label),
    by = "Dataset"
  )

# ============================================================
# 9. PUBLICATION-QUALITY BOXPLOT
# ============================================================

p <- ggplot(all_long_df, aes(x = Group, y = Expression, fill = Group)) +
  
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
    nrow = 1,
    scales = "free_y"
  ) +
  
  # Clear significance label
  geom_label(
    data = anno_df,
    aes(x = x_pos, y = y_pos, label = final_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.7,
    lineheight = 0.88,
    label.size = 0,
    label.padding = unit(0.18, "lines"),
    fill = "white",
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c("Control" = "#4DBBD5", "Disease" = "#E64B35")
  ) +
  
  scale_color_manual(
    values = c("Control" = "#4DBBD5", "Disease" = "#E64B35")
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.35))
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
    
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(
      face = "bold",
      size = 13,
      color = "black"
    ),
    
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# ============================================================
# 10. SAVE FIGURES
# ============================================================

tiff_file <- file.path(
  output_dir,
  paste0(target_gene, "_independent_dataset_validation_boxplot_clear_label.tiff")
)

png_file <- file.path(
  output_dir,
  paste0(target_gene, "_independent_dataset_validation_boxplot_clear_label.png")
)

pdf_file <- file.path(
  output_dir,
  paste0(target_gene, "_independent_dataset_validation_boxplot_clear_label.pdf")
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
# 11. FINAL MESSAGE
# ============================================================

cat("\n============================================================\n")
cat("Analysis completed successfully!\n")
cat("============================================================\n")

cat("\nOutput folder:\n", output_dir, "\n")

cat("\nSaved files:\n")
cat("1. Combined STAT3 expression file:\n", combined_expression_file, "\n")
cat("2. Wilcoxon + BH statistics file:\n", statistics_file, "\n")
cat("3. TIFF figure:\n", tiff_file, "\n")
cat("4. PNG figure:\n", png_file, "\n")
cat("5. PDF figure:\n", pdf_file, "\n")