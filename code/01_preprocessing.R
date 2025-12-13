## 01_preprocessing.R
## Preprocessing of methylomics data for LSABD group project
## This script:
## - loads the raw data
## - cleans and recodes metadata
## - links CpG loci to promoter annotation 
## - explores and filters missingness
## - filters low-variance CpG loci
## - exports a filtered dataset for KNN imputation in Python

## 0. Load packages ----
library(tidyverse)
library(stringr)

## 1. Load raw data ----
raw_data <- read_csv(
  "colon_stomach_lung_cancer_loyfer_zhang_weightedmean_and_metadata.csv")

## 2. Basic metadata cleaning ----

# 2.1 Remove non-informative column
# 'label_organism' is always 'homo sapiens' → not informative
data_clean <- raw_data %>%
  select(-label_organism)

# 2.2 Recode gender to binary: M = 0, F = 1
data_clean <- data_clean %>%
  mutate(
    label_gender = recode(
      label_gender,
      "M" = 0,
      "F" = 1,
      .default = NA_real_
    )
  )

# 2.3 Create new column: label_disease  (Healthy vs Disease)
# First 38 samples are considered Healthy, the rest are Disease
data_clean <- data_clean %>%
  mutate(
    label_disease = if_else(row_number() <= 38, "Healthy", "Disease")
  ) %>%
  # Place label_disease as 6th column for better visualisation of metadata
  relocate(label_disease, .before = 6)

# 2.4 Standardise label_tissue and label_source_of_sample columns
# - trim extra spaces
# - convert to lowercase
# - replace spaces with underscores
# - harmonise different spellings of "lymph node"

data_clean <- data_clean %>%
  mutate(
    label_tissue = label_tissue |>
      str_squish() |>
      str_to_lower() |>
      str_replace_all(" ", "_"),
    label_source_of_sample = label_source_of_sample |>
      str_squish() |>
      str_to_lower() |>
      str_replace_all(" ", "_")
  ) %>%
  mutate(
    label_source_of_sample = case_when(
      label_source_of_sample %in% c("lymph_node", "lymphnode", "lymph_nodes") ~ "lymph_node",
      TRUE ~ label_source_of_sample
    )
  )

# 2.5 Assign clinical stage and site label for Healthy samples
# Healthy samples get clinical stage 0 and site label "Healthy"
data_clean <- data_clean %>%
  mutate(
    `label_clinical stage` = if_else(
      label_disease == "Healthy",
      "0",
      `label_clinical stage`
    ),
    `label_site label` = if_else(
      label_disease == "Healthy",
      "Healthy",
      `label_site label`
    )
  )

## 3. Map CpG loci to promoter names using BED annotation ----
## This step links genomic coordinates in the dataset to gene/promoter names.
## It uses the provided BED file with promoter annotations.

enhancers <- read.table(
  "promoter_hg38_extended.bed",
  header = FALSE,        # BED files typically have no header row
  sep = "\t", #BED is tab-separated
  stringsAsFactors = FALSE
)

# Convert first BED row into column names 
colnames(enhancers) <- enhancers[1, ]
enhancers <- enhancers[-1, ]

# Create key "chrom_chromStart_chromEnd" to match our CpG column names
enhancers$key <- paste(
  enhancers$chrom,
  enhancers$chromStart,
  enhancers$chromEnd,
  sep = "_"
)

# Create mapping: key → promoter/gene name
name_map <- setNames(enhancers$name, enhancers$key)

# Replace CpG locus column names in data_clean using our mapping
old_names <- colnames(data_clean)

new_names <- ifelse(
  old_names %in% names(name_map),
  name_map[old_names],
  old_names
)

colnames(data_clean) <- new_names

# Check that no genomic coordinates remain -> must be zero= mapping worked:
sum(str_detect(colnames(data_clean), "^chr"))

## 4. Define metadata vs CpG columns ----

meta_cols <- c(
  "label_title",
  "label_gender",
  "label_age",
  "label_tissue",
  "label_source_of_sample",
  "label_disease",
  "label_clinical stage",
  "label_site label"
)

# Keep only metadata columns that actually exist
meta_cols <- intersect(meta_cols, names(data_clean))

# CpG / gene columns = all other columns
gene_cols <- setdiff(names(data_clean), meta_cols)

## 5. Explore missingness per CpG locus ----

# Subset to CpG/gene columns only
data_genes <- data_clean[, gene_cols]

# Summarise missing values per CpG locus
missing_summary <- data_genes %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(
    everything(),
    names_to = "CpG_locus",
    values_to = "missing_values"
  ) %>%
  mutate(
    missing_percentage = (missing_values / nrow(data_clean)) * 100
  ) %>%
  arrange(desc(missing_percentage))

# Example plots (for the notebook; can be run here too)
# library(ggplot2)
#
# missing_summary %>%
#   ggplot(aes(x = reorder(CpG_locus, missing_percentage),
#              y = missing_percentage)) +
#   geom_col(fill = "steelblue") +
#   coord_flip() +
#   labs(
#     title = "Percentage of missing values by CpG locus",
#     x = "CpG locus",
#     y = "Percentage of missing values"
#   ) +
#   theme_minimal() +
#   theme(axis.text.y = element_blank())
#
# ggplot(missing_summary, aes(x = missing_percentage)) +
#   geom_histogram(bins = 40, fill = "steelblue", color = "white") +
#   theme_minimal() +
#   labs(
#     title = "Distribution of missingness across CpG loci",
#     x = "% missing",
#     y = "Number of CpG loci"
#   )

## 6. Filter CpG loci with > 25% missingness ----

cutoff <- 0.25  # 25%

na_fraction <- colMeans(is.na(data_genes))  # fraction of NA per CpG locus
genes_keep_missing <- names(na_fraction)[na_fraction <= cutoff]

# Keep all metadata columns + CpG loci that pass the missingness filter
df_filtered <- data_clean[, c(meta_cols, genes_keep_missing)]

## 7. Filter low-variance CpG loci ----

# Recalculate gene columns for the filtered dataset
gene_cols_filtered <- setdiff(names(df_filtered), meta_cols)

# Standard deviation per CpG locus
gene_sd <- sapply(df_filtered[, gene_cols_filtered], sd, na.rm = TRUE)

# Optional: inspect SD distribution (for notebook)
# df_sd <- data.frame(sd = gene_sd)
# ggplot(df_sd, aes(x = sd)) +
#   geom_histogram(bins = 100, fill = "skyblue", color = "white") +
#   labs(
#     title = "Histogram of standard deviations of CpG loci",
#     x = "Standard deviation",
#     y = "Number of loci"
#   ) +
#   theme_minimal()

sd_cutoff <- 0.05
genes_keep_sd   <- names(gene_sd)[gene_sd > sd_cutoff]
genes_remove_sd <- names(gene_sd)[gene_sd <= sd_cutoff]

length(genes_keep_sd)    # how many CpG loci remain
length(genes_remove_sd)  # how many are removed

df_sd_filtered <- df_filtered[, c(meta_cols, genes_keep_sd)]

## 8. Export dataset for KNN imputation in Python ----

# This file will be used in a Python/Colab notebook for KNN imputation.

write_csv(df_sd_filtered, "df_sd_filtered.csv")

```{python}
from google.colab import drive
drive.mount('/content/drive')

import pandas as pd

file_path ='/content/drive/MyDrive/Colab_Notebooks/df_sd_filtered.csv'
df = pd.read_csv(file_path)

metadata_cols = ['label_title', 'label_tissue', 'label_age', 'label_gender',
                 'label_clinical stage', 'label_site label',
                 'label_source_of_sample', 'label_disease']

chrom_cols = [col for col in df.columns if col not in metadata_cols]

metadata_df = df[metadata_cols]
chrom_df = df[chrom_cols]

from sklearn.impute import KNNImputer
imputer = KNNImputer(n_neighbors=5)
chrom_imputed = imputer.fit_transform(chrom_df)

chrom_imputed_df = pd.DataFrame(chrom_imputed, columns=chrom_cols)
df_imputed = pd.concat([metadata_df.reset_index(drop=True), chrom_imputed_df], axis=1)

df_imputed.to_csv('df_imputed.csv', index=False)
```

## 9. (Optional) Post-imputation checks in R ----
## After running KNN imputation in Python, you will get 'df_imputed.csv'.
## The following code can be used to verify that:
## - no NA values remain in the CpG loci
## - the dimensions are unchanged
## - all columns are still present.

df_imputed <- read_csv("df_imputed.csv")

# # Separate gene columns
gene_cols_imputed <- setdiff(names(df_imputed), meta_cols)
genes_only <- df_imputed[, gene_cols_imputed]
#
# 1) Check if any NA remain
sum(is.na(genes_only))

# 2) Check dimensions
dim(df_imputed)

# 3) Check if all original columns are present
setdiff(colnames(df_sd_filtered), colnames(df_imputed))
# # An empty character vector means: all columns are preserved (GOOD).



