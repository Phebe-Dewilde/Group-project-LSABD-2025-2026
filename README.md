# Group-project-LSABD-2025-2026
## Data load
```{r}
library(tidyverse)
methylomics_data <- read_csv("colon_stomach_lung_cancer_loyfer_zhang_weightedmean_and_metadata.csv")
```
## Data pre-processing
### General data cleaning
1. We first deleted the column 'label_organism', because it was all homo sapiens and that did not matter. Then we've made from column 'label_gender' a binary column:
```{r}
datazonderleeftijd <- subset(methylomics_data, select = -label_organism)
binairdata <- datazonderleeftijd %>%
  mutate(label_gender = recode(label_gender, "M" = 0, "F" = 1, .default = NA_real_))
```
2. We added a column 'label_disease' so that the first 38 samples are considered healthy and the rest disease. We've put this column at the 6th place.
```{r}
binairdata$label_disease <- "Healthy"
binairdata$label_disease[39:nrow(binairdata)] <- "Disease"
binairdata<- binairdata %>% 
  relocate(label_disease, .before = 6)
```
3. We gave every value in column 'label_source_of_sample' a capital letter. In retrospect this was not necessary as we did not use it for any of our hypotheses.
```{r}
binairdata$label_source_of_sample <- as.character(binairdata$label_source_of_sample)
library(tools)
binairdata$label_source_of_sample <- toTitleCase(binairdata$label_source_of_sample)
```
4. To link gene promoters to their genomic positions, chromosomal location identifiers were replaced with the corresponding gene names as column headers. We used the bed file.
```{r}
enhancers <-read.table("C:/Users/phebe/OneDrive - UGent/Large Scale Data analysis/files/promoter_hg38_extended.bed", 
                       header = FALSE,    # BED-bestanden hebben meestal geen header
                       sep = "\t",        # BED is tab-gescheiden
                       stringsAsFactors = FALSE)

library(readr)
#tidyverse
write_csv(enhancers,  "promoter_hg38_extended.csv")

colnames(enhancers) <- enhancers[1, ]
enhancers <- enhancers[-1, ]

# Maak een unieke sleutel in de enhancers-tabel
enhancers$key <- paste(enhancers$chrom, enhancers$chromStart, enhancers$chromEnd, sep = "_")

# Maak een naamvertaalschema: sleutel -> enhancer naam
name_map <- setNames(enhancers$name, enhancers$key)

# Huidige kolomnamen van binairdata
old_names <- colnames(binairdata)

# Nieuwe kolomnamen: vervang als ze in het name_map voorkomen
new_names <- ifelse(old_names %in% names(name_map),
                    name_map[old_names],
                    old_names)

# Pas de kolomnamen aan
colnames(binairdata) <- new_names

# Controleer het resultaat
head(colnames(binairdata))
```
5. In columns 'label_tissue' and 'label_source_of_sample' we replaced the spaces with underscores, and also standardised the writing way of Lymphnode.
```{r}
library(dplyr)
library(stringr)

datazonderspaties <- datazonderspaties %>% 
  mutate(
    label_source_of_sample = label_source_of_sample |> 
      str_squish() |>                      # dubbele spaties weg
      str_to_lower() |>                    # alles lowercase
      str_replace_all(" ", "_")            # spaties naar underscores
  ) %>%
  mutate(
    label_source_of_sample = case_when(
      label_source_of_sample %in% c("lymph_node", "lymphnode", "lymph_nodes") ~ "Lymph_node",
      TRUE ~ str_to_title(label_source_of_sample)  # maakt vb lung_epithelium → Lung_epithelium
    )
  )
```
6. We gave the healthy samples in the column 'label_clinical stage', so when 'label_disease == Healthy', gave it stage 0.
```{r}
datazonderspaties$`label_clinical stage`[datazonderspaties$label_disease == "Healthy"] <- "0"

datazonderspaties$`label_site label`[datazonderspaties$label_disease == "Healthy"] <- "Healthy"
```
7. We addressed the missingness.
7.1. We've first made a visualisation of the percentages of missingness for ever CpG loci in our dataset. We've made 2 different visualisations:
```{r}
data_cpg <- datazonderspaties[, -(1:8)] #first 8 columns out because they are not CpG loci

missing_summary <- data_cpg %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = 'CpG_loci', values_to = 'missing_values') %>%
  mutate(missing_percentage = (missing_values / nrow(datazonderspaties)) *100) %>%
  arrange(desc(missing_percentage))
  
missing_summary %>%
  ggplot(aes(x = reorder(CpG_loci, missing_percentage), y = missing_percentage)) +
  geom_col(fill = 'steelblue') +
  coord_flip() +
  labs( title = 'Percentage of missing values by CpG loci', x = 'CpG loci', y = 'Percentage of missing values') +
  theme_minimal() +
    theme(axis.text.y = element_blank())

ggplot(missing_summary, aes(x = missing_percentage)) +
  geom_histogram(bins = 40) +
  theme_minimal() +
  labs(title = "Distribution of missingness",
       x = "% missing", y = "Number of CpG loci")
```
7.2. After that we've decided that the cut off should be at 25% missingness. Hereby were 7800 CpG loci deleted from the database 'datazonderspaties', and we've made a subset of that 'df_filtered'.
```{r}
cutoff <- 0.25  # 25%

na_frac <- colMeans(is.na(datazonderspaties))        # fractie NA per kolom
cols_keep <- names(na_frac)[na_frac <= cutoff]

df_filtered <- datazonderspaties[, cols_keep]
```
7.3. We wanted to an KNN imputation on the remaining NA's for some CpG loci, but the software of R was maxed out. Therefore we filtered the df_filtered dataset on standard variance. The onces with the lowest were filtered out because they could be householde genes and they don't change under certain conditions.
We came to a conclusion that the NAs in the CpG loci were MAR. The missingness is related to another variable. Whereas we stipulated that the NAs and 'unknown' in the metadata were MCAR, therefore we deleted the samples, but only if they would interfere with our hypotheses (see later).
```{r}
#divide the metadata from the CpG loci
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

meta_cols <- intersect(meta_cols, names(df_filtered))
gene_cols <- setdiff(names(df_filtered), meta_cols)

#calculate the stdev for all CpG loci
gene_sd <- sapply(df_filtered[, gene_cols], sd, na.rm = TRUE)

#make a histogram of the stdev, and based on that we can decide a cut off value (<0.05)

```{r}
library(ggplot2)

df_sd <- data.frame(sd = gene_sd)

ggplot(df_sd, aes(x = sd)) +
  geom_histogram(
    bins = 100,
    fill = "skyblue",
    color = "white"
  ) +
  labs(
    title = "Histogram van standaarddeviaties van genen",
    x = "Standaarddeviatie",
    y = "Aantal genen"
  ) +
  theme_minimal()

sd_cutoff <- 0.05
genes_to_keep   <- names(gene_sd)[gene_sd > sd_cutoff]
genes_to_remove <- names(gene_sd)[gene_sd <= sd_cutoff]

length(genes_to_keep)    # ~14.200 genen over
length(genes_to_remove)  # hoeveel eruit gaan (6000-tal)

df_sd_filtered <- df_filtered[, c(meta_cols, genes_to_keep)]
```
7.4. KNN-imputation in Python, because using Google Colab won't max out the services because we use their server.
``` Python
#importing the df_sd_filtered dataset

from google.colab import drive
drive.mount('/content/drive')

import pandas as pd

file_path = '/content/drive/MyDrive/Colab_Notebooks/df_sd_filtered.csv'  # pas mapnaam aan
df = pd.read_csv(file_path, sep=',')  # gebruik sep=',' als komma's gebruikt zijn
df.head()

#dividing the dataset in metadata and CpG loci, so only KNN-imputation occurs on the CpG loci

metadata_cols = ['label_title', 'label_tissue', 'label_age', 'label_gender', 'label_clinical stage', 'label_site label', 'label_source_of_sample', 'label_disease']  # pas aan
chrom_cols = [col for col in df_clean.columns if col not in metadata_cols]

metadata_df = df[metadata_cols]
chrom_df = df[chrom_cols]

#KNN imputation + download the new dataset so we can use it in R
from sklearn.impute import KNNImputer

imputer = KNNImputer(n_neighbors=5)
chrom_imputed = imputer.fit_transform(chrom_df)
chrom_imputed_df = pd.DataFrame(chrom_imputed, columns=chrom_cols)
df_imputed = pd.concat([metadata_df.reset_index(drop=True), chrom_imputed_df], axis=1)
df_imputed.to_csv('df_imputed.csv', index=False)
from google.colab import files
files.download('df_imputed.csv')
```
7.5. We controlled it in R if there were any NAs left and if the size of the dataset remained equal and if all columns are still there:
```{r}
genes_only <- df_imputed[, -(1:8)]   # verwijder de eerste 8 metadata kolommen
sum(is.na(genes_only))

dim(df_imputed)

setdiff(colnames(df_sd_filtered), colnames(df_imputed))
#output lege vector= GOED
```
### Data cleaning for the univariate analysis
### Data cleaning for the multivariate analysis
### Data cleaning for the machine learning algorithm
Non

## Univariate Analysis
## Multivariate Analysis
## Machine Learning Algorithm
