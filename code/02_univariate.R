## 02_univariate.R
## Univariate analysis of methylation in colon samples
## - load imputed dataset
## - subset colon samples (healthy vs cancer)
## - choose an APC2 CpG with highest variance
## - compare methylation between stage 0 and stage IV
## - generate summary statistics and plots

## 0. Load packages ----
library(tidyverse)
library(stringr)

##  1. Load imputed data ----
read.csv("df_imputed.csv", check.names = FALSE)

## 2. Univariate specific metadata cleaning ----

# 2.1 Cleaning up column names by adding underscores 
#-> was overlooked during the general preprocessing

colnames(df_imputed)[colnames(df_imputed) == "label_clinical stage"] <- "label_clinical_stage"
colnames(df_imputed)[colnames(df_imputed) == "label_site label"] <- "label_site_label"

# 2.2 Removal of NAs in label_clinical_stage column - Keep samples with a valid clinical stage
#Since the feature label_clinical_stage is considered for the univariate analysis,
#rows that have unknown values or NAs for label_clinical_stage are removed.

uni_df <- df_imputed[
  !is.na(df_imputed$label_clinical_stage) &                # no NA
  trimws(df_imputed$label_clinical_stage) != "" &          # no empty cells
  tolower(df_imputed$label_clinical_stage) != "unknown",   # no unknown
]
#check distribution of different stages 
table(uni_df$label_clinical_stage)

## 3. Restrict to colon samples (healthy colon + colon cancer) ----
df_filtered_colon <- uni_df %>%
  filter(label_tissue %in% c("Colon_and_rectum_cancer",
                             "Colon"))

## 4. Choose a gene for univariate analysis ----
# Option 1: top variable gene overall (ultimately not chosen to be used in report, but kept for reference)
variances <- apply(genes_only, 2, var)
top_gene_overall <- names(which.max(variances))

# Option 2: restrict to APC-related CpGs= biologial relevance (as described in the report)

apc_cols <- grep("^APC_", names(df_filtered_colon), value = TRUE)
grep("APC", names(df_filtered_colon), value = TRUE)
#multiple CpG loci with APC in name are present --> filter on real APC genes only

apc2_cols <- grep("^APC2_", names(df_filtered_colon), value = TRUE)
apc2_cols
apc2_data <- df_filtered_colon[, apc2_cols]
apc2_variances <- apply(apc2_data, 2, var)
APC2_top_gene <- names(which.max(apc2_variances))
APC2_top_gene

## 5. Build univariate analysis dataset: stage 0 vs stage IV ----
uni_df_stat <- df_filtered_colon %>%
  filter(label_clinical_stage %in% c("0", "Ⅳ")) %>%    # stage 0 and stage IV
  mutate(
    group = factor(
      label_clinical_stage,
      levels = c("0", "Ⅳ"),
      labels = c("Stage_0", "Stage_IV")
    )
  )
## 6. QQ-plot diagnostics ---- 
# QQ-plot included for diagnostics only, not to justify normality.
# Beta-values are bounded (0–1) and cannot be normally distributed.
# Deviations in the tails confirm non-normality → Wilcoxon is appropriate.
qqnorm(uni_df_stat[[APC2_top_gene]])
qqline(uni_df_stat[[APC2_top_gene]], col = "red", lwd = 2)

## 7. Descriptive statistics by group ----
uni_summary <- uni_df_stat %>%
  group_by(group) %>%
  summarise(
    n      = n(),
    mean   = mean(.data[[APC2_top_gene]], na.rm = TRUE),
    sd     = sd(.data[[APC2_top_gene]], na.rm = TRUE),
    median = median(.data[[APC2_top_gene]], na.rm = TRUE),
    IQR    = IQR(.data[[APC2_top_gene]], na.rm = TRUE)
  )

uni_summary

## 8. Outlier inspection (visual only)

# Outliers were inspected visually using boxplots.
# Because methylation beta-values represent biological variability and
# the Wilcoxon/Mann–Whitney test is robust to outliers, no samples were removed.
ggplot(uni_df_stat, aes(y = .data[[APC2_top_gene]], x = "")) +
  geom_boxplot(outlier.colour = "red", alpha = 0.6) +
  labs(title = paste("Outlier visualization for", APC2_top_gene),
       y = "Methylation (beta-value)",
       x = "")

# Outlier check: some outliers are biologically normal.
# Mann–Whitney U (=Wilcoxon rank-sum) is robust → no removal necessary.


## 9.Wilcoxon rank-sum test (robust to non-normality) ----
wilcox_res <- wilcox.test(
  uni_df_stat[[APC2_top_gene]] ~ uni_df_stat$group,
  exact = FALSE
)

wilcox_res

## 10. Visualisations ----
# 10.1 Boxplot with jittered points
ggplot(uni_df_stat, aes(
  x = group,
  y = .data[[APC2_top_gene]],
  fill = group
)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.08, alpha = 0.6, size = 2, color = "black") +
  scale_fill_manual(values = c("#6AAED6", "#CC4C4C")) +
  labs(
    title = paste("Methylation of", APC2_top_gene, "in stage 0 VS stage IV"),
    x = "Clinical stage",
    y = "Methylation (beta-value)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# 10.2 Violin + boxplot + jitter 
ggplot(uni_df_stat, aes(
  x = group,
  y = .data[[APC2_top_gene]],
  fill = group
)) +
  geom_violin(trim = TRUE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.08, alpha = 0.6, size = 2, color = "black") +
  scale_fill_manual(values = c("#6AAED6", "#CC4C4C")) +
  labs(
    title = paste("Methylation of", APC2_top_gene, "in Stage 0 vs Stage IV"),
    x = "Clinical stage",
    y = "Methylation (beta-value)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


