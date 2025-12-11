
```{r}
## 1. Load preprocessed dataframe
df_imputed <- read_csv("df_imputed.csv")

## 2. Make column names more uniform 
colnames(df_imputed)[colnames(df_imputed) == "label_clinical stage"] <- "label_clinical_stage"
colnames(df_imputed)[colnames(df_imputed) == "label_site label"] <- "label_site_label"

## 3. Remove non-informative rows
df_clean <- df_imputed[
  !(df_imputed$label_age == "unknown" |
    is.na(df_imputed$`label_site_label`)), 
]

## 4. Filter for colon samples
df_filtered_colon <- df_clean %>%
  filter(label_tissue == "Colon_and_rectum_cancer")

## 5. Convert age to numeric value
df_filtered_colon$label_age <- as.numeric(as.character(df_filtered_colon$label_age))

## 6. Select gene with maximum methylation variance across gene promotors
genes_only <- df_filtered_colon[, -(1:8)]
dim(genes_only)
variances <- apply(genes_only, 2, var)
top_gene <- names(which.max(variances))
top_gene

# --> HOXB6_1

## 7. Change site_label to binary value
df_finaal <- df_filtered_colon%>%
  mutate(label_site_label = recode(label_site_label, "Primary" = 0, "Metastatic" = 1, .default = NA_real_))

## 8. Linear Regression Models

## 8.1 Individual (No interactions)
model_1 <- glm(
  label_site_label ~ label_age + label_gender + HOXB6_1,
  data = df_finaal,
  family = binomial
)
summary(model_1)

## 8.2 Without_gender (Not significant)
model_2 <- glm(
  label_site_label ~ label_age  + HOXB6_1,
  data = df_finaal,
  family = binomial
)
summary(model_2)

## 8.3 All interactions
model_interactie <- glm(
  label_site_label ~ label_age * HOXB6_1 + label_age + HOXB6_1 + label_gender + label_gender*label_age + label_gender*HOXB6_1,
  data = df_finaal,
  family = binomial
)
```
summary(model_interactie)
