

## 1. Load preprocessed dataframe
df_imputed <- read_csv("df_imputed.csv")

## 2. Make column names more uniform 
colnames(df_imputed)[colnames(df_imputed) == "label_clinical stage"] <- "label_clinical_stage"
colnames(df_imputed)[colnames(df_imputed) == "label_site label"] <- "label_site_label"

## 3. Remove non-informative rows
df_clean <- df_imputed[
  !(df_imputed$label_age == "unknown" |
    is.na(df_imputed$label_site_label) |
    is.na(df_imputed$label_gender)), 
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

## 8. Linearity assumption check
library(car)

df_tested <- df_finaal %>% 
  select(label_site_label, label_age, label_gender, HOXB6_1)

df_tested$label_age <- df_tested$label_age + abs(min(df_tested$label_age)) + 0.01
df_tested$HOXB6_1 <- df_tested$HOXB6_1 + abs(min(df_tested$HOXB6_1)) + 0.01

model_tested <- glm(
  label_site_label ~ label_age + HOXB6_1 +
    label_age:log(label_age) + HOXB6_1:log(HOXB6_1),
  family = binomial,
  data = df_tested
)
summary(model_tested)

## 9. Logistic Regression Models

## 9.1 Individual (No interactions)
model_1 <- glm(
  label_site_label ~ label_age + label_gender + HOXB6_1,
  data = df_finaal,
  family = binomial
)
summary(model_1)
vif(model_1)

## 9.2 Without_gender (Not significant)
model_2 <- glm(
  label_site_label ~ label_age  + HOXB6_1,
  data = df_finaal,
  family = binomial
)
summary(model_2)
vif(model_2)

## 9.3 All interactions
model_interactie <- glm(
  label_site_label ~ label_age * HOXB6_1 + label_age + HOXB6_1 + label_gender + label_gender*label_age + label_gender*HOXB6_1,
  data = df_finaal,
  family = binomial
)
summary(model_interactie)
vif(model_interactie)

## 10. Visualisations

## 10.1 Violin Plot HOXB6_1 methylation
ggplot(df_finaal, aes(x = factor(label_site_label), y = HOXB6_1, fill = factor(label_site_label))) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black") +
  scale_x_discrete(labels = c("0" = "Primary", "1" = "Metastatic")) +
  scale_fill_manual(values = c("0" = "lightcoral", "1" = "turquoise"),
                    labels = c("Primary", "Metastatic"),
                    name = 'Metastasis status') +
  
  labs(
    title = "HOXB6_1 methylation per metastasis status",
    x = "Metastasis status",
    y = "HOXB6_1 methylation"
  ) +
  theme_minimal()

## 10.2 Violin Plot Age
ggplot(df_finaal, aes(x = factor(label_site_label), y = label_age, fill = factor(label_site_label))) +
geom_violin(alpha = 0.6) +
geom_boxplot(width = 0.1, color = "black") +
scale_x_discrete(labels = c("0" = "Primary", "1" = "Metastatic")) +
scale_fill_manual(values = c("0" = "lightcoral", "1" = "turquoise"),
labels = c("Primary", "Metastatic"),
name = 'Metastasis status') +
labs(
title = "Age per metastasis status",
x = "Metastasis status",
y = "Age"
) +
theme_minimal()

## 10.3 Heatmap
plotdata <- expand.grid(
  HOXB6_1 = seq(min(df_finaal$HOXB6_1), max(df_finaal$HOXB6_1), length.out = 100),
  label_age = seq(min(df_finaal$label_age), max(df_finaal$label_age), length.out = 100)
)

plotdata$pred_prob <- predict(
  glm(label_site_label ~ label_age + HOXB6_1, 
      family = binomial, data = df_finaal),
  newdata = newdata, type = "response"
)

ggplot(plotdata, aes(x = HOXB6_1, y = label_age, z = pred_prob)) +
  geom_tile(aes(fill = pred_prob)) +
  geom_contour(color = "white", alpha = 0.7) +
  scale_fill_viridis_c(name = "Predicted probability of metastasis") +
  labs(
    x = "HOXB6 methylation (beta value)",
    y = "Age",
    title = "Predicted probability of metastasis based on HOXB6 and age"
  ) +
  theme_minimal()
