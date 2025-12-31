
#---------Pearson correlation of 10 features---------

#Data loads
data <- read.table("Dual-sgRNA_design_input.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


library(dplyr)

# min 
data_min <- data %>%
  mutate(
    Microhomology = pmin(sg1_MH, sg2_MH),
    GC_content     = pmin(sg1_GC, sg2_GC),
    RuleSet2       = pmin(sg1_Doench, sg2_Doench),
    DeepSpCas9     = pmin(sg1_Deep, sg2_Deep),
    CRISPRscan     = pmin(sg1_Mor, sg2_Mor),
    RuleSet3       = pmin(sg1_RuleSet3, sg2_RuleSet3),
    CRISPRon       = pmin(sg1_CRISPRon, sg2_CRISPRon)
  )


# correlation for features
selected_vars <- data_min %>%
  select(Efficiency, Deletion_bp, Strand,
         Microhomology, GC_content, RuleSet2, DeepSpCas9,
         CRISPRscan, RuleSet3, CRISPRon)

# cor cal 
cor_mat <- cor(selected_vars %>% select(-Strand), use = "complete.obs")

# features seq reset
var_order <- c("Efficiency", "RuleSet2", "RuleSet3", "CRISPRscan",
               "DeepSpCas9", "CRISPRon", "Deletion_bp", "Microhomology", "GC_content")
cor_mat_reordered <- cor_mat[var_order, var_order]


# print
eff_cor <- cor_mat["Efficiency", ]
eff_cor <- eff_cor[names(eff_cor) != "Efficiency"]
print(round(eff_cor, 3))



#---------Spearman correlation-(DeepSpCas9 featrue)---------


spearman_r <- cor(
  data_min$DeepSpCas9,
  data_min$Efficiency,
  method = "spearman",
  use = "complete.obs"
)
spearman_rsq <- spearman_r^2

# print
cat("Spearman rho:", round(spearman_r, 3), "\n")
cat("Spearman rho^2:", round(spearman_rsq, 3), "\n")



#--------- Distribution  analysis of 10 features --------------

library(dplyr)
library(tidyr)

vars <- c(
  "Efficiency", "RuleSet2", "RuleSet3", "CRISPRscan",
  "DeepSpCas9", "CRISPRon", "Deletion_bp",
  "Microhomology", "GC_content"
)

# long form 
long_dist <- data_min %>%
  select(all_of(vars)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")


# median cal 
medians <- long_dist %>%
  group_by(Variable) %>%
  summarise(Median = median(Value, na.rm = TRUE))


print(mutate(medians, Median = round(Median, 2)))

#Efficiency median
medians$Median[medians$Variable == "Efficiency"]




#---------Comparative statistical analysis of 10 features --------------


library(dplyr)


# group setting
data_min$eff_group <- ifelse(data_min$Efficiency >= median(data_min$Efficiency, na.rm = TRUE),
                             "High", "Low")
data_min$eff_group <- factor(data_min$eff_group, levels = c("Low", "High"))

# feature list 
vars <- c("RuleSet2", "RuleSet3", "DeepSpCas9", "CRISPRon",
          "Deletion_bp", "CRISPRscan", "Microhomology", "GC_content")


# total t-test 
t_test_table <- lapply(vars, function(v) {
  t <- t.test(data_min[[v]] ~ data_min$eff_group)
  data.frame(
    Variable = v,
    Mean_Low = mean(data_min[[v]][data_min$eff_group == "Low"], na.rm = TRUE),
    Mean_High = mean(data_min[[v]][data_min$eff_group == "High"], na.rm = TRUE),
    p_value = signif(t$p.value, 4),
    Method = t$method
  )
}) %>% bind_rows()

print(t_test_table)


#---------sgRNA strand orientation --------------

#  Wilcoxon 
wilcox_result <- pairwise.wilcox.test(
  x = data_min$Efficiency,
  g = data_min$Strand,
  p.adjust.method = "BH"
)

# result collection 
wilcox_pvals <- as.data.frame(as.table(wilcox_result$p.value))
colnames(wilcox_pvals) <- c("group1", "group2", "p_adj")
wilcox_pvals <- na.omit(wilcox_pvals)
wilcox_pvals$p_label <- sprintf("p = %.3f", wilcox_pvals$p_adj)

# print
print(wilcox_result)
print(wilcox_pvals)

# N per group
print(table(data_min$Strand))



#---------Design feature based modeling --------------
#---------Design feature based modeling --------------

library(dplyr)
library(randomForest)
library(pROC)

median_eff <- median(data$Efficiency, na.rm = TRUE)

model_data <- data %>%
  mutate(
    label = factor(ifelse(Efficiency >= median_eff, 1, 0)),
    DeepSpCas9 = pmin(sg1_Deep, sg2_Deep, na.rm = TRUE)
  ) %>%
  select(label, DeepSpCas9) %>%
  filter(complete.cases(.))

# Design_GLM
glm_model_design <- glm(label ~ DeepSpCas9, data = model_data, family = binomial())
glm_probs_design <- predict(glm_model_design, type = "response")
glm_roc_design <- roc(model_data$label, glm_probs_design, quiet = TRUE)

cat("Design_GLM AUC:", round(auc(glm_roc_design), 3), "\n")

# Design_RF
set.seed(20250618)
rf_model_design <- randomForest(label ~ DeepSpCas9, data = model_data, importance = TRUE, ntree = 500)
rf_probs_design <- predict(rf_model_design, type = "prob")[, 2]
rf_roc_design <- roc(model_data$label, rf_probs_design, quiet = TRUE)

cat("Design_RF AUC:", round(auc(rf_roc_design), 3), "\n")


#---------K-mer based modeling --------------
#---------K-mer based modeling --------------

library(dplyr)
library(randomForest)
library(pROC)

# count matrix read 
data <- read.table("Kmer_countmatrix.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

model_data <- data %>%
  mutate(label = factor(Label)) %>%
  select(-`No.`, -Label, -Class) %>%
  filter(complete.cases(.))

# KMER_GLM
glm_model_kmer <- glm(label ~ ., data = model_data, family = binomial())
glm_probs_kmer <- predict(glm_model_kmer, type = "response")
glm_roc_kmer <- roc(model_data$label, glm_probs_kmer, quiet = TRUE)

cat("KMER_GLM AUC:", round(auc(glm_roc_kmer), 3), "\n")

# KMER_RF
set.seed(20250618)
rf_model_kmer <- randomForest(label ~ ., data = model_data, importance = TRUE, ntree = 500)
rf_probs_kmer <- predict(rf_model_kmer, type = "prob")[, 2]
rf_roc_kmer <- roc(model_data$label, rf_probs_kmer, quiet = TRUE)

cat("KMER_RF AUC:", round(auc(rf_roc_kmer), 3), "\n")


#---------KMER GLM & RF importance --------------

library(ggplot2)
library(tibble)

# KMER GLM
glm_imp_df <- summary(glm_model_kmer)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  filter(Feature != "(Intercept)") %>%
  mutate(Importance = abs(Estimate)) %>%
  arrange(desc(Importance)) %>%
  slice(1:10)  

print(glm_imp_df)


# KMER RF
rf_imp_df <- importance(rf_model_kmer, type = 2) %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice(1:10)

print(rf_imp_df)



#---------frequencies of statistically significant K-mer(2mers) --------------

library(dplyr)
library(tidyr)

# read count matrix
df <- read.table("Kmer_countmatrix.txt",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# label 
df <- df %>%
  mutate(label = as.character(Label))

# k-mer (k=2)
kmer_cols <- grep("^[ACGT]{2}$", colnames(df), value = TRUE)

# Group-wise k-mer counts
kmer_group_sum <- df %>%
  select(label, all_of(kmer_cols)) %>%
  group_by(label) %>%
  summarise(across(all_of(kmer_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = -label, names_to = "kmer", values_to = "count") %>%
  pivot_wider(names_from = label, values_from = count,
              names_prefix = "group_", values_fill = 0)

# Fisher’s exact test with Benjamini–Hochberg correction

tot1 <- sum(kmer_group_sum$group_1, na.rm = TRUE)
tot0 <- sum(kmer_group_sum$group_0, na.rm = TRUE)

kmer_stats <- kmer_group_sum %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(
      matrix(c(group_1, group_0,
               tot1 - group_1, tot0 - group_0), nrow = 2)
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(fdr)

# Important k-mer

sig_kmer_table <- kmer_stats %>%
  filter(fdr < 0.05) %>%
  mutate(
    FDR = signif(fdr, 3),
    p_value = signif(p_value, 3)
  ) %>%
  select(
    kmer,
    Low  = group_0,
    High = group_1,
    p_value,
    FDR
  )

print(sig_kmer_table)



#---------Venn diagram of important k-mer in KMER GLM & RF --------------

library(VennDiagram)
library(grid)

# set setup 
top_glm_kmers <- unique(glm_imp_df$Feature)
top_rf_kmers  <- unique(rf_imp_df$Feature)
sig_fdr_kmers <- unique(sig_kmer_table$kmer)

# venn production
venn.plot <- venn.diagram(
  x = list(
    `GLM top10`  = top_glm_kmers,
    `RF top10`   = top_rf_kmers,
    `FDR < 0.05` = sig_fdr_kmers
  ),
  filename = NULL,
  fill = c("#002366", "#4CAF50", "#e31a1c"),
  alpha = 0.5,
  cex = 1.3,
  cat.cex = 1.3,
  cat.fontface = "bold",
  margin = 0.1,
  main = NULL
)

tiff("kmer_venn_diagram.tiff", width = 4, height = 4, units = "in", res = 300, compression = "lzw")
grid.newpage()
grid.draw(venn.plot)
dev.off()

# intersections 
cat("GLM ∩ RF:\n"); print(intersect(top_glm_kmers, top_rf_kmers))
cat("\nGLM ∩ FDR:\n"); print(intersect(top_glm_kmers, sig_fdr_kmers))
cat("\nRF ∩ FDR:\n");  print(intersect(top_rf_kmers, sig_fdr_kmers))
cat("\nGLM ∩ RF ∩ FDR:\n"); print(Reduce(intersect, list(top_glm_kmers, top_rf_kmers, sig_fdr_kmers)))



#---------Combined modeling --------------
#---------Combined modeling --------------

combined_glm_prob <- (glm_probs_design + glm_probs_kmer) / 2
combined_rf_prob <- (rf_probs_design + rf_probs_kmer) / 2

combined_glm_roc <- roc(model_data$label, combined_glm_prob)
cat("Combined GLM AUC:", round(auc(combined_glm_roc), 3), "\n")

combined_rf_roc <- roc(model_data$label, combined_rf_prob)
cat("Combined RF AUC:", round(auc(combined_rf_roc), 3), "\n")



#---------Combined GLM & RF importance --------------

library(dplyr)
library(tibble)

#Design GLM importance
glm_imp_design <- summary(glm_model_design)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  filter(Feature != "(Intercept)") %>%
  mutate(Importance_design = abs(Estimate)) %>%
  select(Feature, Importance_design)

# K-mer GLM importance
glm_imp_kmer <- summary(glm_model_kmer)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  filter(Feature != "(Intercept)") %>%
  mutate(Importance_kmer = abs(Estimate)) %>%
  select(Feature, Importance_kmer)

# Combined importance
combined_glm_imp <- full_join(
  glm_imp_design,
  glm_imp_kmer,
  by = "Feature"
) %>%
  mutate(
    CombinedImportance = rowMeans(
      select(., Importance_design, Importance_kmer),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(CombinedImportance)) %>%
  slice(1:10)

# Print 
cat("\nTop 10 features by Combined GLM importance:\n")
print(combined_glm_imp)


# Design RF importance
rf_imp_design <- importance(rf_model_design, type = 2) %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  rename(Importance_design = MeanDecreaseGini) %>%
  select(Feature, Importance_design)

# K-mer RF importance
rf_imp_kmer <- importance(rf_model_kmer, type = 2) %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  rename(Importance_kmer = MeanDecreaseGini) %>%
  select(Feature, Importance_kmer)

# Combined importance
combined_rf_imp <- full_join(
  rf_imp_design,
  rf_imp_kmer,
  by = "Feature"
) %>%
  mutate(
    CombinedImportance = rowMeans(
      select(., Importance_design, Importance_kmer),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(CombinedImportance)) %>%
  slice(1:10)

# Print 
cat("\nTop 10 features by Combined RF importance:\n")
print(combined_rf_imp)



#---------Hybrid modeling --------------
#---------Hybrid modeling --------------

# hybrid1: Design GLM + KMER RF
hybrid1_prob <- (glm_probs_design + rf_probs_kmer) / 2
hybrid1_roc  <- roc(model_data$label, hybrid1_prob)

cat("Hybrid1 (Design GLM + KMER RF) AUC:", round(auc(hybrid1_roc), 3), "\n")

# hybrid2: Design RF + KMER GLM
hybrid2_prob <- (rf_probs_design + glm_probs_kmer) / 2
hybrid2_roc  <- roc(model_data$label, hybrid2_prob)

cat("Hybrid2 (Design RF + KMER GLM) AUC:", round(auc(hybrid2_roc), 3), "\n")




#-----------10-fold cross-validation across different models-----------


library(pROC)
library(dplyr)

# label vector
label_vec <- factor(model_data$label, levels = c(0, 1)) 

# prediction probablility
cv_data <- data.frame(
  label           = label_vec,
  Design_GLM      = glm_probs_design,
  Design_RF       = rf_probs_design,
  KMER_GLM        = glm_probs_kmer,
  KMER_RF         = rf_probs_kmer,
  Combined_GLM    = combined_glm_prob,
  Combined_RF     = combined_rf_prob,
  Hybrid1         = hybrid1_prob,
  Hybrid2         = hybrid2_prob
)

# justice for model name
model_names <- colnames(cv_data)[-1] 

# repeat and fold N selection 
set.seed(20250625)
k <- 10
repeats <- 10

# result
auc_results <- lapply(model_names, function(x) numeric(k * repeats))
names(auc_results) <- model_names

# repeat kfold operation
counter <- 1
for (r in 1:repeats) {
  repeat {
    folds <- sample(rep(1:k, length.out = nrow(cv_data)))
    valid_folds <- all(sapply(1:k, function(i) length(unique(cv_data$label[folds == i])) == 2))
    if (valid_folds) break
  }
  
  for (i in 1:k) {
    test_idx <- which(folds == i)
    test_label <- cv_data$label[test_idx]
    
    for (model in model_names) {
      test_pred <- cv_data[[model]][test_idx]
      auc_val <- auc(test_label, test_pred)
      auc_results[[model]][counter] <- auc_val
    }
    
    counter <- counter + 1
  }
}

# mean of AUC and SD 
cv_auc_df <- data.frame(
  Model    = model_names,
  Mean_AUC = sapply(auc_results, mean),
  SD       = sapply(auc_results, sd)
)

print(cv_auc_df)



# t-test

# AUC long form transformation 
auc_long <- bind_rows(lapply(names(auc_results), function(m) {
  data.frame(Model = m, AUC = auc_results[[m]])
}))


comparisons <- list(
  c("Combined_GLM", "KMER_GLM"),
  c("Combined_GLM", "Hybrid1"),
  c("Combined_GLM", "Hybrid2")
)

# Helper to extract vectors
get_auc <- function(m) auc_long$AUC[auc_long$Model == m]

# Paired t-tests
paired_df <- do.call(rbind, lapply(comparisons, function(cp) {
  x <- get_auc(cp[1])
  y <- get_auc(cp[2])
  stopifnot(length(x) == length(y))
  
  tt <- t.test(x, y, paired = TRUE)
  data.frame(
    Comparison = paste(cp[1], "vs", cp[2]),
    Mean_diff = mean(x - y),
    SD_diff   = sd(x - y),
    p_value   = tt$p.value,
    stringsAsFactors = FALSE
  )
}))

paired_df$p_adj_BH <- p.adjust(paired_df$p_value, method = "BH")

cat("\nPaired t-tests:\n")
print(transform(paired_df,
                Mean_diff = round(Mean_diff, 3),
                SD_diff   = round(SD_diff, 3),
                p_value   = signif(p_value, 4),
                p_adj_BH  = signif(p_adj_BH, 4)))





#--------cross validated pearson correlation 

# Pearson 
all_preds_list <- lapply(model_names, function(x) numeric(nrow(cv_data)))
names(all_preds_list) <- model_names
all_labels <- as.numeric(cv_data$label)  

# Kfold set
set.seed(20250625)
k <- 10
repeats <- 10

# Operation for kfold
for (r in 1:repeats) {
  repeat {
    folds <- sample(rep(1:k, length.out = nrow(cv_data)))
    valid_folds <- all(sapply(1:k, function(i) length(unique(cv_data$label[folds == i])) == 2))
    if (valid_folds) break
  }
  
  for (i in 1:k) {
    test_idx <- which(folds == i)
    
    for (model in model_names) {
      test_pred <- cv_data[[model]][test_idx]
      all_preds_list[[model]][test_idx] <- test_pred
    }
  }
}

# Pearson R calculation 
pearson_r_results <- sapply(model_names, function(model) {
  cor(all_preds_list[[model]], all_labels, method = "pearson")
})

# Result
pearson_r_df <- data.frame(
  Model     = model_names,
  Pearson_R = round(pearson_r_results, 4)
)

print(pearson_r_df)


#--------cross validated spearman correlation 

# Spearman
all_preds_list_spearman <- lapply(model_names, function(x) numeric(nrow(cv_data)))
names(all_preds_list_spearman) <- model_names
all_labels <- as.numeric(cv_data$label)  

# kfold set
set.seed(20250625)
k <- 10
repeats <- 10

# Operation for kfold
for (r in 1:repeats) {
  repeat {
    folds <- sample(rep(1:k, length.out = nrow(cv_data)))
    valid_folds <- all(sapply(1:k, function(i) length(unique(cv_data$label[folds == i])) == 2))
    if (valid_folds) break
  }
  
  for (i in 1:k) {
    test_idx <- which(folds == i)
    
    for (model in model_names) {
      test_pred <- cv_data[[model]][test_idx]
      all_preds_list_spearman[[model]][test_idx] <- test_pred
    }
  }
}

# spearman R calculation 
spearman_r_results <- sapply(model_names, function(model) {
  cor(all_preds_list_spearman[[model]], all_labels, method = "spearman")
})

# result
spearman_r_df <- data.frame(
  Model      = model_names,
  Spearman_R = round(spearman_r_results, 4)
)

print(spearman_r_df)




#----------confusion matrix----------

library(caret)

threshold <- 0.5

# labeling
true_class <- factor(ifelse(model_data$label == 1, "High", "Low"),
                     levels = c("High", "Low"))

pred_class <- factor(ifelse(combined_glm_prob >= threshold, "High", "Low"),
                     levels = c("High", "Low"))

cm <- confusionMatrix(pred_class, true_class, positive = "High")

cat("\n Combined GLM Confusion Matrix \n")

# reorder
cm_table_ordered <- cm$table[c("Low", "High"), c("High", "Low")]
print(cm_table_ordered)

# print
cat("\n[Overall]\n"); print(cm$overall)
cat("\n[ByClass]\n"); print(cm$byClass)
cm$overall["Accuracy"]



# ---------------- Precision-Recall------------

library(PRROC)

labels01 <- as.numeric(as.character(model_data$label))

pr <- pr.curve(
  scores.class0 = combined_glm_prob[labels01 == 1],
  scores.class1 = combined_glm_prob[labels01 == 0],
  curve = TRUE
)

cat("\n[PR Curve] Combined_GLM\n")
cat("PR-AUC (integral) =", sprintf("%.4f", pr$auc.integral), "\n")
cat("PR-AUC (Davis-Goadrich) =", sprintf("%.4f", pr$auc.davis.goadrich), "\n")
cat("Prevalence (baseline precision) =", sprintf("%.4f", mean(labels01 == 1)), "\n")

curve <- pr$curve
recall <- curve[,1]; precision <- curve[,2]; thr <- curve[,3]

den <- precision + recall
f1 <- ifelse(den == 0, NA_real_, 2 * precision * recall / den)
idx <- which.max(f1)

cat("Best F1 =", sprintf("%.4f", f1[idx]),
    "at threshold =", sprintf("%.4f", thr[idx]),
    " (Precision =", sprintf("%.4f", precision[idx]),
    ", Recall =", sprintf("%.4f", recall[idx]), ")\n")

# console-friendly points
grid <- seq(0, 1, 0.1)
pick <- sapply(grid, function(g) which.min(abs(recall - g)))
print(data.frame(Recall=recall[pick], Precision=precision[pick], Threshold=thr[pick]),
      row.names = FALSE)



#-------Delong test-------

library(pROC)

# justice of model list
roc_list <- list(
  Combined_GLM = combined_glm_roc,
  Hybrid1 = hybrid1_roc,
  Hybrid2 = hybrid2_roc,
  KMER_GLM = glm_roc_kmer
)

# model pair prodcution 
model_names <- names(roc_list)
combs <- combn(model_names, 2, simplify = FALSE)

# test
delong_results <- lapply(combs, function(pair) {
  res <- roc.test(roc_list[[pair[1]]], roc_list[[pair[2]]], method = "delong")
  data.frame(
    Model1 = pair[1],
    Model2 = pair[2],
    Z = as.numeric(res$statistic),
    P_value = as.numeric(res$p.value),
    stringsAsFactors = FALSE
  )
})

# result 
delong_df <- do.call(rbind, delong_results)
print(delong_df)


# confirm
sapply(roc_list, function(r) length(r$response))
sapply(roc_list, function(r) length(r$predictor))
sapply(roc_list, function(r) identical(names(r$predictor), names(roc_list$Combined_GLM$predictor)))



#---- Brier Score  ----------

labels01 <- as.numeric(as.character(model_data$label))

brier_score_combined <- mean((combined_glm_prob - labels01)^2, na.rm = TRUE)

cat("\n[Brier Score] Combined_GLM\n")
cat("Brier Score =", sprintf("%.4f", brier_score_combined), "\n")


#------Hosmer–Lemeshow---------

library(ResourceSelection)

labels01 <- as.numeric(as.character(model_data$label))

hl_combined <- hoslem.test(
  labels01,
  combined_glm_prob,
  g = 10
)

cat("\n[Hosmer–Lemeshow Test] Combined_GLM\n")
cat("Chi-square =", round(hl_combined$statistic, 4), "\n")
cat("df         =", hl_combined$parameter, "\n")
cat("p-value    =", round(hl_combined$p.value, 4), "\n")




#---------Y-scrambling test --------------

library(dplyr)
library(pROC)

# Prepare design-feature dataset (label + DeepSpCas9)
# The design input file is reloaded here to avoid overwriting

data_design <- read.table("Dual-sgRNA_design_input.txt",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

median_eff <- median(data_design$Efficiency, na.rm = TRUE)

model_vars <- data_design %>%
  mutate(
    label = ifelse(Efficiency >= median_eff, 1, 0),
    DeepSpCas9 = pmin(sg1_Deep, sg2_Deep, na.rm = TRUE)
  ) %>%
  select(label, DeepSpCas9) %>%
  filter(complete.cases(.))

# Prepare k-mer feature dataset (label + k-mer features)

data_kmer <- read.table("Kmer_countmatrix.txt",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

model_data <- data_kmer %>%
  mutate(label = as.integer(as.character(Label))) %>%   # Label이 "0"/"1"이라고 가정
  select(-`No.`, -Label, -Class) %>%
  filter(complete.cases(.))

# sanity check
stopifnot(nrow(model_vars) == nrow(model_data))

# Construct combined dataset for ensemble modeling
X <- model_vars %>%
  select(label, DeepSpCas9) %>%
  bind_cols(model_data %>% select(-label))

X$label <- factor(X$label, levels = c(0,1))

# Y-scrambling

make_folds <- function(y, k = 5, seed = 1) {
  set.seed(seed)
  y <- factor(y, levels = c(0,1))
  idx0 <- which(y == 0); idx1 <- which(y == 1)
  f0 <- split(sample(idx0), rep(1:k, length.out = length(idx0)))
  f1 <- split(sample(idx1), rep(1:k, length.out = length(idx1)))
  lapply(1:k, function(i) c(f0[[i]], f1[[i]]))
}

cv_combined_glm <- function(dat, y_col = "label", k = 5, seed = 1, y_override = NULL) {
  y <- dat[[y_col]]
  if (!is.null(y_override)) y <- y_override
  y <- factor(y, levels = c(0,1))
  
  folds <- make_folds(y, k = k, seed = seed)
  
  oof_design <- rep(NA_real_, nrow(dat))
  oof_kmer   <- rep(NA_real_, nrow(dat))
  
  design_cols <- c("DeepSpCas9")
  kmer_cols <- setdiff(colnames(dat), c(y_col, design_cols))
  
  for (i in seq_len(k)) {
    te <- folds[[i]]
    tr <- setdiff(seq_len(nrow(dat)), te)
    
    train <- dat[tr, , drop = FALSE]
    test  <- dat[te, , drop = FALSE]
    train[[y_col]] <- y[tr]
    test[[y_col]]  <- y[te]
    
    m_design <- glm(label ~ ., data = train[, c(y_col, design_cols), drop = FALSE], family = binomial())
    m_kmer   <- glm(label ~ ., data = train[, c(y_col, kmer_cols),   drop = FALSE], family = binomial())
    
    oof_design[te] <- predict(m_design, newdata = test[, c(y_col, design_cols), drop = FALSE], type = "response")
    oof_kmer[te]   <- predict(m_kmer,   newdata = test[, c(y_col, kmer_cols),   drop = FALSE], type = "response")
  }
  
  combined <- (oof_design + oof_kmer) / 2
  auc_val <- as.numeric(auc(roc(y, combined, quiet = TRUE)))
  list(auc = auc_val, combined = combined)
}

y_scramble_test <- function(dat, n_perm = 200, k = 5, seed = 42) {
  base <- cv_combined_glm(dat, k = k, seed = seed)
  base_auc <- base$auc
  
  set.seed(seed)
  perm_auc <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    y_perm <- sample(as.character(dat$label))
    y_perm <- factor(y_perm, levels = c(0,1))
    perm_auc[p] <- cv_combined_glm(dat, k = k, seed = seed + p, y_override = y_perm)$auc
  }
  
  pval <- (sum(perm_auc >= base_auc) + 1) / (n_perm + 1)
  
  cat("Base (CV) combined GLM AUC:", round(base_auc, 3), "\n")
  cat("Permuted mean AUC:", round(mean(perm_auc), 3), " sd:", round(sd(perm_auc), 3), "\n")
  cat("Permutation p-value:", signif(pval, 3), "\n")
  
  invisible(list(base_auc = base_auc, perm_auc = perm_auc, p_value = pval))
}

res_perm <- y_scramble_test(X, n_perm = 200, k = 5, seed = 20250618)



#---------Feature stabiliy (KMER GLM) --------------


feature_stability_kmer <- function(dat, top_k = 50, k = 5, seed = 42) {
  
  y <- factor(dat$label, levels = c(0, 1))
  folds <- make_folds(y, k = k, seed = seed)
  
  # Exclude design feature(s) and keep only k-mer predictors
  design_cols <- c("DeepSpCas9")
  kmer_cols <- setdiff(colnames(dat), c("label", design_cols))
  
  top_by_fold <- vector("list", k)
  
  for (i in seq_len(k)) {
    te <- folds[[i]]
    tr <- setdiff(seq_len(nrow(dat)), te)
    
    train <- dat[tr, c("label", kmer_cols), drop = FALSE]
    train$label <- y[tr]
    
    # Fit k-mer GLM (binomial logistic regression)
    m <- suppressWarnings(glm(label ~ ., data = train, family = binomial()))
    
    b <- coef(m)
    b <- b[names(b) != "(Intercept)"]
    b <- b[!is.na(b)]
    
    if (length(b) == 0) {
      top_by_fold[[i]] <- character(0)
      next
    }
    
    top_features <- names(sort(abs(b), decreasing = TRUE))[1:min(top_k, length(b))]
    top_by_fold[[i]] <- top_features
  }
  
  freq <- sort(table(unlist(top_by_fold)), decreasing = TRUE)
  
  stability <- data.frame(
    kmer = names(freq),
    folds_in_topk = as.integer(freq),
    prop = as.integer(freq) / k,
    stringsAsFactors = FALSE
  )
  
  stability
}

# Run feature stability

stab <- feature_stability_kmer(X, top_k = 50, k = 5, seed = 20250618)

cat("\nTop 30 stable k-mers (by fold frequency):\n")
print(head(stab, 30))


colnames(model_data)

summary(model_data$TT)

