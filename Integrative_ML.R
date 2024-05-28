library(limma)
library(dplyr)
library(edgeR)
library(randomForest)
library(pROC)
library(ggplot2)

# Read data
cluster_results <- read.csv("cluster_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
response <- read.csv("Womac_score_pain_function.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
response <- response %>% filter(ID %in% cluster_results$ID[cluster_results$Cluster_Labels == 0])
response <- response[complete.cases(response$Womac.Pain), ]
###Strong vs Mild+No
response$Womac.Function_up <- ifelse(response$Womac.Pain == 3, 0, 1)

# Process metabolites data
process_data <- function(file_path, response) {
  count_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, sep = ",")
  data <- count_data %>% filter(ID %in% response$ID) %>% select(-ID)
  return(data)
}

# Differential expression analysis
diff_expression_analysis <- function(data, response) {
  design <- model.matrix(~ factor(response$Womac.Function_up))
  colnames(design)[2] <- "NewName"
  y <- DGEList(counts = t(data), group = response$Womac.Function_up)
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = TRUE)
  fit <- lmFit(v, design)
  contrast <- makeContrasts(NewName = "NewName", levels = design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "NewName", number = Inf)
  return(results)
}

# Process and analyze data
metabolites_data <- process_data("/Users/divvi/Documents/PMH/Kapoor lab/metabolites_raw.csv", response)
results_meta <- diff_expression_analysis(metabolites_data, response)
significant_results_meta <- results_meta[abs(results_meta$logFC) >= 0.138, ]

# Random forest classification and ROC
random_forest_analysis <- function(filtered_expression_values, response) {
  rf_classifier <- randomForest(response$Womac.Function_up ~ ., data = filtered_expression_values, ntree = 500, mtry = 10, importance = TRUE)
  predictions <- predict(rf_classifier, data = filtered_expression_values, type = "prob")
  roc_curve <- roc(response$Womac.Function_up, predictions[,1])
  auc_value <- auc(roc_curve)
  return(list(rf_classifier = rf_classifier, roc_curve = roc_curve, auc_value = auc_value))
}

# Filter significant results
filter_significant_results <- function(results, logFC_threshold) {
  return(results[abs(results$logFC) >= logFC_threshold, ])
}

# Get top 20 variables
get_top_variables <- function(rf_classifier, significant_results) {
  importance_values <- as.data.frame(importance(rf_classifier))
  importance_values_df <- as.data.frame(cbind(rownames(importance_values), importance_values$MeanDecreaseGini))
  colnames(importance_values_df) <- c("Variable", "Importance")
  importance_values_df <- importance_values_df[order(importance_values_df$Importance, decreasing = TRUE), ]
  top_variables <- as.data.frame(importance_values_df[1:min(20, nrow(importance_values_df)), ])
  colnames(top_variables) <- c("Variables", "Importance")
  matched_indices <- match(top_variables$Variables, rownames(significant_results))
  top_variables$logFC <- significant_results$logFC[matched_indices]
  top_variables$P.Value <- significant_results$P.Value[matched_indices]
  min_value <- min(as.numeric(top_variables$Importance), na.rm = TRUE)
  max_value <- max(as.numeric(top_variables$Importance), na.rm = TRUE)
  top_variables$Normalized_Importance <- ((as.numeric(top_variables$Importance) - min_value) / (max_value - min_value)) * 100
  return(top_variables)
}

# Process metabolites data
metabolite_data <- process_data("metabolites.csv", response)
results_metabolite <- diff_expression_analysis(metabolite_data, response)
significant_results_metabolite <- filter_significant_results(results_metabolite, 0.138)
expression_values_metabolite <- as.data.frame(t(voom(DGEList(counts = t(metabolite_data), group = response$Womac.Function_up), model.matrix(~ factor(response$Womac.Function_up)), plot = FALSE)$E))
filtered_expression_values_metabolite <- expression_values_metabolite[, rownames(significant_results_metabolite)]
rf_metabolite <- random_forest_analysis(filtered_expression_values_metabolite, response)
top_20_variables_metabolite <- get_top_variables(rf_metabolite$rf_classifier, significant_results_metabolite)
predictions_metabolite <- predict(rf_metabolite, data = filtered_expression_values_metabolite, type = "prob")
roc_curve_metabolite <- roc(response, predictions_metabolite[, 1])
# Calculate AUC
auc_value_metabolite <- auc(roc_curve_metabolite)

# Process RNA plasma data
rna_plasma_data <- process_data("RNA_plasma.csv", response)
results_rna_plasma <- diff_expression_analysis(rna_plasma_data, response)
significant_results_rna_plasma <- filter_significant_results(results_rna_plasma, 0.585)
expression_values_rna_plasma <- as.data.frame(t(voom(DGEList(counts = t(rna_plasma_data), group = response$Womac.Function_up), model.matrix(~ factor(response$Womac.Function_up)), plot = FALSE)$E))
filtered_expression_values_rna_plasma <- expression_values_rna_plasma[, rownames(significant_results_rna_plasma)]
rf_rna_plasma <- random_forest_analysis(filtered_expression_values_rna_plasma, response)
top_20_variables_rna_plasma <- get_top_variables(rf_rna_plasma$rf_classifier, significant_results_rna_plasma)
predictions_RNA_plasma <- predict(rf_rna_plasma,data=filtered_expression_values_rna_plasma,type = "prob")
roc_curve_rnaplasma <- roc(response, predictions_RNA_plasma[,1])
# Calculate AUC
auc_value_rnaplasma <- auc(roc_curve_rnaplasma)


# Process RNA synovial data
rna_syn_data <- process_data("RNA_synovial.csv", response)
results_rna_syn <- diff_expression_analysis(rna_syn_data, response)
significant_results_rna_syn <- filter_significant_results(results_rna_syn, 0.585)
expression_values_rna_syn <- as.data.frame(t(voom(DGEList(counts = t(rna_syn_data), group = response$Womac.Function_up), model.matrix(~ factor(response$Womac.Function_up)), plot = FALSE)$E))
filtered_expression_values_rna_syn <- expression_values_rna_syn[, rownames(significant_results_rna_syn)]
rf_rna_syn <- random_forest_analysis(filtered_expression_values_rna_syn, response)
top_20_variables_rna_syn <- get_top_variables(rf_rna_syn$rf_classifier, significant_results_rna_syn)
predictions_RNA_syn <- predict(rf_rna_syn,data=filtered_expression_values_rna_syn,type = "prob")
roc_curve_rnasyn <- roc(response, predictions_RNA_syn[,1])
# Calculate AUC
auc_value_rnasyn <- auc(roc_curve_rnasyn)

# Process RNA urine data
rna_urine_data <- process_data("RNA_urine.csv", response)
results_rna_urine <- diff_expression_analysis(rna_urine_data, response)
significant_results_rna_urine <- filter_significant_results(results_rna_urine, 0.585)
expression_values_rna_urine <- as.data.frame(t(voom(DGEList(counts = t(rna_urine_data), group = response$Womac.Function_up), model.matrix(~ factor(response$Womac.Function_up)), plot = FALSE)$E))
filtered_expression_values_rna_urine <- expression_values_rna_urine[, rownames(significant_results_rna_urine)]
rf_rna_urine <- random_forest_analysis(filtered_expression_values_rna_urine, response)
top_20_variables_rna_urine <- get_top_variables(rf_rna_urine$rf_classifier, significant_results_rna_urine)
predictions_RNA_urine <- predict(rf_rna_urine,data=filtered_expression_values_rna_urine,type = "prob")
roc_curve_rnaurine <- roc(response, predictions_RNA_urine[,1])
# Calculate AUC
auc_value_rnaurine <- auc(roc_curve_rnaurine)


###Clinical Data

data <- read.csv("Clinical_data.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
data <- data[, 1:7]

# Convert specific columns to factors
data$anxiety <- factor(data$anxiety)
data$depression <- factor(data$depression)
data$pain_detect_categories <- factor(data$pain_detect_categories)

# Impute missing data
imp <- mice(data, maxit = 10, m = 1, seed = 1)
imputed_data <- complete(imp, 1)
imputed_data$sex <- factor(imputed_data$sex)

rf_clinical <- random_forest_analysis(imputed_data, response)
top_20_variables_rf_clinical <- get_top_variables(rf_clinicale$rf_classifier, imputed_data)
predictions_clinical <- predict(rf_clinical,data=imputed_data,type = "prob")
roc_curve_clinical <- roc(response, predictions_clinical[,1])
# Calculate AUC
auc_value_clinical <- auc(roc_curve_clinical)


# Combine top variables from all analyses
top_20_variables_rna_plasma$Variables <- paste0(top_20_variables_rna_plasma$Variables, "_pla")
top_20_variables_rna_syn$Variables <- paste0(top_20_variables_rna_syn$Variables, "_syn")
top_20_variables_rna_urine$Variables <- paste0(top_20_variables_rna_urine$Variables, "_uri")

final <- rbind(top_20_variables_metabolites, top_20_variables_rna_plasma, top_20_variables_rna_syn, top_20_variables_rna_urine)
final <- arrange(final, desc(Normalized_Importance))

# Extract the top 20 rows
top_20

 #Combine predictions into a single feature matrix (assuming these variables are defined elsewhere)
 combined_features <- as.data.frame(cbind(predictions_meta[, 1], predictions_RNA_plasma[, 1], predictions_RNA_syn[, 1], predictions_RNA_uri[, 1], predictions_clinical[,1]))

 #Fit the meta-classifier
 library(e1071) # Loading the library for Naive Bayes classifier

 # Fit Naive Bayes classifier
 meta_classifier <- naiveBayes(response$Womac.Function_up ~ ., data = combined_features)
 meta_predictions <- predict(meta_classifier, combined_features)
 roc_curve_combined <- roc(response, meta_predictions)
 auc_value <- auc(roc_curve_combined)

# Plot ROC curves
auc_values <- c(auc_value_metabolite, auc_value_rnaplasma, auc_value_rnasyn, auc_rnaurine, auc_value_clinical, auc_combined) 
curve_names <- c("Metabolite (AUC =", "miRNA plasma (AUC =", "miRNA synovial (AUC =", "miRNA urine (AUC =", "Clinical (AUC =", "Integrated (AUC =")

plot(1 - roc_curve_metabolite[["specificities"]], roc_curve_metabolite[["sensitivities"]], col = "purple", 
     main = "ROC Curves", type = "l", lwd = 1, 
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)")

lines(1 - roc_curve_rnaplasma[["specificities"]], roc_curve_rnaplasma[["sensitivities"]], col = "blue", lwd = 1)
lines(1 - roc_curve_rnasyn[["specificities"]], roc_curve_rnasyn[["sensitivities"]], col = "green", lwd = 1)
lines(1 - roc_curve_rnaurine[["specificities"]], roc_curve_rnaurine[["sensitivities"]], col = "darkorange3", lwd = 1)
lines(1 - roc_curve_clinical[["specificities"]], roc_curve_clinical[["sensitivities"]], col = "black", lwd = 1)
lines(1 - roc_curve_combined[["specificities"]], roc_curve_combined[["sensitivities"]], col = "red", lwd = 3)

abline(a = 0, b = 1, col = "black", lty = "dashed")
grid()

# Add legend
legend("bottomright", legend = paste(curve_names, round(auc_values, 3), ")", sep = ""), 
       col = c("purple", "blue", "green", "darkorange3", "black", "red"), lwd = c(1, 1, 1, 1, 1, 3),
       xjust = 0.5, yjust = 0.5, box.lwd = 0)


