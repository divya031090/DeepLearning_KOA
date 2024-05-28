library(limma) 
library(dplyr)
library(edgeR)

cluster_results<-read.csv("cluster_results.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
cluster_results$Cluster_Labels_up<-ifelse(cluster_results$Cluster_Labels==2,0,1)
##change to other data domains accordingly
count_data <- read.csv("RNA_plasma.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
data <- count_data %>% filter(ID %in% cluster_results$ID)
data<-data[2:length(data)]

design <- model.matrix(~ factor(cluster_results$Cluster_Labels_up))
colnames(design)[colnames(design) == "factor(cluster_results$Cluster_Labels_up)1"] <- "NewName"

y <- DGEList(counts = t(data), group = cluster_results$Cluster_Labels_up)

# Calculate normalizastion factors
y <- calcNormFactors(y)

v <- voom(y, design, plot = TRUE)

# Specify contrasts of interest
contrast <- c(-1, 1)  # Assuming a contrast between two groups

# Fit linear model
fit <- lmFit(v, design)

contrasts <- makeContrasts(contrast, levels = design)

# Empirical Bayes moderation
fit <- eBayes(fit, trend = TRUE)

# Extract differential expression results
results <- topTable(fit, coef = "NewName", number = Inf)

results_anova_p<-as.data.frame(p.adjust(results$P.Value, method = "BH"))
significant_results <- results[results$logFC < 0 & results$P.Value < 0.05, ]


