library(limma) 
library(dplyr)

library(edgeR)
cluster_results<-read.csv("cluster_results.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
cluster_results_1 <- cluster_results %>% 
  filter(Cluster_Labels == 0)

response<-read.csv("Womac_score_pain_function.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
response <- response %>% filter(ID %in% cluster_results_1$ID)

##just for removing A from womac func

response <- response[complete.cases(response$Womac.Function), ]

count_data <- read.csv("metabolites.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
data <- count_data %>% filter(ID %in% response$ID)
data<-data[2:length(data)]



###string and mild response are 1 and 2, no response is 3. Therefre, we are binarizing
response$Womac.Function_up<-ifelse(response$Womac.Function==3,0,1)

design <- model.matrix(~ factor(response$Womac.Function_up))
colnames(design)[colnames(design) == "factor(response$Womac.Function_up)1"] <- "NewName"

y <- DGEList(counts = t(data), group = response$Womac.Function_up)

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
significant_results <- results[results$P.Value < 0.05, ]

