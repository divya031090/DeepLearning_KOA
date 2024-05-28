#Comparing standardized mean differences/diff analysis/mean values

library(tableone)
cluster_results<-read.csv("cluster_results.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")


##comparing cluster 1 vs 2 + 3 and so on (clustering is in labels 0,1,2)
library(edgeR)
data<-data1 ### RNA or metabolite data
cluster_results$Cluster_Labels_up<-ifelse(cluster_results$Cluster_Labels==2,0,1)

design <- model.matrix(~ factor(cluster_results$Cluster_Labels_up))
colnames(design)[colnames(design) == "factor(cluster_results$Cluster_Labels_up)1"] <- "NewName"

y <- DGEList(counts = t(data), group = cluster_results$Cluster_Labels_up)

# Calculate normalizastion factors
y <- calcNormFactors(y)

v <- voom(y, design, plot = TRUE)

data<-as.data.frame(t(v$E))
data$cluster <- cluster_results$Cluster_Labels_up
vars=colnames(data)

tabUnmatched <- CreateTableOne(vars = vars, strata = "cluster", data = data,smd=TRUE)

result_smd<-as.data.frame(print(tabUnmatched,smd=TRUE,exactLocation = FALSE))
colnames(result_smd)<-c("Cluster0","Cluster1","p","test","SMD")
new_rownames <- gsub("\\(.*?\\))", "", rownames(result_smd))

# Assign the modified row names back to the data frame
rownames(result_smd) <- new_rownames

result_smd$Cluster0 <- gsub("\\(.*?\\)", "", result_smd$Cluster0)
result_smd$Cluster1 <- gsub("\\(.*?\\)", "", result_smd$Cluster1)
result_smd<-result_smd[2:length(result_smd[,1]),1:5]
result_smd_fil <- result_smd[result_smd$Cluster0 > result_smd$Cluster1, ]
result_smd_fil$p<-as.numeric(as.character(result_smd_fil$p))
result_smd_fil<-result_smd_fil[result_smd_fil$p<0.05,]

