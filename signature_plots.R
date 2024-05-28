library(dplyr)
library(tidyr)
library(ggplot2)
cluster_results<-read.csv("cluster_results.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")

data1<-read.csv("metabolites.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
data1 <- data1 %>% filter(ID %in% cluster_results$ID)
data1<-data1[2:length(data1)]
colnames(data1) <- paste(colnames(data1), "pla", sep = "_")
data2<-read.csv("RNA_plasma.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
data2 <- data2 %>% filter(ID %in% cluster_results$ID)
data2<-data2[2:length(data2)]
colnames(data2) <- paste(colnames(data2), "pla", sep = "_")
data3<-read.csv("RNA_synovial.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
data3 <- data3 %>% filter(ID %in% cluster_results$ID)
data3<-data3[2:length(data3)]
colnames(data3) <- paste(colnames(data3), "syn", sep = "_")
data4<-read.csv("RNA_urine.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")
data4 <- data4 %>% filter(ID %in% cluster_results$ID)
data4<-data4[2:length(data4)]
colnames(data4) <- paste(colnames(data4), "uri", sep = "_")

variables<-read.csv("Variables_smd_DE.csv",header = TRUE, stringsAsFactors = FALSE,sep=",")

full_data<-cbind(data1,data2,data3,data4)


###For cluster 1


variables$Cluster.1_normalized <- tolower(trimws(variables$Cluster.1))
full_data_names_normalized <- tolower(trimws(names(full_data)))



matching_names <- full_data_names_normalized %in% variables$Cluster.1_normalized

# Subset full_data using matching variable names
filtered_data <- full_data[, matching_names]


posterior.states<-cluster_results
posterior.states$state<-factor(cluster_results$Cluster_Labels)

mydata_new<-filtered_data

plot.data <- cbind(mydata_new, posterior.states) %>% 
  gather(key="measure", value="value", colnames(mydata_new[1]):colnames(mydata_new[length(mydata_new)]))

summary.plot.data <- plot.data %>% 
  group_by(state, measure) %>% 
  summarize(z=mean(value))

summary.plot.data$state<-factor(summary.plot.data$state)


#### order based on decreasing values

state_0 <- summary.plot.data[summary.plot.data$state == 0, ]
state_1 <- summary.plot.data[summary.plot.data$state == 1, ]
state_2 <- summary.plot.data[summary.plot.data$state == 2, ]
fulldata<-cbind(state_0,state_1,state_2)
full_data_sorted <- fulldata[order(-fulldata$z), ]
summary.plot.data_sorted<-rbind(full_data_sorted[,1:4],full_data_sorted[,5:8],full_data_sorted[,9:12])



summary.plot.data <- summary.plot.data %>%
  mutate(state = case_when(
    state == 0 ~ "Cluster1",
    state == 1 ~ "Cluster2",
    state == 2 ~ "Cluster3",
    TRUE ~ as.character(state)  # Keep other values unchanged
  ))

summary.plot.data$measure <- factor(summary.plot.data$measure, levels = unique(summary.plot.data$measure))

###
p <- ggplot(summary.plot.data, aes(y = z, x = measure, group = state, color = state)) + 
  geom_point() + 
  geom_line() +
  scale_y_continuous(breaks = 0) +
  
  theme_minimal()+
  labs( x = "Feature", y = "Mean Value", color = "Cluster")+
  theme(axis.text.x = element_text(size=12,face="bold",angle = 90, vjust = 0.5, hjust = 1))

p

###For cluster 2


variables$Cluster.2_normalized <- tolower(trimws(variables$Cluster.2))
full_data_names_normalized <- tolower(trimws(names(full_data)))



matching_names <- full_data_names_normalized %in% variables$Cluster.2_normalized

# Subset full_data using matching variable names
filtered_data <- full_data[, matching_names]


posterior.states<-cluster_results
posterior.states$state<-factor(cluster_results$Cluster_Labels)

mydata_new<-filtered_data

plot.data <- cbind(mydata_new, posterior.states) %>% 
  gather(key="measure", value="value", colnames(mydata_new[1]):colnames(mydata_new[length(mydata_new)]))

summary.plot.data <- plot.data %>% 
  group_by(state, measure) %>% 
  summarize(z=mean(value))

summary.plot.data$state<-factor(summary.plot.data$state)


#### order based on decreasing values

state_0 <- summary.plot.data[summary.plot.data$state == 0, ]
state_1 <- summary.plot.data[summary.plot.data$state == 1, ]
state_2 <- summary.plot.data[summary.plot.data$state == 2, ]
fulldata<-cbind(state_0,state_1,state_2)
colnames(fulldata)[8]<-"Cluster2z"
full_data_sorted <- fulldata[order(-fulldata$Cluster2z), ]
colnames(full_data_sorted)[8]<-"z"
summary.plot.data_sorted<-rbind(full_data_sorted[,1:4],full_data_sorted[,5:8],full_data_sorted[,9:12])



summary.plot.data <- summary.plot.data_sorted %>%
  mutate(state = case_when(
    state == 0 ~ "Cluster1",
    state == 1 ~ "Cluster2",
    state == 2 ~ "Cluster3",
    TRUE ~ as.character(state)  # Keep other values unchanged
  ))

summary.plot.data$measure <- factor(summary.plot.data$measure, levels = unique(summary.plot.data$measure))

###

p <- ggplot(summary.plot.data, aes(y = z, x = measure, group = state, color = state)) + 
  geom_point() + 
  geom_line() +
  scale_y_continuous(breaks = 0) +
  
  theme_minimal()+
  labs( x = "Feature", y = "Mean Value", color = "Cluster")+
  theme(axis.text.x = element_text(size=12,face="bold",angle = 90, vjust = 0.5, hjust = 1))

p

###For cluster 3


variables$Cluster.3_normalized <- tolower(trimws(variables$Cluster.3))
full_data_names_normalized <- tolower(trimws(names(full_data)))



matching_names <- full_data_names_normalized %in% variables$Cluster.3_normalized

# Subset full_data using matching variable names
filtered_data <- full_data[, matching_names]


posterior.states<-cluster_results
posterior.states$state<-factor(cluster_results$Cluster_Labels)

mydata_new<-filtered_data

plot.data <- cbind(mydata_new, posterior.states) %>% 
  gather(key="measure", value="value", colnames(mydata_new[1]):colnames(mydata_new[length(mydata_new)]))

summary.plot.data <- plot.data %>% 
  group_by(state, measure) %>% 
  summarize(z=mean(value))

summary.plot.data$state<-factor(summary.plot.data$state)
#### order based on decreasing values

state_0 <- summary.plot.data[summary.plot.data$state == 0, ]
state_1 <- summary.plot.data[summary.plot.data$state == 1, ]
state_2 <- summary.plot.data[summary.plot.data$state == 2, ]
fulldata<-cbind(state_0,state_1,state_2)
colnames(fulldata)[12]<-"Cluster3z"
full_data_sorted <- fulldata[order(-fulldata$Cluster3z), ]
colnames(full_data_sorted)[12]<-"z"
summary.plot.data_sorted<-rbind(full_data_sorted[,1:4],full_data_sorted[,5:8],full_data_sorted[,9:12])



summary.plot.data <- summary.plot.data_sorted %>%
  mutate(state = case_when(
    state == 0 ~ "Cluster1",
    state == 1 ~ "Cluster2",
    state == 2 ~ "Cluster3",
    TRUE ~ as.character(state)  # Keep other values unchanged
  ))

summary.plot.data$measure <- factor(summary.plot.data$measure, levels = unique(summary.plot.data$measure))

###
p <- ggplot(summary.plot.data, aes(y = z, x = measure, group = state, color = state)) + 
  geom_point() + 
  geom_line() +
  scale_y_continuous(breaks = 0) +
  
  theme_minimal()+
  labs( x = "Feature", y = "Mean Value", color = "Cluster")+
  theme(axis.text.x = element_text(size=12,face="bold",angle = 90, vjust = 0.5, hjust = 1))

p




