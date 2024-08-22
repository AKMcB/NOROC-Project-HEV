
library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(grid)
library(gridExtra)
library(data.table)
#-------------------------------------------------------------------------------------
#Subsetting the TCGA dataset based on Cluster

expr <- read.csv2("example_expr_file.csv", sep = ";", as.is = T, check.names = F)
expr <- expr[,-c(2,3)]
expr$`Gene names` <- gsub("\\;.*","",expr$`Gene names`)

expr <- distinct(expr,expr$`Gene names`,.keep_all = T)
expr$`expr$\`Gene names\`` <- NULL
colnames(expr)[1] <- "id"
rownames(expr) <- expr[,1]
expr$id <- NULL
expr <- as.data.frame(t(expr))
expr <- tibble::rownames_to_column(expr, "id")


patient <- read.csv2("example_cluster_info_from_heatmap.csv")
patient <- patient[,-1]

merged <- merge(expr, patient, by.x="id", by.y = "Patient_ID")

info <- read.csv2("example_clinical_file.csv")
info <- info[,c(1,52:53)]

info<- info[complete.cases(info$HEV_score_Ibrahim),]
info$HEV_score_LowHigh_mean <- sub("2", "High", info$HEV_score_LowHigh_mean)
info$HEV_score_LowHigh_mean <- sub("1", "Low",info$HEV_score_LowHigh_mean)

merged <- merge(merged, info, by.x="id", by.y = "PID")

nonMatch_Uniquedf2 <- info %>% 
  filter(!info$PID %in% expr$id)


#Divide data based on either high or low HEV level 
#2=High 1=Low
high <- subset(merged, merged$HEV_score_LowHigh_mean == "2")
low <- subset(merged, merged$HEV_score_LowHigh_mean == "1")

cluster1 <- subset(merged, merged$Cluster == "cluster1")
cluster2 <- subset(merged, merged$Cluster == "cluster2")


merged_1 <- subset(merged, merged$Cluster == "cluster1")
high_1 <- subset(merged_1, merged_1$HEV_score_LowHigh_mean == "High")
low_1 <- subset(merged_1, merged_1$HEV_score_LowHigh_mean == "Low")

merged_2 <- subset(merged, merged$Cluster == "cluster2")
high_2 <- subset(merged_2, merged_2$HEV_score_LowHigh_mean == "High")
low_2 <- subset(merged_2, merged_2$HEV_score_LowHigh_mean == "Low")

#The patient id must include either High or low 
#High/low included so to compare these two groups
#This places High/Low before patient id

rownames(high) <- high[,1]
high <- high[,-1]
high <- as.data.frame(t(high))
colnames(high) <- paste("High",colnames(high), sep= "_")

rownames(low) <- low[,1]
low <- low[,-1]
low <- as.data.frame(t(low))
colnames(low) <- paste("Low",colnames(low), sep= "_")

rownames(cluster1) <- cluster1[,1]
cluster1 <- cluster1[,-1]
cluster1 <- as.data.frame(t(cluster1))
colnames(cluster1) <- paste("cluster1",colnames(cluster1), sep= "_")

rownames(cluster2) <- cluster2[,1]
cluster2 <- cluster2[,-1]
cluster2 <- as.data.frame(t(cluster2))
colnames(cluster2) <- paste("cluster2",colnames(cluster2), sep= "_")

rownames(high_1) <- high_1[,1]
high_1 <- high_1[,-1]
high_1 <- as.data.frame(t(high_1))
colnames(high_1) <- paste("high_1",colnames(high_1), sep= "_")

rownames(low_1) <- low_1[,1]
low_1 <- low_1[,-1]
low_1 <- as.data.frame(t(low_1))
colnames(low_1) <- paste("low_1",colnames(low_1), sep= "_")

rownames(high_2) <- high_2[,1]
high_2 <- high_2[,-1]
high_2 <- as.data.frame(t(high_2))
colnames(high_2) <- paste("high_2",colnames(high_2), sep= "_")

rownames(low_2) <- low_2[,1]
low_2 <- low_2[,-1]
low_2 <- as.data.frame(t(low_2))
colnames(low_2) <- paste("low_2",colnames(low_2), sep= "_")

#----------------------------------------------------------------------------------------------
#Merging the expression files for GSEA
#Merge the files 
comb <- cbind(cluster1, cluster2)

comb_1 <- cbind(high_1, low_1)

comb_2 <- cbind(high_2, low_2)

#Create a phenotype file-include only high/Low  
pheno <- comb[4651, ]
comb <- comb[-c(4651:4653),]

pheno_1 <- comb_1[4653, ]
comb_1 <- comb_1[-c(4651:4653),]

pheno_2 <- comb_2[4653, ]
comb_2 <- comb_2[-c(4651:4653),]

#Remove everything after ; in gene names in rownames
#Just keeping the first gene
#comb <- tibble::rownames_to_column(comb, "gene")
#comb$gene <- gsub("\\;.*","",comb$gene)

#check for duplications
dup <- comb_2$gene[duplicated(comb_2$gene)]

#save file 
#Edit the files in excel to fit the requirements for GSEA
write.csv2(comb_2, "GSEA/NOROC_cluster2_hev_high_low_all_proteins_GSEA.csv")

write.csv2(pheno_2, "GSEA/NOROC_cluster2_hev_high_low_all_proteins_pheno_GSEA.csv") 
