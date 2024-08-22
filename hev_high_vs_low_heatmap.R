
###############
## Libraries ##
###############
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)
library(circlize)
library(dendextend)
library(tidyverse)
library(clValid)
library(dendsort)
library(fgsea)

##########################
## Read expression file ##
##########################

expr <- read.csv2("example_expr_file", sep = ";", as.is = T,check.names = F)
expr <- expr[,-c(2,3)]
expr <- distinct(expr,expr$`Gene names`,.keep_all = T) #to remove duplicates 
expr$`Gene names` <- gsub("\\;.*","",expr$`Gene names`)
#be aware that the distinct function creates a column at the end of the df 
#remove it before continuing
expr$`expr$\`Gene names\`` <- NULL 

gens <- as.data.frame(gmtPathways("gene_lists/HALLMARK_ANGIOGENESIS.v2023.2.Hs (2).gmt"))

expr <- subset(expr, expr$`Gene names` %in% gens$HALLMARK_ANGIOGENESIS)

rownames(expr) <- expr[,1]
expr$`Gene names` <- NULL
expr <- as.data.frame(t(expr))
expr <- rownames_to_column(expr, "id")

########################
## Read clinical info ##
########################

info <- read.csv2("example_clinical_file.csv", as.is = T, check.names = F)
info <- info[,c(1,53)]
info<- info[complete.cases(info$HEV_score_LowHigh_mean),]

colnames(info)[1] <- "PID"

merged <- merge(expr, info, by.x = "id", by.y = "PID")

#create annotation
ann <- list(merged$id, merged$HEV_score_LowHigh_mean)
ann <- as.data.frame(ann)
colnames(ann) <- c("id", "HEV_level")
rownames(ann) <- ann[,1]
ann$id <- NULL

#2=high 1=low
ann$HEV_level <- sub("2", "High", ann$HEV_level)
ann$HEV_level <- sub("1", "Low", ann$HEV_level)

merged$HEV_score_LowHigh_mean <- NULL
rownames(merged) <- merged$id
merged$id <- NULL
merged <- as.data.frame(t(merged))
str(merged)
merged <- merged %>% mutate_all(as.numeric)


#################################
## Calculate cluster stability ##
#################################

#Patients has to be in rownames therefore use t(merged) in the function 
#This will not change the df and the genes will still be in rownames when 
#converting the df into a matrix
internal <- clValid(t(merged), method = "complete", metric = "correlation", clMethods = "hierarchical", nClust = 2:10, validation = "internal")
plot(internal, legend = FALSE)


####################
## Create Heatmap ##
####################

#The genes should be in the rownames
h <- as.matrix(merged) #Create a matrix of the main expression file 
str(h)

#Check skewness of the data
#Should be a sharp peak if the data is normalized
dat <- data.frame(values = as.numeric(h))

ggplot(dat, aes(values)) + geom_density(bw = "sj")


##We will plot based on quantiles of the expression values
#First find the breaks- 10 breaks from the lowest to the highest expression 
h_breaks <- seq(min(h), max(h), length.out = 10)
h_breaks



#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

h_breaks <- quantile_breaks(h, n = 11)
h_breaks

#Make clusters based on their pearson correlation
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")

#We can use the dendsort package and reorder the clustering
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

#Choose the color gradient here. You can change the colors by googling R color codes.
#Run these code to add the option 
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))


a = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))


#Get the annotation files. Make sure this file is in the same order as the expression matrix file.
#For example, order of genes/patients in expr == order of genes/patients in annotation file.
#Otherwise the annotation will show misleading information.


ha <- HeatmapAnnotation ("HEV Score"= ann$HEV_level,
                         col = list("HEV Score"= c("High"= "#ff0000",
                                                   "Low" = "#1e90ff")),
                         simple_anno_size = unit(0.30, "cm"),
                         border = T,
                         annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                         show_legend = T,
                         annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))


ht <- Heatmap(h,col = col_fun,
              cluster_columns = Colv,
              name = "Expression Values",
              show_heatmap_legend = T ,
              top_annotation = ha,
              #left_annotation = hr,
              show_column_names = F,
              show_row_names = F,cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=10),
              column_title_gp = gpar(fontsize=10),
              height = unit(12, "cm"),
              width = unit(10.55, "cm"),
              column_dend_reorder = T,
              row_dend_reorder = T,
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 90,
              #legend_height = unit(4, "cm"),
              row_names_gp = grid::gpar(fontsize= 7.5),
              column_names_gp = grid::gpar(fontsize = 10,fontface="bold"),
              heatmap_legend_param = list(title="Scaled Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)),
              row_split = 3, column_split = 2,
              column_gap = unit(c(0.7,0.7), "mm"),
              row_gap = unit(c(0.7,0.7,0.7), "mm"))

#Print the heatmap as a pdf in your local drive. Use padding to make it look better (in my opinion)
draw(ht, merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))

#You have to run dev.off() before checking the pdf in your local drive

dev.off()

#To find the gene and patient orders, draw the heatmap object first so that it will not change with every run
ht <- draw(ht)

#For getting the column orders
for (i in 1:length(column_order(ht)))   if (i == 1) {
  clu <- t(t(colnames(merged[,column_order(ht)[[i]]])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("Patient_ID", "Cluster")   } else {
    clu <- t(t(colnames(merged[,column_order(ht)[[i]]])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)   } 

out

write.csv2(out, "cluster_column_order_high_low_hev_updated_patients.csv")


#For getting the gene orders
for (i in 1:length(row_order(ht))){   if (i == 1) {
  clu <- t(t(row.names(merged[row_order(ht)[[i]],])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("GeneID", "Cluster")   } else {
    clu <- t(t(row.names(merged[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)   } 
}

out

write.csv2(out, file = "cluster_rows_order_high_low_hev_updated_patients.csv")

pdf("heatmap_high_low_hev_updated_patients.pdf", height = 8, width = 8)
print(ht)
dev.off()

png("heatmap_high_low_hev_updated_patients.png", res=200, height = 2000, width = 2000)
print(ht)
dev.off() 
