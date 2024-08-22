#############
# Libraries #
#############
library(tidyverse)
library(data.table)

#Read cluster information
id <- as.data.frame(fread("example_cluster_info_from_heatmap.csv", header = TRUE))
id$V1 <- NULL

info <- read.csv2("example_clinical_file.csv")
info <- info[,c(1,52:53)]

merged <- merge(id, info, by.x = "Patient_ID", by.y="PID")
merged$HEV_score_LowHigh_mean <- sub("2", "High", merged$HEV_score_LowHigh_mean )
merged$HEV_score_LowHigh_mean  <- sub("1", "Low", merged$HEV_score_LowHigh_mean )

merged$Cluster <- gsub("cluster1", "Cluster 1", merged$Cluster)
merged$Cluster <- gsub("cluster2", "Cluster 2", merged$Cluster)

counts_merged <- table(merged$Cluster,merged$HEV_score_LowHigh_mean)
chi_square <- chisq.test(counts_merged)
chi_square 

chi_square_label <- paste("Chi-Square =", round(chi_square$statistic, 2), 
                          "\nDegrees of Freedom =", chi_square$parameter,
                          "\nP-Value =", format.pval(chi_square$p.value, digits = 3)) 

counts_merged <- as.data.frame(counts_merged)

counts_merged <- counts_merged %>%
  group_by(Var1) %>% 
  mutate(proportion = (Freq/sum(Freq)*100))

counts_merged$proportion <- round(counts_merged$proportion, 2)

p <- ggplot(counts_merged, aes(x = Var1, y = proportion, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(proportion, "%"), y = proportion), 
            position = position_stack(vjust = 0.5), size = 3, color = "black") +
  labs(x = "", y = "Percentage", title = "HEV Level in Heatmap Clusters") +
  theme_classic() +
  guides(fill = guide_legend(title = "HEV Level"))+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5), # Italic if it is a gene. 
        axis.text.x = element_text(size = 10), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))+ 
  scale_fill_manual(values = c(
    "High" = "#ff0000",
    "Low" = "#1e90ff")) + 
  annotate("text", x = Inf, y = Inf, label = chi_square_label,
           hjust = 1, vjust = 1, size = 3)
p

pdf("hevlevel_heatmap_cluster_noroc_updated_normalized.pdf", height = 11, width = 10)
print(p)
dev.off()

png("hevlevel_heatmap_cluster_noroc_updated_normalized.png", res = 200, height = 2200, width = 1800)
print(p)
dev.off()
