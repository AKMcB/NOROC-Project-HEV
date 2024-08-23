#############
# Libraries #
#############

library(data.table)
library(fgsea)
library(ggpubr)
library(tidyverse)

#####################
# Read GSEA results #
#####################

high_all<-as.data.frame(fread("hallmark_cluster1_vs_cluster2_all_proteins_noroc.Gsea.1707483329832/gsea_report_for_cluster1_1707483329832.tsv"))
high_all <- high_all[1:20,]

high_1 <- as.data.frame(fread("hallmark_cluster1_vs_cluster2_all_proteins_noroc.Gsea.1707483329832/gsea_report_for_cluster2_1707483329832.tsv"))
high_1 <- high_1[1:20,]

high_2 <- as.data.frame(fread("kegg_legacy_cluster2_high_vs_low_hev_all_proteins_noroc.Gsea.1707488233830/gsea_report_for_High_1707488233830.tsv"))
high_2 <- high_2[1:20,]

high_all <- high_all[,c(1,6,8)]
high_1 <- high_1[,c(1,6,8)]
high_2 <- high_2[,c(1,6,8)]

high_all <- high_all %>% mutate(source = "Cluster 1")
high_1 <- high_1 %>% mutate(source = "Cluster 2")
high_2 <- high_2 %>% mutate(source = "HEV High Cluster 2")

merged <- as.data.frame(bind_rows(high_all, high_1))

merged <- merged[complete.cases(merged$NAME), ]

merged$NAME <- gsub("HALLMARK_", "", merged$NAME)
merged$NAME <- gsub("_", " ", merged$NAME)

# Order the combinations of results
pathway_counts <- merged %>%
  group_by(NAME) %>%
  summarise(count = n_distinct(source),
            source_combination = paste(sort(unique(source)), collapse = ", ")) %>%
  mutate(priority = case_when(
    source_combination == "HEV High All Patients, HEV High Cluster 1, HEV High Cluster 2" ~ 1,
    source_combination == "HEV High All Patients, HEV High Cluster 1" ~ 2,
    source_combination == "HEV High All Patients, HEV High Cluster 2" ~ 3,
    source_combination == "HEV High Cluster 1,HEV High Cluster 2" ~ 4,
    source_combination == "HEV High All Patients" ~ 5,
    source_combination == "HEV High Cluster 1" ~ 6,
    source_combination == "HEV High Cluster 2" ~ 7,
    count == 1 ~ 7,  # Unique pathways last
    TRUE ~ 8  # Catch all for any other combinations
  )) %>%
  arrange(priority,(count))

merged_1 <- merge(merged, pathway_counts, by = "NAME" )

average_scores <- merged_1 %>%
  group_by(NAME) %>%
  summarize(avg_normalized_score = mean(NES))

ordered_pathways <- average_scores %>%
  arrange((avg_normalized_score)) %>%
  pull(NAME)

merged_1$NAME <- factor(merged_1$NAME, levels = ordered_pathways)

#############
# Tileplots #
#############

nes <- ggplot(merged_1, aes(x = source, y = reorder(NAME, priority, decreasing = T) , fill = NES)) +
  geom_tile(color = "black", size = 0.6)  +
  scale_fill_gradient(low = "lightgrey", high = "darkred") +
  theme_transparent() +
  labs(title = "Enrichment Score of Hallmarks Between Patient Clusters", 
       x = "", 
       y = "Pathways") + 
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),  
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -1))+
  theme(aspect.ratio = 3) 

nes

fdr <- ggplot(merged, aes(x = source, y = NAME, fill = `FDR q-val`)) +
  geom_tile(color = "black", size = 0.6)  +
  scale_fill_gradient(low = "lightgrey", high = "darkred") +
  theme_minimal() +
  labs(title = "Significance of Shared Hallmarks", 
       x = "", 
       y = "Pathways") + 
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),  
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -1))+
  theme(aspect.ratio = 3)


dev.off()

pdf("shared_hallmarks_patient_clusters.pdf", width = 8, height = 8, onefile = F)
print(nes)
dev.off()

pdf("shared_hallmarks_hev_high_fdr_v2.pdf", width = 8, height = 8, onefile = F)
print(fdr)
dev.off()
