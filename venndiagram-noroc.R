#############
# Libraries #
#############

library(data.table)
library(VennDetail)
library(ggvenn)

#####################
# Read gsea results #
#####################

high_all <- as.data.frame(fread("GSEA/hallmark_leg_hev_high_vs_low_all_proteins.Gsea.1705919994552/gsea_report_for_High_1705919994552.tsv"))
high_all <- high_all[1:20,] #only if there are 20 more pathways


high_clus1 <- as.data.frame(fread("GSEA/hallmark_legacy_cluster1_high_vs_low_hev_all_proteins_noroc.Gsea.1707487969762/gsea_report_for_High_1707487969762.tsv"))
high_clus1 <- high_clus1[1:20,]


high_clus2 <- as.data.frame(fread("GSEA/hallmark_legacy_cluster2_high_vs_low_hev_all_proteins_noroc.Gsea.1707488233830/gsea_report_for_High_1707488233830.tsv"))
high_clus2 <- high_clus2[1:20,]

all <- high_all
cluster1 <- high_clus1
cluster2 <- high_clus2


all <- all[,"NAME"]
cluster1 <- cluster1[,"NAME"]
cluster2 <- cluster2[,"NAME"]

################
# Venn Diagram #
################

overlap_list <- list("All Patients" = all, "Cluster 1" = cluster1, "Cluster 2" = cluster2)
overlap_plot <- venndetail(overlap_list) #identify overlapping terms
overlap_df <- result(overlap_plot)

fwrite(overlap_df, "overlap_hev_high_all_cluster1_cluster2_gsea_hallmark_only_high.csv", sep = ",")


x = list(all,cluster1, cluster2)

names(x) <- c("All Patients", "Cluster1", "Cluster2")

plot <- ggvenn(x,set_name_size = 3,
      fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),stroke_size = 0.5)

plot

dev.off()

ggexport(plot, filename = "GSEA_NOROC_High_hev_venndiagram_hallmark.png",res=200,width = 2500, height = 2000)
ggexport(plot, filename = "GSEA_NOROC_High_hev_venndiagram_hallmark.pdf",width = 15, height = 15)
