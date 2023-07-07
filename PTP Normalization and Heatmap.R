library(preprocessCore)
library(circlize)
library(ComplexHeatmap)
library(viridis)

ptp = read.csv("PTP Genes_DEGs_LogFC1 - Copy.csv")
ptp_norm = normalize.quantiles(as.matrix(ptp))

max(ptp_norm)
min(ptp_norm)
ptp_norm.scaled = scale(ptp_norm, center = T, scale = T)
max(ptp_norm.scaled)
min(ptp_norm.scaled)
write.csv(ptp_norm.scaled, "ptp_normalized.csv")

ptp_norm_file = read.csv("ptp_normalized.csv", row.names = 1)
ptp_mat_norm = as.matrix(ptp_norm_file)
f4 = colorRamp2(seq(min(ptp_norm.scaled), max(ptp_norm.scaled),length = 3), c("blue", "white", "red"), space = "LAB")
Heatmap(ptp_mat_norm, 
        col = f4, 
        cluster_rows = T, 
        cluster_columns = F, 
        show_row_names = T,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),
        )
