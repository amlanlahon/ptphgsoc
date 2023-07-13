library(scales)0000000000

vol_data_logFC1 = read.csv("volcano_data_logFC1.csv")

with(vol_data_logFC1, plot(-1*logFC, -log10(adj.P.Val), pch=20, main="Volcano plot", col = alpha("black", 0.2), xlim=c(-6.2,7), xlab = "logFC", ylab = "-log10(pVal)"))

with(subset(vol_data_logFC1, adj.P.Val<0.01 & abs(logFC)>1), points(-1*logFC, -log10(adj.P.Val), pch=1, col=alpha("orange", 0.8)))

with(res_table, plot(-1*logFC, -log10(adj.P.Val), pch=20, main="Volcano plot", col = alpha("black", 0.2), xlim=c(-6.2,7), xlab = "logFC", ylab = "-log10(pVal)"))

with(subset(res_table, adj.P.Val<0.01 & abs(logFC)>1), points(-1*logFC, -log10(adj.P.Val), pch=1, col=alpha("orange", 0.8)))
