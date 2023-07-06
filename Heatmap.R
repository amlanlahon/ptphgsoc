`#######################
#HEATMAP

heatmap_data = read.csv("annotated_DEGs.csv", row.names = 2)
heatmap_data$X.1 = NULL
heatmap_data_mat = as.matrix(heatmap_data)
f1 = colorRamp2(seq(min(heatmap_data_mat), max(heatmap_data_mat),length = 3), c("blue", "white", "red"), space = "LAB")
Heatmap(heatmap_data_mat, 
        col = f2, 
        cluster_rows = T, 
        cluster_columns = T, 
        show_row_names = F,
        column_names_gp = gpar(fontsize = 5))

######################
#INTERACTIVE HEATMAP

f2 = colorRamp2(seq(min(exp_data.norm.scaled_logFC1), max(exp_data.norm.scaled_logFC1), length = 3), 
                c("blue", "white", "red"),
                space = "RGB")
Heatmap(exp_data.norm.scaled_logFC1,
        col = f2, 
        cluster_rows = T, 
        cluster_columns = F, 
        show_row_names = F,
        column_names_gp = gpar(fontsize = 5),
        show_column_names = T)
library(heatmaply)
dir.create("interactive_hm_new")
heatmaply(exp_data.norm.scaled_logFC1, file = "interactive_hm_new/hm_int_logFC1.html")

#########################`
