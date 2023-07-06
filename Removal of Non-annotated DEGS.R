annot_genes = read.csv("annotated_genes.csv", row.names = 1)
degs = read.csv("expression_values_of_degs_logFC1.csv", row.names = 1)
write.csv(degs[match(rownames(annot_genes), rownames(degs)),], "annotated_DEGs_with_NA.csv")
Degs_NA = read.csv("annotated_DEGs_with_NA.csv", row.names = 1)
annot_degs = na.omit(Degs_NA)
write.csv(annot_degs, "annotated_degs.csv")
annot_degs = read.csv("annotated_degs.csv")
annot_degs$X.1 = NULL
write.csv(annot_degs, "annotated_degs.csv")
