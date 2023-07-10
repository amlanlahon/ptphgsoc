#n=10868
net_dat_degs = read.csv("C:/Users/Chocobunbhai/Documents/Network_data_woduplicates.csv", row.names = 1, stringsAsFactors = F)
net_dat_degs = t(scale(t(net_dat_degs), center = T, scale = T))

datExpr = t(net_dat_degs)
library(WGCNA)

# Choose a set of soft-thresholding powers

powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 1.5;

# Scale-free topology fit index as a function of the soft-thresholding power

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h

abline(h=0.88,col="red")

# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

POW = 9
write.csv(net_dat_degs, file = "data_network_pids_10_12_2019.csv")

# here we define the adjacency matrix using soft thresholding with beta=6

ADJ1=abs(cor(datExpr,use="p"))^POW

# When you have relatively few genes (<5000) use the following code

#k=as.vector(apply(ADJ1,2,sum, na.rm=T))

# When you have a lot of genes use the following code

k=softConnectivity(datE=datExpr,power=6)

# Plot a histogram of k and a scale free topology plot

sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

# Turn adjacency into a measure of dissimilarity

dissADJ=1-ADJ1
dissTOM=TOMdist(ADJ1)
collectGarbage()

hierADJ=hclust(as.dist(dissADJ), method="average" )

# Plot the resulting clustering tree together with the true color assignment

sizeGrWindow(10,5);

#plotDendroAndColors(hierADJ, dendroLabels = FALSE, hang = 0.03,

#main = "Gene hierarchical clustering dendrogram and simulated module colors" )

plot(hierADJ, xlab="", sub="", main = "Gene clustering on hierarchical clustering",
     labels = FALSE, hang = 0.001)

branch.number=cutreeDynamic(hierADJ,method="tree")

# This function transforms the branch numbers into colors

colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))

colorDynamicADJ=labels2colors(branch.number )

colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ,
                                                  cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))

# Plot results of all module detection methods together:

sizeGrWindow(10,5)
plotDendroAndColors(dendro = hierADJ,
                    colors=data.frame(colorStaticADJ,
                                      colorDynamicADJ, colorDynamicHybridADJ),
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Gene dendrogram and module colors")

### TOM ####

# Calculate the dendrogram

hierTOM = hclust(as.dist(dissTOM),method="average");

# The reader should vary the height cut-off parameter h1

# (related to the y-axis of dendrogram) in the following

#colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))

colorDynamicTOM = labels2colors(cutreeDynamic(hierTOM,method="tree", minClusterSize = 20, deepSplit = 4))

#colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,

#deepSplit=4, pamRespectsDendro = FALSE))

# Now we plot the results

# sizeGrWindow(10,5)

# plotDendroAndColors(hierTOM,

#                     colors=data.frame(colorStaticTOM,

#                                       colorDynamicTOM, colorDynamicHybridTOM),

#                     dendroLabels = FALSE, marAll = c(1, 8, 3, 1),

#                     main = "Gene dendrogram and module colors, TOM dissimilarity")

sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")


dynamicMods = cutreeDynamic(hierTOM,method="tree", minClusterSize = 20, deepSplit = 4)
table(dynamicMods)
mods<- table(dynamicMods)
write.csv(mods, "modules.csv")
dynamicColors = labels2colors(dynamicMods)
moduleColors = dynamicColors
table(dynamicColors)
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
ME_unmerged = MEList
ME_unmerged = MEList$eigengenes

##go to line 177

datTraits = trait_table
moduleTraitCor = cor(ME_unmerged, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3, 3, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(ME_unmerged),
               ySymbols = names(ME_unmerged),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
"mod_trait_unmerged_big_font"
MEs = MEList$eigengenes

##############################################

MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
MEDissThres = 0.25
par(mfrow=c(1,1))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
#moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(mergedColors, colorOrder)-1;
#MEs = mergedMEs;


sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,
                    colors=(cbind(moduleColors,mergedColors)),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and merged module colors, TOM dissimilarity")


nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
moduleColors = mergedColors
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

trait_table = data.frame(rownames(datExpr))
trait_table$normal = c(rep(1,12), rep(0,57))
trait_table$cancer = c(rep(0,12), rep(1,57))
rownames(trait_table) = trait_table$rownames.datExpr.
trait_table = trait_table[,-1]
datTraits = trait_table
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 9, 3, 3));

# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
cancerstatus = as.data.frame(datTraits$cancer)
names(cancerstatus) = "Cancerstatus"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, cancerstatus, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(cancerstatus), sep="");
names(GSPvalue) = paste("p.GS.", names(cancerstatus), sep="");
genes =data.frame(rownames(net_dat_degs), stringsAsFactors = F)

module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module
red_genes = genes[moduleGenes,]
column = match("brown", modNames);
moduleGenes = moduleColors=="brown"
brown_genes = genes[moduleGenes,]
column = match("darkolivegreen", modNames);
moduleGenes = moduleColors=="darkolivegreen"
darkolivegreen_genes = genes[moduleGenes,]
column = match("lightgreen", modNames);
moduleGenes = moduleColors=="lightgreen"
lightgreen_genes = genes[moduleGenes,]
column = match("grey", modNames);
moduleGenes = moduleColors=="grey"
grey_genes = genes[moduleGenes,]


write.csv(c(red_genes, brown_genes, darkolivegreen_genes, lightgreen_genes, grey_genes), "network_genes.csv")
write.csv(GSPvalue, "gene_Sig.csv")
write.csv(geneModuleMembership,"modmem.csv")

###################################3

# module membership vs gene significance p value

library(ggplot2)
scatter_data = data.frame(cbind(abs(geneModuleMembership$MMbrown),abs(geneTraitSignificance$GS.Cancerstatus)))

#scatter_data = data.frame(cbind((geneModuleMembership$MMgrey),(geneTraitSignificance$GS.Cancerstatus)))

brown_g = (as.data.frame(brown_genes, header = F))
x = geneModuleMembership[match(brown_g$V1, rownames(geneModuleMembership)),]
y = geneTraitSignificance[match(brown_g$V1, rownames(geneTraitSignificance)),]
scatter_genes = data.frame(cbind(abs(x$MMbrown),abs(y)))


ggplot(scatter_data, aes(x=X1, y=X2)) +
  geom_point(size=2, shape = 16, alpha = 1/3) +
  xlab("Module membership for Brown module") +
  ylab("Gene significance for cancer") +
  geom_point(mapping = aes(X1, X2) ,data = scatter_genes, colour = "orange", shape = 16)

#####################################################################################