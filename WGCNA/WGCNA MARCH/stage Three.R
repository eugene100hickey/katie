
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ABAData, httr, readxl, tidyverse, patchwork, ggplotify, janitor, ggrepel, ggtext, glue, dplyr, WGCNA)


df <- read_csv("Pardinas.csv")

all_sz_genes <- df$`Gene(s) tagged` %>% 
  str_split(",") %>% 
  unlist() %>% 
  str_trim() %>% 
  as.data.frame() %>% 
  distinct()
names(all_sz_genes) <- "genes"

data(dataset_5_stages)



sz_stage_3 <-
  dataset_5_stages %>%
  dplyr::filter(hgnc_symbol %in% all_sz_genes$genes) %>% 
  dplyr::filter(age_category == 3) %>% 
  dplyr::select(-c(entrezgene, ensembl_gene_id, age_category)) %>% 
  pivot_wider(names_from = structure, values_from = signal) %>% 
  column_to_rownames(var = "hgnc_symbol")%>%
  t() %>% 
  scale() %>% 
  t()


kThree <- sz_stage_3 %>% t()

powersThree = c(c(1:10), seq(from = 12, to=20, by=2))

sftThree = pickSoftThreshold(kThree, powerVector = powersThree, verbose = 5)

cexOne = 0.9;

pThree = sftThree$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = -sign(slope)*sft_r_sq, label = powersThree)) +
  labs(x = "Soft Threshold (power)",
       y = "Scale Free Topology Model Fit, signed R<sup>2</sup>",
       title = "Scale independence") +
  geom_point() +
  geom_text_repel() +
  theme(axis.title.y = element_textbox_simple(
    orientation = "left-rotated",
    width = NULL)) +
  geom_hline(yintercept = 0.83, col = "red" )

pThree2 = sftThree$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = mean_k, label = powersThree)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  geom_point() +
  geom_text_repel()

pThree + pThree2
dev.off()	
softPowerThree = 6
adjacencyThree = adjacency(kThree, power = softPowerThree);

# Turn adjacency into topological overlap 
TOMThree = TOMsimilarity(adjacencyThree); 

dissTOMThree = 1-TOMThree

# Call the hierarchical clustering function 
geneTreeThree = hclust(as.dist(dissTOMThree), method = "average")
# Plot the resulting clustering tree (dendrogram) 
plot(geneTreeThree, xlab="", sub="", main = glue("Gene clustering on TOM-based dissimilarity\nfor schizophrenia genes during stage Three"), labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = 10; 
# Module identification using dynamic tree cut: 
dynamicModsThree = cutreeDynamic(dendro = geneTreeThree, 
                                 distM = dissTOMThree, 
                                 deepSplit = 2, 
                                 pamRespectsDendro = FALSE, 
                                 minClusterSize = minModuleSize,
                                 cutHeight = 0.95); 

table(dynamicModsThree)

dynamicColorsThree = labels2colors(dynamicModsThree) 
table(dynamicColorsThree) 

plotDendroAndColors(geneTreeThree, 
                    dynamicColorsThree, 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = glue("Age category Three colour dendrogram"))

MEListThree = moduleEigengenes(sz_stage_3 %>% t(), colors = dynamicModsThree)


MEsThree = MEListThree$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissThree = 1-cor(MEsThree); 
# Cluster module eigengenes 
METreeThree = hclust(as.dist(MEDissThree), method = "average"); 
# Plot the result 
plot(METreeThree, main = glue("Clustering of module eigengenes for stage Three"), xlab = "", sub = "")
MEDissThresThree = 0.65
# Plot the cut line into the dendrogram
abline(h=MEDissThresThree, col = "red")

mergeThree = mergeCloseModules(kThree, dynamicColorsThree, cutHeight = MEDissThresThree, verbose = 3)

mergedColorsThree = mergeThree$colors;
# Eigengenes of the new merged modules:
mergedMEsThree = mergeThree$newMEs

plot(geneTreeThree, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 0.5, 
     cex.axis = 1.5, 
     cex.main = 2, 
     cex = 0.6)

plotDendroAndColors(geneTreeThree, 
                    main = "Developmental Stage Three Cluster Dendrogram",
                    cbind(dynamicColorsThree, mergedColorsThree),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



MEListThree = moduleEigengenes(sz_stage_3 %>% t(), colors = dynamicModsThree)


#MEDissThresOne = 1
# Plot the cut line into the dendrogram
#abline(h=MEDissThresOne, col = "red")

mergeThree = mergeCloseModules(kThree, dynamicColorsThree, cutHeight = MEDissThresThree, verbose = 3)

mergedColorsThree = mergeThree$colors;
geneInfoALLThree <- data.frame(geneIDs = colnames(kThree), 
                               moduleColor=mergedColorsThree)

mergedMEsThree = mergeThree$newMEs
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreeThree, cbind(dynamicColorsThree, mergedColorsThree),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()
# Rename to moduleColors
moduleColorsThree = mergedColorsThree
# Construct numerical labels corresponding to the colors
colorOrderThree = c("grey", standardColors(50));
moduleLabelsThree = match(moduleColorsThree, colorOrderThree)-1;
MEsThree = mergedMEsThree;

MEsOne3 = moduleEigengenes(kThree, moduleColorsThree)$eigengenes
MEsOne33 = orderMEs(MEsOne3)
#View(MEs50)
MEs3Three = MEListThree$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissThree = 1-cor(MEs3Three); 
# Cluster module eigengenes 
METreeThree = hclust(as.dist(MEDissThree), method = "average"); 
# Plot the result 
plot(METreeThree, main = glue("Clustering of module eigengenes for age category  Three"), xlab = "", sub = "")

dissTOMThree = 1-TOMsimilarityFromExpr(kThree, power =softPowerThree);

plotTOMThree = dissTOMThree^softPowerThree;
# Set diagonal to NA for a nicer plot
diag(plotTOMThree) = NA;
# Call the plot function
TOMplot(plotTOMThree, geneTreeThree, moduleColorsThree, main = glue("Network heatmap plot, all sz genes, stage Three"))

par(cex = 1.0)
if(ncol(MEsThree) > 2) plotEigengeneNetworks(MEsThree, glue("Eigengene dendrogram for stage Three"), marDendro = c(0,4,2,0), plotHeatmaps = FALSE)

par(cex = 1.0)
if(ncol(MEsThree) > 2) plotEigengeneNetworks(MEsThree, glue("Eigengene adjacency heatmap for Stage Three"), marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

moduleEigengenes(kThree,
                 moduleColorsThree,
                 impute = TRUE,
                 nPC = 1,
                 align = "along average",
                 excludeGrey = FALSE,
                 grey = if (is.numeric(moduleColorsThree)) 0 else "grey",
                 subHubs =TRUE,
                 softPower =softPowerThree,
                 scale = TRUE,
                 verbose = 0, indent = 0)

chooseTopHubInEachModule(
  kThree,
  moduleColorsThree,
  omitColors = "grey",
  power = softPowerThree,
  type = "unsigned")

geneModuleMembershipOne = as.data.frame(cor(kThree, MEsOne33, use = "p"))

## New

setwd("C:/Users/elena/Documents")

library(anRichment)
library(org.Hs.eg.db)
GOcollection = buildGOcollection(organism = "human")


allLLIDs<-read_csv("ABALocusLinked26-03.csv")
GOenrThree = enrichmentAnalysis(
  classLabels = moduleColorsThree, identifiers = allLLIDs$LOCUSLINK_ID,
  refCollection = GOcollection,
  useBackground = "allOrgGenes",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")


GO_per_setThree <- cbind(GOenrThree$enrichmentTable$class, GOenrThree$enrichmentTable$dataSetID,  GOenrThree$enrichmentTable$dataSetName, GOenrThree$enrichmentTable$FDR)
GO_per_set_Threecont <- cbind(GOenrThree$enrichmentTable$class, GOenrThree$enrichmentTable$dataSetID,  GOenrThree$enrichmentTable$dataSetName, GOenrThree$enrichmentTable$FDR,GOenrThree$enrichmentTable$overlapGenes)
tail(GO_per_setThree)

colnames(GO_per_setThree) <- c("Module", "GO term", "GO process", "FDR")
colnames(GO_per_set_Threecont) <- c("Module", "GO Term", "GO Process", "FDR", "Genes")

write.csv(GO_per_set_Threecont, file="GO_per_setThreecont.csv")
# Pull the top (most significant) for each module

top_GO_per_moduleThree <- GO_per_setThree[!duplicated(GO_per_setThree[,1]), ]
colnames(top_GO_per_moduleThree) <- c("Module", "GO term", "GO process", "FDR")

# Recalculate module eigengenes
MEsThree = moduleEigengenes(kThree, moduleColorsThree)$eigengenes

plotEigengeneNetworks(MEsThree, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 

write.csv(GO_per_setThree, file="GO_per_setThree.csv")

GO_per_setThree <- as.data.frame(GO_per_setThree)
GO_per_setThreecont <- as.data.frame(GO_per_set_Threecont)


#Brown
GO_per_setThreecont_Brown <- GO_per_setThreecont[GO_per_setThreecont$Module == "brown",]
write.csv(GO_per_setThreecont_Brown, file="GO_per_setThree_brown.csv")
#Turq
GO_per_setThreecont_Turquoise <- GO_per_setThreecont[GO_per_setThreecont$Module == "turquoise",]
write.csv(GO_per_setThreecont_Turquoise, file="GO_per_setThree_Turquoise.csv")
#Yellow
GO_per_setThreecont_Yell <- GO_per_setThreecont[GO_per_setThreecont$Module == "yellow",]
write.csv(GO_per_setThreecont_Yell, file="GO_per_setThree_Yellow.csv")



pdf("plotEigeneNetworks.pdf")
plotEigengeneNetworks(MEsThree, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(kThree), moduleColor=moduleColorsThree)




#Turq
TurquoiseModule3<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(TurquoiseModule3, "TurquoiseModule3.csv")

geneModuleMembership3 = as.data.frame(cor(kThree, MEsOne33, use = "p")) #module membership for all genes all modules
geneModuleMembership3[TurquoiseModule3$geneIDs,]->geneModuleMembership3_Turquoise#get modulemembership for genes in blue module
geneModuleMembership3_Turquoise$METurquoise->geneModuleMembershipTurquoise #correlation of each gene with the with module eigengene for genes in blue module

kThree[,TurquoiseModule3$geneIDs] #get expression values for genes in turquoise module only

cor((kThree[,TurquoiseModule3$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership3_Turquoise)->cytoTurquoise3


write.csv(cytoTurquoise3$edgeData, "cyto_TurquoiseEDGE3.csv",quote=FALSE)
write.csv(cytoTurquoise3$nodeData, "cyto_TurquoiseNODE3.csv",quote=FALSE)
save(cytoTurquoise3, file="cytoTurquoise3.Rdata")

TurquoiseModule<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(TurquoiseModule, "TurquoiseModule3.csv")

dataTurquoise3 <- read.csv("cyto_TurquoiseEDGE3.csv")

dataTurquoise3$fromNode2 = dataTurquoise3$fromNode

names(dataTurquoise3)[8] <- "genes"

cytoTurquoise3 <- dataTurquoise3 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoTurquoise3, "cyto_Turquoise3.csv",quote=FALSE)

cytoTurquoise33= cytoTurquoise3 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoTurquoise33, "cyto_TurquoiseEDGE33.csv",quote=FALSE)
#Brown
geneInfoALL<-data.frame(geneIDs = colnames(kThree), moduleColor=moduleColorsThree)

brownModule<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(brownModule, "brownModule3.csv")
geneModuleMembership3 = as.data.frame(cor(kThree, MEsOne33, use = "p")) #module membership for all genes all modules
geneModuleMembership3[brownModule$geneIDs,]->geneModuleMembership3_brown#get modulemembership for genes in blue module
geneModuleMembership3_brown$MEbrown->geneModuleMembershipbrown #correlation of each gene with the with module eigengene for genes in blue module

kThree[,brownModule$geneIDs] #get expression values for genes in turquoise module only

cor((kThree[,brownModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership3_brown)->cytobrown3


write.csv(cytobrown3$edgeData, "cyto_brownEDGE3.csv",quote=FALSE)
write.csv(cytobrown3$nodeData, "cyto_brownNODE3.csv",quote=FALSE)
save(cytobrown3, file="cytobrown3.Rdata")

databrown3 <- read.csv("cyto_brownEDGE3.csv")

databrown3$fromNode2 = databrown3$fromNode

names(databrown3)[8] <- "genes"

cytobrown3 <- databrown3 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytobrown3, "cyto_brownEDGE3.csv",quote=FALSE)

cytobrown33 = cytobrown3 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytobrown33, "cyto_brownEDGE33.csv",quote=FALSE)

#Yellow
geneInfoALL<-data.frame(geneIDs = colnames(kThree), moduleColor=moduleColorsThree)

YellModule3<-subset(geneInfoALL, moduleColor == "yellow", select = c("geneIDs", "moduleColor"))
write.csv(YellModule3, "YellModule3.csv")
geneModuleMembership3 = as.data.frame(cor(kThree, MEsOne33, use = "p")) #module membership for all genes all modules
geneModuleMembership3[YellModule3$geneIDs,]->geneModuleMembership3_Yell#get modulemembership for genes in blue module
geneModuleMembership3_Yell$MEYell->geneModuleMembershipYell #correlation of each gene with the with module eigengene for genes in blue module

kThree[,YellModule3$geneIDs] #get expression values for genes in turquoise module only

cor((kThree[,YellModule3$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership3_Yell)->cytoYell3


write.csv(cytoYell3$edgeData, "cyto_YellEDGE3.csv",quote=FALSE)
write.csv(cytoYell3$nodeData, "cyto_YellNODE3.csv",quote=FALSE)
save(cytoYell3, file="cytoYell3.Rdata")

dataYell3 <- read.csv("cyto_YellEDGE3.csv")

dataYell3$fromNode2 = dataYell3$fromNode

names(dataYell3)[8] <- "genes"

cytoYell3 <- dataYell3 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytoYell3, "cyto_YellEDGE3.csv",quote=FALSE)
cytoYell33 = cytoYell3 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoYell33, "cyto_YellEDGE33.csv",quote=FALSE)

