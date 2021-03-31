if (!require("pacman")) install.packages("pacman")
pacman::p_load(WGCNA, ABAData, httr, readxl, tidyverse, patchwork, ggplotify, janitor, ggrepel, ggtext, glue, dplyr)

df <- read_csv("Pardinas.csv")

all_sz_genes <- df$`Gene(s) tagged` %>% 
  str_split(",") %>% 
  unlist() %>% 
  str_trim() %>% 
  as.data.frame() %>% 
  distinct()
names(all_sz_genes) <- "genes"

data(dataset_5_stages)

sz_stage_1 <-
  dataset_5_stages %>%
  dplyr::filter(hgnc_symbol %in% all_sz_genes$genes) %>% 
  dplyr::filter(age_category == 1) %>% 
  dplyr::select(-c(entrezgene, ensembl_gene_id, age_category)) %>% 
  pivot_wider(names_from = structure, values_from = signal) %>% 
  column_to_rownames(var = "hgnc_symbol")%>%
  t() %>% 
  scale() %>% 
  t()


kOne <- sz_stage_1 %>% t()

powersOne = c(c(1:10), seq(from = 12, to=20, by=2))

sftOne = pickSoftThreshold(kOne, powerVector = powersOne, verbose = 5)

cexOne = 0.9;

pOne = sftOne$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = -sign(slope)*sft_r_sq, label = powersOne)) +
  labs(x = "Soft Threshold (power)",
       y = "Scale Free Topology Model Fit, signed R<sup>2</sup>",
       title = "Scale independence") +
  geom_point() +
  geom_text_repel() +
  theme(axis.title.y = element_textbox_simple(
    orientation = "left-rotated",
    width = NULL)) +
  geom_hline(yintercept = 0.77, col = "red" )

p2 = sftOne$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = mean_k, label = powersOne)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  geom_point() +
  geom_text_repel()

pOne + p2

dev.off()
softPowerOne = 7
adjacencyOne = adjacency(kOne, power = softPowerOne);

# Turn adjacency into topological overlap 
TOMOne = TOMsimilarity(adjacencyOne); 

dissTOMOne = 1-TOMOne

# Call the hierarchical clustering function 
geneTreeOne = hclust(as.dist(dissTOMOne), method = "average")
# Plot the resulting clustering tree (dendrogram) 
plot(geneTreeOne, xlab="", sub="", main = glue("Gene clustering on TOM-based dissimilarity\nfor schizophrenia genes during stage One"), labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = 10; 
# Module identification using dynamic tree cut: 
dynamicModsOne = cutreeDynamic(dendro = geneTreeOne, 
                               distM = dissTOMOne, 
                               deepSplit = 2, 
                               pamRespectsDendro = FALSE, 
                               minClusterSize = minModuleSize,
                               cutHeight = 0.95); 

table(dynamicModsOne)

dynamicColorsOne = labels2colors(dynamicModsOne) 
table(dynamicColorsOne) 

plotDendroAndColors(geneTreeOne, 
                    dynamicColorsOne, 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = glue("Age category One colour dendrogram"))

MEListOne = moduleEigengenes(sz_stage_1 %>% t(), colors = dynamicModsOne)


MEsOne = MEListOne$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissOne = 1-cor(MEsOne); 
# Cluster module eigengenes 
METreeOne = hclust(as.dist(MEDissOne), method = "average"); 
# Plot the result 
plot(METreeOne, main = glue("Clustering of module eigengenes for stage One"), xlab = "", sub = "")
MEDissThresOne = 0.4
# Plot the cut line into the dendrogram
abline(h=MEDissThresOne, col = "red")

mergeOne = mergeCloseModules(kOne, dynamicColorsOne, cutHeight = MEDissThresOne, verbose = 3)

mergedColorsOne = mergeOne$colors;
# Eigengenes of the new merged modules:
mergedMEsOne = mergeOne$newMEs

plot(geneTreeOne, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 0.5, 
     cex.axis = 1.5, 
     cex.main = 2, 
     cex = 0.6)

plotDendroAndColors(geneTreeOne,
                    main = "Developmental Stage One Merged Cluster Dendrogram",
                    cbind(dynamicColorsOne, mergedColorsOne),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Tue May 12 11:35:46 2020 ------------------------------
MEListOne = moduleEigengenes(sz_stage_1 %>% t(), colors = dynamicModsOne)
# Tue May 12 11:35:52 2020 ------------------------------

#MEDissThresOne = 1
# Plot the cut line into the dendrogram
#abline(h=MEDissThresOne, col = "red")

mergeOne = mergeCloseModules(kOne, dynamicColorsOne, cutHeight = MEDissThresOne, verbose = 3)

mergedColorsOne = mergeOne$colors;
geneInfoALLOne <- data.frame(geneIDs = colnames(kOne), 
                             moduleColor=mergedColorsOne)

mergedMEsOne = mergeOne$newMEs
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreeOne, cbind(dynamicColorsOne, mergedColorsOne),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Rename to moduleColors
moduleColorsOne = mergedColorsOne
# Construct numerical labels corresponding to the colors
colorOrderone = c("grey", standardColors(50));
moduleLabelsOne = match(moduleColorsOne, colorOrderone)-1;
MEsOne = mergedMEsOne;

MEsOne1 = moduleEigengenes(kOne, moduleColorsOne)$eigengenes
MEsOne11 = orderMEs(MEsOne1)
#View(MEs50)
MEs1One = MEListOne$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissOne = 1-cor(MEs1One); 
# Cluster module eigengenes 
METreeOne = hclust(as.dist(MEDissOne), method = "average"); 
# Plot the result 
plot(METreeOne, main = glue("Clustering of module eigengenes for age category  One"), xlab = "", sub = "")

dissTOMOne = 1-TOMsimilarityFromExpr(kOne, power =softPowerOne);

plotTOMOne = dissTOMOne^softPowerOne;
# Set diagonal to NA for a nicer plot
diag(plotTOMOne) = NA;
# Call the plot function
TOMplot(plotTOMOne, geneTreeOne, moduleColorsOne, main = glue("Network heatmap plot, all sz genes, stage One"))

par(cex = 1.0)
if(ncol(MEsOne) > 2) plotEigengeneNetworks(MEsOne, glue("Eigengene dendrogram for stage One"), marDendro = c(0,4,2,0), plotHeatmaps = FALSE)

par(cex = 1.0)
if(ncol(MEsOne) > 2) plotEigengeneNetworks(MEsOne, glue("Eigengene adjacency heatmap for Developmental Stage  One"), marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

moduleEigengenes(kOne,
                 moduleColorsOne,
                 impute = TRUE,
                 nPC = 1,
                 align = "along average",
                 excludeGrey = FALSE,
                 grey = if (is.numeric(moduleColorsOne)) 0 else "grey",
                 subHubs =TRUE,
                 softPower =softPowerOne,
                 scale = TRUE,
                 verbose = 0, indent = 0)

chooseTopHubInEachModule(
  kOne,
  moduleColorsOne,
  omitColors = "grey",
  power = softPowerOne,
  type = "unsigned")

geneModuleMembershipOne = as.data.frame(cor(kOne, MEsOne11, use = "p"))

## New

setwd("C:/Users/elena/Documents")

library(anRichment)
library(org.Hs.eg.db)
GOcollection = buildGOcollection(organism = "human")


allLLIDs<-read_csv("ABALocusLinked26-03.csv")
GOenrOne = enrichmentAnalysis(
  classLabels = moduleColorsOne, identifiers = allLLIDs$LOCUSLINK_ID,
  refCollection = GOcollection,
  useBackground = "allOrgGenes",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")


GO_per_setOne <- cbind(GOenrOne$enrichmentTable$class, GOenrOne$enrichmentTable$dataSetID,  GOenrOne$enrichmentTable$dataSetName, GOenrOne$enrichmentTable$FDR)
GO_per_set_Onecont <- cbind(GOenrOne$enrichmentTable$class, GOenrOne$enrichmentTable$dataSetID,  GOenrOne$enrichmentTable$dataSetName, GOenrOne$enrichmentTable$FDR,GOenrOne$enrichmentTable$overlapGenes)
tail(GO_per_setOne)

colnames(GO_per_setOne) <- c("Module", "GO term", "GO process", "FDR")
colnames(GO_per_set_Onecont) <- c("Module", "GO Term", "GO Process", "FDR", "Genes")

write.csv(GO_per_set_Onecont, file="GO_per_setOnecont.csv")
# Pull the top (most significant) for each module

top_GO_per_moduleOne <- GO_per_setOne[!duplicated(GO_per_setOne[,1]), ]
colnames(top_GO_per_moduleOne) <- c("Module", "GO term", "GO process", "FDR")

# Recalculate module eigengenes
MEsOne = moduleEigengenes(kOne, moduleColorsOne)$eigengenes

plotEigengeneNetworks(MEsOne, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 

write.csv(GO_per_setOne, file="GO_per_setOne.csv")

GO_per_setOne <- as.data.frame(GO_per_setOne)
GO_per_setOnecont <- as.data.frame(GO_per_set_Onecont)
#Black
GO_per_setOnecont_Black <- GO_per_setOnecont[GO_per_setOnecont$Module == "black",]
write.csv(GO_per_setOnecont_Black, file="GO_per_setOne_black.csv")
#BlUE
GO_per_setOnecont_Blue <- GO_per_setOnecont[GO_per_setOnecont$Module == "blue",]
write.csv(GO_per_setOnecont_Blue, file="GO_per_setOne_blue.csv")
#Brown
GO_per_setOnecont_Brown <- GO_per_setOnecont[GO_per_setOnecont$Module == "brown",]
write.csv(GO_per_setOnecont_Brown, file="GO_per_setOne_brown.csv")
#Pink
GO_per_setOnecont_Pink <- GO_per_setOnecont[GO_per_setOnecont$Module == "pink",]
write.csv(GO_per_setOnecont_Pink, file="GO_per_setOne_Pink.csv")
#Turquoise
GO_per_setOnecont_Turquoise <- GO_per_setOnecont[GO_per_setOnecont$Module == "turquoise",]
write.csv(GO_per_setOnecont_Turquoise, file="GO_per_setOne_Turquoise.csv")

pdf("plotEigeneNetworks.pdf")

plotEigengeneNetworks(MEsOne, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

blackModule<-subset(geneInfoALL, moduleColor == "black", select = c("geneIDs", "moduleColor"))
write.csv(blackModule, "blackModuleOne.csv")

geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[blackModule$geneIDs,]->geneModuleMembershipOne_black #get modulemembership for genes in blue module
geneModuleMembershipOne_black$MEblack->geneModuleMembership1black #correlation of each gene with the with module eigengene for genes in blue module

kOne[,blackModule$geneIDs] #get expression values for genes in blue module only

cor((kOne[,blackModule$geneIDs]))->cor3 #correlation of each gene expression value with every other gene in the blue module


rownames(cor3)->vec3 # list of gene names in blue module
exportNetworkToCytoscape(cor3, altNodeNames =vec3, nodeAttr=geneModuleMembershipOne_black)->cytoblack1
write.csv(cytoblack1$edgeData, "cyto_blackEDGE1.csv",quote=FALSE)
write.csv(cytoblack1$nodeData, "cyto_blackNODE1.csv",quote=FALSE)
save(cytoblack1, file="cytoblack1.Rdata")

pacman::p_load(tidyverse, httr, readxl, janitor,dplyr)

pf <- df %>%
  clean_names() %>%
  separate_rows(gene_s_tagged) %>%
  dplyr::select(genes = gene_s_tagged, p_value) %>%
  distinct() %>% 
  group_by(genes) %>% 
  summarise(p_value = min(p_value, na.rm = T)) %>% 
  ungroup()


datablack1 <- read.csv("cyto_blackEDGE1.csv")
datablack1$fromNode2 = datablack1$fromNode

names(datablack1)[8] <- "genes"

cytoblack1 <- datablack1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoblack11 = cytoblack1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)
write.csv(cytoblack1, "cyto_blackEDGE1.csv",quote=FALSE)
write.csv(cytoblack11, "cyto_blackEDGE11.csv",quote=FALSE)


#Blue
BlueModule1<-subset(geneInfoALL, moduleColor == "blue", select = c("geneIDs", "moduleColor"))
write.csv(BlueModule1, "BlueModule1.csv")

geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[BlueModule1$geneIDs,]->geneModuleMembership1_Blue#get modulemembership for genes in blue module
geneModuleMembership1_Blue$MEBlue->geneModuleMembershipBlue #correlation of each gene with the with module eigengene for genes in blue module

kOne[,BlueModule1$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,BlueModule1$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_Blue)->cytoBlue1

write.csv(cytoBlue1$edgeData, "cyto_BlueEDGE1.csv",quote=FALSE)
write.csv(cytoBlue1$nodeData, "cyto_BlueNODE1.csv",quote=FALSE)
save(cytoBlue1, file="cytoBlue1.Rdata")

BlueModule<-subset(geneInfoALL, moduleColor == "blue", select = c("geneIDs", "moduleColor"))
write.csv(BlueModule, "BlueModule.csv")

dataBlue1 <- read.csv("cyto_BlueEDGE1.csv")

dataBlue1$fromNode2 = dataBlue1$fromNode

names(dataBlue1)[8] <- "genes"

cytoBlue1 <- dataBlue1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytoBlue1, "cyto_BlueEDGE1.csv",quote=FALSE)

cytoBlue11 = cytoBlue1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoBlue11, "cyto_BlueEDGE11.csv",quote=FALSE)
#Brown
geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

brownModule<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(brownModule, "brownModule.csv")
geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[brownModule$geneIDs,]->geneModuleMembership1_brown#get modulemembership for genes in blue module
geneModuleMembership1_brown$MEbrown->geneModuleMembershipbrown #correlation of each gene with the with module eigengene for genes in blue module

kOne[,brownModule$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,brownModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_brown)->cytobrown1

write.csv(cytobrown1$edgeData, "cyto_brownEDGE1.csv",quote=FALSE)
write.csv(cytobrown1$nodeData, "cyto_brownNODE1.csv",quote=FALSE)
save(cytobrown1, file="cytobrown1.Rdata")

BrownModule1<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(BrownModule1, "BrownModule1.csv")

databrown1 <- read.csv("cyto_brownEDGE1.csv")

databrown1$fromNode2 = databrown1$fromNode

names(databrown1)[8] <- "genes"

cytobrown1 <- databrown1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytobrown1, "cyto_brownEDGE1.csv",quote=FALSE)


cytobrown11 = cytobrown1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytobrown11, "cyto_brownEDGE11.csv",quote=FALSE)
#Pink
geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

PinkModule<-subset(geneInfoALL, moduleColor == "pink", select = c("geneIDs", "moduleColor"))
write.csv(PinkModule, "PinkModule.csv")
geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[PinkModule$geneIDs,]->geneModuleMembership1_Pink#get modulemembership for genes in blue module
geneModuleMembership1_Pink$MEPink->geneModuleMembershipPink #correlation of each gene with the with module eigengene for genes in blue module

kOne[,PinkModule$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,PinkModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_Pink)->cytoPink1

write.csv(cytoPink1$edgeData, "cyto_PinkEDGE1.csv",quote=FALSE)
write.csv(cytoPink1$nodeData, "cyto_PinkNODE1.csv",quote=FALSE)
save(cytoPink1, file="cytoPink1.Rdata")

dataPink1 <- read.csv("cyto_PinkEDGE1.csv")

dataPink1$fromNode2 = dataPink1$fromNode

names(dataPink1)[8] <- "genes"

cytoPink1 <- dataPink1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoPink1, "cyto_PinkEDGE1.csv",quote=FALSE)

cytoPink11 = cytoPink1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoPink11, "cyto_PinkEDGE11.csv",quote=FALSE)

#Turq
geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

TurqModule<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(TurqModule, "TurqModule.csv")
geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[TurqModule$geneIDs,]->geneModuleMembership1_Turq#get modulemembership for genes in blue module
geneModuleMembership1_Turq$METurq->geneModuleMembershipTurq #correlation of each gene with the with module eigengene for genes in blue module

kOne[,TurqModule$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,TurqModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_Turq)->cytoTurq1

write.csv(cytoTurq1$edgeData, "cyto_TurqEDGE1.csv",quote=FALSE)
write.csv(cytoTurq1$nodeData, "cyto_TurqNODE1.csv",quote=FALSE)
save(cytoTurq1, file="cytoTurq1.Rdata")

dataTurq1 <- read.csv("cyto_TurqEDGE1.csv")

dataTurq1$fromNode2 = dataTurq1$fromNode

names(dataTurq1)[8] <- "genes"

cytoTurq1 <- dataTurq1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoTurq1, "cyto_TurqEDGE1.csv",quote=FALSE)

cytoTurq11 = cytoTurq1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoTurq11, "cyto_TurqEDGE11.csv",quote=FALSE)

