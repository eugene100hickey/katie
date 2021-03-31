
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



sz_stage_4 <-
  dataset_5_stages %>%
  dplyr::filter(hgnc_symbol %in% all_sz_genes$genes) %>% 
  dplyr::filter(age_category == 4) %>% 
  dplyr::select(-c(entrezgene, ensembl_gene_id, age_category)) %>% 
  pivot_wider(names_from = structure, values_from = signal) %>% 
  column_to_rownames(var = "hgnc_symbol")%>%
  t() %>% 
  scale() %>% 
  t()


kFour <- sz_stage_4 %>% t()

powersFour = c(c(1:10), seq(from = 12, to=20, by=2))

sftFour = pickSoftThreshold(kFour, powerVector = powersFour, verbose = 5)

cexOne = 0.9;

pFour = sftFour$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = -sign(slope)*sft_r_sq, label = powersFour)) +
  labs(x = "Soft Threshold (power)",
       y = "Scale Free Topology Model Fit, signed R<sup>2</sup>",
       title = "Scale independence") +
  geom_point() +
  geom_text_repel() +
  theme(axis.title.y = element_textbox_simple(
    orientation = "left-rotated",
    width = NULL)) +
  geom_hline(yintercept = 0.85, col = "red" )

pFour2 = sftFour$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = mean_k, label = powersFour)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  geom_point() +
  geom_text_repel()

pFour + pFour2

softPowerFour = 9
adjacencyFour = adjacency(kFour, power = softPowerFour);

# Turn adjacency into topological overlap 
TOMFour = TOMsimilarity(adjacencyFour); 

dissTOMFour = 1-TOMFour

# Call the hierarchical clustering function 
geneTreeFour = hclust(as.dist(dissTOMFour), method = "average")
# Plot the resulting clustering tree (dendrogram) 
plot(geneTreeFour, xlab="", sub="", main = glue("Gene clustering on TOM-based dissimilarity\nfor schizophrenia genes during stage Four"), labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = 10; 
# Module identification using dynamic tree cut: 
dynamicModsFour = cutreeDynamic(dendro = geneTreeFour, 
                                distM = dissTOMFour, 
                                deepSplit = 2, 
                                pamRespectsDendro = FALSE, 
                                minClusterSize = minModuleSize,
                                cutHeight = 0.95); 

table(dynamicModsFour)

dynamicColorsFour = labels2colors(dynamicModsFour) 
table(dynamicColorsFour) 

plotDendroAndColors(geneTreeFour, 
                    dynamicColorsFour, 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = glue("Age category Four colour dendrogram"))

MEListFour = moduleEigengenes(sz_stage_4 %>% t(), colors = dynamicModsFour)


MEsFour = MEListFour$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissFour = 1-cor(MEsFour); 
# Cluster module eigengenes 
METreeFour = hclust(as.dist(MEDissFour), method = "average"); 
# Plot the result 
plot(METreeFour, main = glue("Clustering of module eigengenes for stage Four"), xlab = "", sub = "")
MEDissThresFour = 0.35
# Plot the cut line into the dendrogram
abline(h=MEDissThresFour, col = "red")

mergeFour = mergeCloseModules(kFour, dynamicColorsFour, cutHeight = MEDissThresFour, verbose = 3)

mergedColorsFour = mergeFour$colors;
# Eigengenes of the new merged modules:
mergedMEsFour = mergeFour$newMEs

plot(geneTreeFour, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 0.5, 
     cex.axis = 1.5, 
     cex.main = 2, 
     cex = 0.6)

plotDendroAndColors(geneTreeFour, 
                    main = "Developmental Stage Four Cluster Dendrogram",
                    cbind(dynamicColorsFour, mergedColorsFour),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



MEListFour = moduleEigengenes(sz_stage_4 %>% t(), colors = dynamicModsFour)


#MEDissThresOne = 1
# Plot the cut line into the dendrogram
#abline(h=MEDissThresOne, col = "red")

mergeFour = mergeCloseModules(kFour, dynamicColorsFour, cutHeight = MEDissThresFour, verbose = 3)

mergedColorsFour = mergeFour$colors;
geneInfoALLFour <- data.frame(geneIDs = colnames(kFour), 
                              moduleColor=mergedColorsFour)

mergedMEsFour = mergeFour$newMEs
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreeFour, cbind(dynamicColorsFour, mergedColorsFour),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Rename to moduleColors
moduleColorsFour = mergedColorsFour
# Construct numerical labels corresponding to the colors
colorOrderFour = c("grey", standardColors(50));
moduleLabelsFour = match(moduleColorsFour, colorOrderFour)-1;
MEsFour = mergedMEsFour;

MEsOne4 = moduleEigengenes(kFour, moduleColorsFour)$eigengenes
MEsOne44 = orderMEs(MEsOne4)
#View(MEs50)
MEs4Four = MEListFour$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissFour = 1-cor(MEs4Four); 
# Cluster module eigengenes 
METreeFour = hclust(as.dist(MEDissFour), method = "average"); 
# Plot the result 
plot(METreeFour, main = glue("Clustering of module eigengenes for age category  Four"), xlab = "", sub = "")

dissTOMFour = 1-TOMsimilarityFromExpr(kFour, power =softPowerFour);

plotTOMFour = dissTOMFour^softPowerFour;
# Set diagonal to NA for a nicer plot
diag(plotTOMFour) = NA;
# Call the plot function
TOMplot(plotTOMFour, geneTreeFour, moduleColorsFour, main = glue("Network heatmap plot, all sz genes, stage Four"))

par(cex = 1.0)
if(ncol(MEsFour) > 2) plotEigengeneNetworks(MEsFour, glue("Eigengene dendrogram for stage Four"), marDendro = c(0,4,2,0), plotHeatmaps = FALSE)

par(cex = 1.0)
if(ncol(MEsFour) > 2) plotEigengeneNetworks(MEsFour, glue("Eigengene adjacency heatmap for stage  Four"), marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

moduleEigengenes(kFour,
                 moduleColorsFour,
                 impute = TRUE,
                 nPC = 1,
                 align = "along average",
                 excludeGrey = FALSE,
                 grey = if (is.numeric(moduleColorsTwo)) 0 else "grey",
                 subHubs =TRUE,
                 softPower =softPowerFour,
                 scale = TRUE,
                 verbose = 0, indent = 0)

chooseTopHubInEachModule(
  kFour,
  moduleColorsFour,
  omitColors = "grey",
  power = softPowerFour,
  type = "unsigned")

geneModuleMembershipFour = as.data.frame(cor(kFour, MEsOne44, use = "p"))

## New

setwd("C:/Users/elena/Documents")

library(anRichment)
library(org.Hs.eg.db)
GOcollection = buildGOcollection(organism = "human")


allLLIDs<-read_csv("ABALocusLinked26-03.csv")
GOenrFour = enrichmentAnalysis(
  classLabels = moduleColorsFour, identifiers = allLLIDs$LOCUSLINK_ID,
  refCollection = GOcollection,
  useBackground = "allOrgGenes",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")


GO_per_setFour <- cbind(GOenrFour$enrichmentTable$class, GOenrFour$enrichmentTable$dataSetID,  GOenrFour$enrichmentTable$dataSetName, GOenrFour$enrichmentTable$FDR)
GO_per_set_Fourcont <- cbind(GOenrFour$enrichmentTable$class, GOenrFour$enrichmentTable$dataSetID,  GOenrFour$enrichmentTable$dataSetName, GOenrFour$enrichmentTable$FDR,GOenrFour$enrichmentTable$overlapGenes)
tail(GO_per_setFour)

colnames(GO_per_setFour) <- c("Module", "GO term", "GO process", "FDR")
colnames(GO_per_set_Fourcont) <- c("Module", "GO Term", "GO Process", "FDR", "Genes")

write.csv(GO_per_set_Fourcont, file="GO_per_setFourcont.csv")
# Pull the top (most significant) for each module

top_GO_per_moduleFour <- GO_per_setFour[!duplicated(GO_per_setFour[,1]), ]
colnames(top_GO_per_moduleFour) <- c("Module", "GO term", "GO process", "FDR")

# Recalculate module eigengenes
MEsFour = moduleEigengenes(kFour, moduleColorsFour)$eigengenes

plotEigengeneNetworks(MEsFour, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 

write.csv(GO_per_setFour, file="GO_per_setFour.csv")

GO_per_setFour <- as.data.frame(GO_per_setFour)
GO_per_setFourcont <- as.data.frame(GO_per_set_Fourcont)
#Blue
GO_per_setFourcont_Blue <- GO_per_setFourcont[GO_per_setFourcont$Module == "blue",]
write.csv(GO_per_setFourcont_Blue, file="GO_per_setFour_blue.csv")
#Brown
GO_per_setFourcont_Brown <- GO_per_setFourcont[GO_per_setFourcont$Module == "brown",]
write.csv(GO_per_setFourcont_Brown, file="GO_per_setFour_brown.csv")
#Green
GO_per_setFourcont_Green <- GO_per_setFourcont[GO_per_setFourcont$Module == "green",]
write.csv(GO_per_setFourcont_Green, file="GO_per_setFour_Green.csv")
#Red
GO_per_setFourcont_Red <- GO_per_setFourcont[GO_per_setFourcont$Module == "red",]
write.csv(GO_per_setFourcont_Red, file="GO_per_setFour_Red.csv")
#Magenta
GO_per_setFourcont_Magenta <- GO_per_setFourcont[GO_per_setFourcont$Module == "magenta",]
write.csv(GO_per_setFourcont_Magenta, file="GO_per_setFour_Magenta.csv")
#Purple
GO_per_setFourcont_Purple <- GO_per_setFourcont[GO_per_setFourcont$Module == "purple",]
write.csv(GO_per_setFourcont_Purple, file="GO_per_setFour_Purple.csv")
#Turq
GO_per_setFourcont_Turq <- GO_per_setFourcont[GO_per_setFourcont$Module == "turquoise",]
write.csv(GO_per_setFourcont_Turq, file="GO_per_setFour_Turq.csv")
#Yellow
GO_per_setFourcont_Yellow <- GO_per_setFourcont[GO_per_setFourcont$Module == "yellow",]
write.csv(GO_per_setFourcont_Yellow, file="GO_per_setFour_Yellow.csv")


pdf("plotEigeneNetworks.pdf")
plotEigengeneNetworks(MEsFour, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(kFour), moduleColor=moduleColorsFour)

BlueModule4<-subset(geneInfoALL, moduleColor == "blue", select = c("geneIDs", "moduleColor"))
write.csv(BlueModule4, "BlueModuleFour.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[BlueModule$geneIDs,]->geneModuleMembershipFour_Blue #get modulemembership for genes in blue module
geneModuleMembershipFour_Blue$MEblue->geneModuleMembership4Blue #correlation of each gene with the with module eigengene for genes in blue module

kFour[,BlueModule$geneIDs] #get expression values for genes in blue module only

cor((kFour[,BlueModule$geneIDs]))->cor3 #correlation of each gene expression value with every other gene in the blue module


rownames(cor3)->vec3 # list of gene names in blue module
exportNetworkToCytoscape(cor3, altNodeNames =vec3, nodeAttr=geneModuleMembershipFour_Blue)->cytoBlue4
write.csv(cytoBlue4$edgeData, "cyto_BlueEDGE4.csv",quote=FALSE)
write.csv(cytoBlue4$nodeData, "cyto_BlueNODE4.csv",quote=FALSE)
save(cytoBlue4, file="cytoBlue4.Rdata")


pf <- df %>%
  clean_names() %>%
  separate_rows(gene_s_tagged) %>%
  dplyr::select(genes = gene_s_tagged, p_value) %>%
  distinct() %>% 
  group_by(genes) %>% 
  summarise(p_value = min(p_value, na.rm = T)) %>% 
  ungroup()


dataBlue4 <- read.csv("cyto_BlueEDGE4.csv")

dataBlue4$fromNode2 = dataBlue4$fromNode

names(dataBlue4)[8] <- "genes"

cytoBlue4 <- dataBlue4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoBlue44 = cytoBlue4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)
write.csv(cytoBlue4, "cyto_BlueEDGE4.csv",quote=FALSE)
write.csv(cytoBlue44, "cyto_BlueEDGE44.csv",quote=FALSE)


#Green
GreenModule4<-subset(geneInfoALL, moduleColor == "green", select = c("geneIDs", "moduleColor"))
write.csv(GreenModule4, "GreenModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[GreenModule4$geneIDs,]->geneModuleMembership4_Green#get modulemembership for genes in blue module
geneModuleMembership4_Green$MEGreen->geneModuleMembershipGreen #correlation of each gene with the with module eigengene for genes in blue module

kFour[,GreenModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,GreenModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module

rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_Green)->cytoGreen4

write.csv(cytoGreen4$edgeData, "cyto_GreenEDGE4.csv",quote=FALSE)
write.csv(cytoGreen4$nodeData, "cyto_GreenNODE4.csv",quote=FALSE)
save(cytoGreen4, file="cytoGreen4.Rdata")

GreenModule<-subset(geneInfoALL, moduleColor == "green", select = c("geneIDs", "moduleColor"))
write.csv(GreenModule, "GreenModule4.csv")

dataGreen4 <- read.csv("cyto_GreenEDGE4.csv")

dataGreen4$fromNode2 = dataGreen4$fromNode

names(dataGreen4)[8] <- "genes"

cytoGreen4 <- dataGreen4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 
write.csv(cytoGreen4, "cyto_Green4.csv",quote=FALSE)
cytoGreen44 = cytoGreen4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)
write.csv(cytoGreen44, "cyto_GreenEDGE44.csv",quote=FALSE)

#Brown

BrownModule4<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(BrownModule4, "brownModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[BrownModule4$geneIDs,]->geneModuleMembership4_brown#get modulemembership for genes in blue module
geneModuleMembership4_brown$MEbrown->geneModuleMembershipbrown #correlation of each gene with the with module eigengene for genes in blue module

kFour[,BrownModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,BrownModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_brown)->cytobrown4


write.csv(cytobrown4$edgeData, "cyto_brownEDGE4.csv",quote=FALSE)
write.csv(cytobrown4$nodeData, "cyto_brownNODE4.csv",quote=FALSE)
save(cytobrown4, file="cytobrown4.Rdata")

brownModule<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(brownModule, "brownModule4.csv")

databrown4 <- read.csv("cyto_brownEDGE4.csv")

databrown4$fromNode2 = databrown4$fromNode

names(databrown4)[8] <- "genes"

cytobrown4 <- databrown4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytobrown4, "cyto_brown4.csv",quote=FALSE)

cytobrown44 = cytobrown4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytobrown44, "cyto_brownEDGE44.csv",quote=FALSE)

#Turquoise

TurqModule4<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(TurqModule4, "TurqModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[TurqModule4$geneIDs,]->geneModuleMembership4_Turq#get modulemembership for genes in blue module
geneModuleMembership4_Turq$METurq->geneModuleMembershipTurq #correlation of each gene with the with module eigengene for genes in blue module

kFour[,TurqModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,TurqModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_Turq)->cytoTurq4


write.csv(cytoTurq4$edgeData, "cyto_TurqEDGE4.csv",quote=FALSE)
write.csv(cytoTurq4$nodeData, "cyto_TurqNODE4.csv",quote=FALSE)
save(cytoTurq4, file="cytoTurq4.Rdata")

TurqModule4<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(TurqModule4, "TurqModule4.csv")

dataTurq4 <- read.csv("cyto_TurqEDGE4.csv")

dataTurq4$fromNode2 = dataTurq4$fromNode

names(dataTurq4)[8] <- "genes"

cytoTurq4 <- dataTurq4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoTurq4, "cyto_Turq4.csv",quote=FALSE)

cytoTurq44 = cytoTurq4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoTurq44, "cyto_TurqEDGE44.csv",quote=FALSE)

#Magenta

MagentaModule4<-subset(geneInfoALL, moduleColor == "magenta", select = c("geneIDs", "moduleColor"))
write.csv(MagentaModule4, "MagentaModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[MagentaModule4$geneIDs,]->geneModuleMembership4_Magenta#get modulemembership for genes in blue module
geneModuleMembership4_Magenta$MEMagenta->geneModuleMembershipMagenta #correlation of each gene with the with module eigengene for genes in blue module

kFour[,MagentaModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,MagentaModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module
rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_Magenta)->cytoMagenta4

write.csv(cytoMagenta4$edgeData, "cyto_MagentaEDGE4.csv",quote=FALSE)
write.csv(cytoMagenta4$nodeData, "cyto_MagentaNODE4.csv",quote=FALSE)
save(cytoMagenta4, file="cytoMagenta4.Rdata")

MagentaModule<-subset(geneInfoALL, moduleColor == "magenta", select = c("geneIDs", "moduleColor"))
write.csv(MagentaModule, "MagentaModule4.csv")

dataMagenta4 <- read.csv("cyto_MagentaEDGE4.csv")

dataMagenta4$fromNode2 = dataMagenta4$fromNode

names(dataMagenta4)[8] <- "genes"

cytoMagenta4 <- dataMagenta4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoMagenta4, "cyto_Magenta4.csv",quote=FALSE)

cytoMagenta44 = cytoMagenta4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoMagenta44, "cyto_MagentaEDGE44.csv",quote=FALSE)


#Purple
PurpleModule4<-subset(geneInfoALL, moduleColor == "purple", select = c("geneIDs", "moduleColor"))
write.csv(PurpleModule4, "PurpleModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[PurpleModule4$geneIDs,]->geneModuleMembership4_Purple#get modulemembership for genes in blue module
geneModuleMembership4_Purple$MEPurple->geneModuleMembershipPurple #correlation of each gene with the with module eigengene for genes in blue module

kFour[,PurpleModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,PurpleModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_Purple)->cytoPurple4


write.csv(cytoPurple4$edgeData, "cyto_PurpleEDGE4.csv",quote=FALSE)
write.csv(cytoPurple4$nodeData, "cyto_PurpleNODE4.csv",quote=FALSE)
save(cytoPurple4, file="cytoPurple4.Rdata")

PurpleModule<-subset(geneInfoALL, moduleColor == "purple", select = c("geneIDs", "moduleColor"))
write.csv(PurpleModule, "PurpleModule4.csv")

dataPurple4 <- read.csv("cyto_PurpleEDGE4.csv")

dataPurple4$fromNode2 = dataPurple4$fromNode

names(dataPurple4)[8] <- "genes"

cytoPurple4 <- dataPurple4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoPurple4, "cyto_Purple4.csv",quote=FALSE)

cytoPurple44 = cytoPurple4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoPurple44, "cyto_PurpleEDGE44.csv",quote=FALSE)

#Red
RedModule4<-subset(geneInfoALL, moduleColor == "red", select = c("geneIDs", "moduleColor"))
write.csv(RedModule4, "RedModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[RedModule4$geneIDs,]->geneModuleMembership4_Red#get modulemembership for genes in blue module
geneModuleMembership4_Red$MERed->geneModuleMembershipRed #correlation of each gene with the with module eigengene for genes in blue module

kFour[,RedModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,RedModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_Red)->cytoRed4


write.csv(cytoRed4$edgeData, "cyto_RedEDGE4.csv",quote=FALSE)
write.csv(cytoRed4$nodeData, "cyto_RedNODE4.csv",quote=FALSE)
save(cytoRed4, file="cytoRed4.Rdata")

RedModule4<-subset(geneInfoALL, moduleColor == "red", select = c("geneIDs", "moduleColor"))
write.csv(RedModule4, "RedModule4.csv")

dataRed4 <- read.csv("cyto_RedEDGE4.csv")

dataRed4$fromNode2 = dataRed4$fromNode

names(dataRed4)[8] <- "genes"

cytoRed4 <- dataRed4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoRed4, "cyto_Red4.csv",quote=FALSE)

cytoRed44 = cytoRed4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoRed44, "cyto_RedEDGE44.csv",quote=FALSE)

#Yellow
YellModule4<-subset(geneInfoALL, moduleColor == "yellow", select = c("geneIDs", "moduleColor"))
write.csv(YellModule4, "YellModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[YellModule4$geneIDs,]->geneModuleMembership4_Yell#get modulemembership for genes in blue module
geneModuleMembership4_Yell$MEYell->geneModuleMembershipYell #correlation of each gene with the with module eigengene for genes in blue module

kFour[,YellModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,YellModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_Yell)->cytoYell4


write.csv(cytoYell4$edgeData, "cyto_YellEDGE4.csv",quote=FALSE)
write.csv(cytoYell4$nodeData, "cyto_YellNODE4.csv",quote=FALSE)
save(cytoYell4, file="cytoYell4.Rdata")

YellModule4<-subset(geneInfoALL, moduleColor == "yellow", select = c("geneIDs", "moduleColor"))
write.csv(YellModule4, "YellModule4.csv")

dataYell4 <- read.csv("cyto_YellEDGE4.csv")

dataYell4$fromNode2 = dataYell4$fromNode

names(dataYell4)[8] <- "genes"

cytoYell4 <- dataYell4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoYell4, "cyto_Yell4.csv",quote=FALSE)

cytoYell44 = cytoYell4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.8)

write.csv(cytoYell44, "cyto_YellEDGE44.csv",quote=FALSE)
