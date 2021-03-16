library(WGCNA)



if (!require("pacman")) install.packages("pacman")
pacman::p_load(ABAData, httr, readxl, tidyverse, patchwork, ggplotify, janitor, ggrepel, ggtext, glue, dplyr)

url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx")))

df <- read_excel(temp_file, sheet = 4, skip = 3) 
unlink(temp_file)

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
  geom_hline(yintercept = 0.77, col = "red" )

pThree2 = sftThree$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = mean_k, label = powersThree)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  geom_point() +
  geom_text_repel()

pThree + pThree2

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
if(ncol(MEsThree) > 2) plotEigengeneNetworks(MEsThree, glue("Eigengene adjacency heatmap for stage Three"), marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

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


allLLIDs<-read_csv("ABALocusLinked.csv")
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

#Blue
GO_per_setThreecont_Blue <- GO_per_setThreecont[GO_per_setThreecont$Module == "blue",]
write.csv(GO_per_setThreecont_Blue, file="GO_per_setThree_bluecont.csv")
#Brown
GO_per_setThreecont_Brown <- GO_per_setThreecont[GO_per_setThreecont$Module == "brown",]
write.csv(GO_per_setThreecont_Brown, file="GO_per_setThree_brown.csv")
#Turq
GO_per_setThreecont_Turquoise <- GO_per_setThreecont[GO_per_setThreecont$Module == "turquoise",]
write.csv(GO_per_setThreecont_Turquoise, file="GO_per_setThree_Turquoise.csv")
#Black
GO_per_setThreecont_Black <- GO_per_setThreecont[GO_per_setThreecont$Module == "black",]
write.csv(GO_per_setThreecont_Black, file="GO_per_setThree_Black.csv")
#Green
GO_per_setThreecont_Black <- GO_per_setThreecont[GO_per_setThreecont$Module == "black",]
write.csv(GO_per_setThreecont_Black, file="GO_per_setThree_Black.csv")


pdf("plotEigeneNetworks.pdf")
plotEigengeneNetworks(MEsThree, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(kThree), moduleColor=moduleColorsThree)

blueModule<-subset(geneInfoALL, moduleColor == "blue", select = c("geneIDs", "moduleColor"))
write.csv(blueModule, "blueModuleThree.csv")

geneModuleMembership3 = as.data.frame(cor(kThree, MEsOne33, use = "p")) #module membership for all genes all modules
geneModuleMembership3[blueModule$geneIDs,]->geneModuleMembershipThree_blue #get modulemembership for genes in blue module
geneModuleMembershipThree_blue$MEblue->geneModuleMembership3blue #correlation of each gene with the with module eigengene for genes in blue module

kThree[,blueModule$geneIDs] #get expression values for genes in blue module only

cor((kOne[,blueModule$geneIDs]))->cor3 #correlation of each gene expression value with every other gene in the blue module


rownames(cor3)->vec3 # list of gene names in blue module
exportNetworkToCytoscape(cor3, altNodeNames =vec3, nodeAttr=geneModuleMembershipThree_blue)->cytoblue3
write.csv(cytoblue3$edgeData, "cyto_blueEDGE3.csv",quote=FALSE)
write.csv(cytoblue3$nodeData, "cyto_blueNODE3.csv",quote=FALSE)
save(cytoblue3, file="cytoblue3.Rdata")

pacman::p_load(tidyverse, httr, readxl, janitor,dplyr)


url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
df <- read_excel(temp_file, sheet = 4, skip = 3) # reads into a dataframe. First six rows of the excel file are just header
unlink(temp_file)     # deletes the temporary file

pf <- df %>%
  clean_names() %>%
  separate_rows(gene_s_tagged) %>%
  dplyr::select(genes = gene_s_tagged, p_value) %>%
  distinct() %>% 
  group_by(genes) %>% 
  summarise(p_value = min(p_value, na.rm = T)) %>% 
  ungroup()


datablue3 <- read.csv("cyto_blueEDGE3.csv")

datablue3$fromNode2 = datablue3$fromNode

names(datablue3)[8] <- "genes"

cytoblue3 <- datablue3 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoblue33 = cytoblue3 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)
write.csv(cytoblue3, "cyto_blueEDGE3.csv",quote=FALSE)
write.csv(cytoblue33, "cyto_blueEDGE33.csv",quote=FALSE)


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
  dplyr::filter(weight > 0.7)

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
  dplyr::filter(weight > 0.7)

write.csv(cytobrown33, "cyto_brownEDGE33.csv",quote=FALSE)

#Black
geneInfoALL<-data.frame(geneIDs = colnames(kThree), moduleColor=moduleColorsThree)

BlackModule<-subset(geneInfoALL, moduleColor == "black", select = c("geneIDs", "moduleColor"))
write.csv(BlackModule, "BlackModule3.csv")
geneModuleMembership3 = as.data.frame(cor(kThree, MEsOne33, use = "p")) #module membership for all genes all modules
geneModuleMembership3[BlackModule$geneIDs,]->geneModuleMembership3_Black#get modulemembership for genes in blue module
geneModuleMembership3_Black$MEBlack->geneModuleMembershipBlack #correlation of each gene with the with module eigengene for genes in blue module

kThree[,BlackModule$geneIDs] #get expression values for genes in turquoise module only

cor((kThree[,BlackModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership3_Black)->cytoBlack3


write.csv(cytobBlack3$edgeData, "cyto_BlackEDGE3.csv",quote=FALSE)
write.csv(cytoBlack3$nodeData, "cyto_BlackNODE3.csv",quote=FALSE)
save(cytoBlack3, file="cytoBlack3.Rdata")

dataBlack3 <- read.csv("cyto_BlackEDGE3.csv")

dataBlack3$fromNode2 = dataBlack3$fromNode

names(dataBlack3)[8] <- "genes"

cytoBlack33 <- dataBlack3 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytoBlack33, "cyto_BlackEDGE33.csv",quote=FALSE)

#Green
geneInfoALL<-data.frame(geneIDs = colnames(kThree), moduleColor=moduleColorsThree)

GreenModule<-subset(geneInfoALL, moduleColor == "green", select = c("geneIDs", "moduleColor"))
write.csv(GreenModule, "GreenModule3.csv")
geneModuleMembership3 = as.data.frame(cor(kThree, MEsOne33, use = "p")) #module membership for all genes all modules
geneModuleMembership3[GreenModule$geneIDs,]->geneModuleMembership3_Green#get modulemembership for genes in blue module
geneModuleMembership3_Green$MEGreen->geneModuleMembershipGreen #correlation of each gene with the with module eigengene for genes in blue module

kThree[,GreenModule$geneIDs] #get expression values for genes in turquoise module only

cor((kThree[,GreenModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership3_Green)->cytoGreen3


write.csv(cytobGreen3$edgeData, "cyto_GreenEDGE3.csv",quote=FALSE)
write.csv(cytoGreen3$nodeData, "cyto_GreenNODE3.csv",quote=FALSE)
save(cytoGreen3, file="cytoGreen3.Rdata")

dataGreen3 <- read.csv("cyto_GreenEDGE3.csv")

dataGreen3$fromNode2 = dataGreen3$fromNode

names(dataGreen3)[8] <- "genes"

cytoGreen33 <- dataGreen3 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytoGreen33, "cyto_GreenDGE33.csv",quote=FALSE)

cytobrown33 = cytobrown3 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytobrown33, "cyto_brownEDGE33.csv",quote=FALSE)
