library(WGCNA)



if (!require("pacman")) install.packages("pacman")
pacman::p_load(ABAData, httr, readxl, tidyverse, patchwork, ggplotify, janitor, ggrepel, ggtext, glue, dplyr)
theme_set(theme_minimal())

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
  geom_hline(yintercept = 0.77, col = "red" )

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


allLLIDs<-read_csv("ABALocusLinked.csv")
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
#Black
GO_per_setFourcont_Black <- GO_per_setFourcont[GO_per_setFourcont$Module == "black",]
write.csv(GO_per_setFourcont_Black, file="GO_per_setFour_blackcont.csv")
#Brown
GO_per_setFourcont_Brown <- GO_per_setFourcont[GO_per_setFourcont$Module == "brown",]
write.csv(GO_per_setFourcont_Brown, file="GO_per_setFour_brown.csv")
#GreenYellow
GO_per_setFourcont_GreenYellow <- GO_per_setFourcont[GO_per_setFourcont$Module == "greenyellow",]
write.csv(GO_per_setFourcont_GreenYellow, file="GO_per_setFour_GreenYellow.csv")
#Grey60
GO_per_setFourcont_Grey60 <- GO_per_setFourcont[GO_per_setFourcont$Module == "grey60",]
write.csv(GO_per_setFourcont_Grey60, file="GO_per_setFour_Grey60.csv")
#Magenta
GO_per_setFourcont_Magenta <- GO_per_setFourcont[GO_per_setFourcont$Module == "magenta",]
write.csv(GO_per_setFourcont_Magenta, file="GO_per_setFour_Magenta.csv")
#Purple
GO_per_setFourcont_Purple <- GO_per_setFourcont[GO_per_setFourcont$Module == "purple",]
write.csv(GO_per_setFourcont_Purple, file="GO_per_setFour_Purple.csv")

pdf("plotEigeneNetworks.pdf")
plotEigengeneNetworks(MEsFour, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(kFour), moduleColor=moduleColorsFour)

BlackeModule<-subset(geneInfoALL, moduleColor == "black", select = c("geneIDs", "moduleColor"))
write.csv(BlackModule, "BlackModuleFour.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[BlackModule$geneIDs,]->geneModuleMembershipFour_Black #get modulemembership for genes in blue module
geneModuleMembershipFour_Black$MEblue->geneModuleMembership4Black #correlation of each gene with the with module eigengene for genes in blue module

kFour[,BlackModule$geneIDs] #get expression values for genes in blue module only

cor((kOne[,BlackModule$geneIDs]))->cor3 #correlation of each gene expression value with every other gene in the blue module


rownames(cor3)->vec3 # list of gene names in blue module
exportNetworkToCytoscape(cor3, altNodeNames =vec3, nodeAttr=geneModuleMembershipFour_Black)->cytoBlack4
write.csv(cytoBlack4$edgeData, "cyto_BlackEDGE4.csv",quote=FALSE)
write.csv(cytoBlack4$nodeData, "cyto_BlackNODE4.csv",quote=FALSE)
save(cytoBlack4, file="cytoBlack4.Rdata")

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


dataBlack4 <- read.csv("cyto_BlackEDGE4.csv")

dataBlack4$fromNode2 = dataBlack4$fromNode

names(dataBlack4)[8] <- "genes"

cytoBlack4 <- databBlack4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoBlack44 = cytoBlack4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)
write.csv(cytoBlack4, "cyto_BlackEDGE4.csv",quote=FALSE)
write.csv(cytoBlack44, "cyto_BlackEDGE44.csv",quote=FALSE)


#Greenyellow
GreenyellowModule4<-subset(geneInfoALL, moduleColor == "greenyellow", select = c("geneIDs", "moduleColor"))
write.csv(GreenyellowModule4, "GreenyellowModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[GreenyellowModule4$geneIDs,]->geneModuleMembership4_Greenyellow#get modulemembership for genes in blue module
geneModuleMembership4_Greenyellow$MEGreenyellow->geneModuleMembershipGreenyellow #correlation of each gene with the with module eigengene for genes in blue module

kFour[,GreenyellowModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,GreenyellowModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module

rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_Greenyellow)->cytoGreenyellow4

write.csv(cytoGreenyellow4$edgeData, "cyto_GreenyellowEDGE4.csv",quote=FALSE)
write.csv(cytoGreenyellow4$nodeData, "cyto_GreenyellowNODE4.csv",quote=FALSE)
save(cytoGreenyellow4, file="cytoGreenyellow4.Rdata")

GreenyellowModule<-subset(geneInfoALL, moduleColor == "greenyellow", select = c("geneIDs", "moduleColor"))
write.csv(GreenyellowModule, "GreenyellowModule4.csv")

dataGreenyellow4 <- read.csv("cyto_GreenyellowEDGE4.csv")

dataGreenyellow4$fromNode2 = dataGreenyellow4$fromNode

names(dataGreenyellow4)[8] <- "genes"

cytoGreenyellow4 <- dataGreenyellow4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 
write.csv(cytoGreenyellow4, "cyto_Greenyellow4.csv",quote=FALSE)
cytoGreenyellow44 = cytoGreenyellow4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)
write.csv(cytoGreenyellow44, "cyto_GreenyellowEDGE44.csv",quote=FALSE)

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
  dplyr::filter(weight > 0.7)

write.csv(cytobrown44, "cyto_brownEDGE44.csv",quote=FALSE)

#Grey60

Grey60Module4<-subset(geneInfoALL, moduleColor == "grey", select = c("geneIDs", "moduleColor"))
write.csv(Grey60Module4, "Grey60Module4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[Grey60Module4$geneIDs,]->geneModuleMembership4_Grey60#get modulemembership for genes in blue module
geneModuleMembership4_Grey60$MEGrey60->geneModuleMembershipGrey60 #correlation of each gene with the with module eigengene for genes in blue module

kFour[,Grey60Module4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,Grey60Module4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_Grey60)->cytoGrey604


write.csv(cytoGrey604$edgeData, "cyto_Grey60EDGE4.csv",quote=FALSE)
write.csv(cytoGrey604$nodeData, "cyto_Grey60NODE4.csv",quote=FALSE)
save(cytoGrey604, file="cytoGrey604.Rdata")

Grey60Module<-subset(geneInfoALL, moduleColor == "grey60", select = c("geneIDs", "moduleColor"))
write.csv(Grey60Module, "Grey60Module4.csv")

dataGrey604 <- read.csv("cyto_Grey60EDGE4.csv")

dataGrey604$fromNode2 = dataGrey604$fromNode

names(dataGrey604)[8] <- "genes"

cytoGrey604 <- dataGrey604 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoGrey604, "cyto_Grey604.csv",quote=FALSE)

cytoGrey6044 = cytoGrey604 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoGrey6044, "cyto_Grey60EDGE44.csv",quote=FALSE)

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
  dplyr::filter(weight > 0.7)

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
  dplyr::filter(weight > 0.7)

write.csv(cytoPurple44, "cyto_PurpleEDGE44.csv",quote=FALSE)
