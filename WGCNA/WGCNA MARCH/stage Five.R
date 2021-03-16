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



sz_stage_5 <-
  dataset_5_stages %>%
  dplyr::filter(hgnc_symbol %in% all_sz_genes$genes) %>% 
  dplyr::filter(age_category == 5) %>% 
  dplyr::select(-c(entrezgene, ensembl_gene_id, age_category)) %>% 
  pivot_wider(names_from = structure, values_from = signal) %>% 
  column_to_rownames(var = "hgnc_symbol")%>%
  t() %>% 
  scale() %>% 
  t()


kFive <- sz_stage_5 %>% t()

powersFive = c(c(1:10), seq(from = 12, to=20, by=2))

sftFive = pickSoftThreshold(kFive, powerVector = powersFive, verbose = 5)

cexOne = 0.9;

pFive = sftFive$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = -sign(slope)*sft_r_sq, label = powersFive)) +
  labs(x = "Soft Threshold (power)",
       y = "Scale Free Topology Model Fit, signed R<sup>2</sup>",
       title = "Scale independence") +
  geom_point() +
  geom_text_repel() +
  theme(axis.title.y = element_textbox_simple(
    orientation = "left-rotated",
    width = NULL)) +
  geom_hline(yintercept = 0.77, col = "red" )

pFive2 = sftFive$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = mean_k, label = powersFive)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  geom_point() +
  geom_text_repel()

pFive + pFive2

softPowerFive = 10
adjacencyFive = adjacency(kFive, power = softPowerFive);

# Turn adjacency into topological overlap 
TOMFive = TOMsimilarity(adjacencyFive); 

dissTOMFive = 1-TOMFive

# Call the hierarchical clustering function 
geneTreeFive = hclust(as.dist(dissTOMFive), method = "average")
# Plot the resulting clustering tree (dendrogram) 
plot(geneTreeFive, xlab="", sub="", main = glue("Gene clustering on TOM-based dissimilarity\nfor schizophrenia genes during stage Five"), labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = 10; 
# Module identification using dynamic tree cut: 
dynamicModsFive = cutreeDynamic(dendro = geneTreeFive, 
                                distM = dissTOMFive, 
                                deepSplit = 2, 
                                pamRespectsDendro = FALSE, 
                                minClusterSize = minModuleSize,
                                cutHeight = 0.95); 

table(dynamicModsFive)

dynamicColorsFive = labels2colors(dynamicModsFive) 
table(dynamicColorsFive) 

plotDendroAndColors(geneTreeFive, 
                    dynamicColorsFive, 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = glue("Age category Five colour dendrogram"))

MEListFive = moduleEigengenes(sz_stage_5 %>% t(), colors = dynamicModsFive)


MEsFive = MEListFive$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissFive = 1-cor(MEsFive); 
# Cluster module eigengenes 
METreeFive = hclust(as.dist(MEDissFive), method = "average"); 
# Plot the result 
plot(METreeFive, main = glue("Clustering of module eigengenes for stage Five"), xlab = "", sub = "")
MEDissThresFive = 0.3
# Plot the cut line into the dendrogram
abline(h=MEDissThresFive, col = "red")

mergeFive = mergeCloseModules(kFive, dynamicColorsFive, cutHeight = MEDissThresFive, verbose = 3)

mergedColorsFive = mergeFive$colors;
# Eigengenes of the new merged modules:
mergedMEsFive = mergeFive$newMEs

plot(geneTreeFive, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 0.5, 
     cex.axis = 1.5, 
     cex.main = 2, 
     cex = 0.6)

plotDendroAndColors(geneTreeFive, 
                    cbind(dynamicColorsFive, mergedColorsFive),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



MEListFive = moduleEigengenes(sz_stage_5 %>% t(), colors = dynamicModsFive)


#MEDissThresOne = 1
# Plot the cut line into the dendrogram
#abline(h=MEDissThresOne, col = "red")

mergeFive = mergeCloseModules(kFive, dynamicColorsFive, cutHeight = MEDissThresFive, verbose = 3)

mergedColorsFive = mergeFive$colors;
geneInfoALLFive <- data.frame(geneIDs = colnames(kFive), 
                              moduleColor=mergedColorsFive)

mergedMEsFive = mergeFive$newMEs
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreeFive, cbind(dynamicColorsFive, mergedColorsFive),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Rename to moduleColors
moduleColorsFive = mergedColorsFive
# Construct numerical labels corresponding to the colors
colorOrderFive = c("grey", standardColors(50));
moduleLabelsFive = match(moduleColorsFive, colorOrderFive)-1;
MEsFive = mergedMEsFive;

MEsOne5 = moduleEigengenes(kFive, moduleColorsFive)$eigengenes
MEsOne55 = orderMEs(MEsOne5)
#View(MEs50)
MEs4Five = MEListFive$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissFive = 1-cor(MEs4Five); 
# Cluster module eigengenes 
METreeFive = hclust(as.dist(MEDissFive), method = "average"); 
# Plot the result 
plot(METreeFive, main = glue("Clustering of module eigengenes for age category  Five"), xlab = "", sub = "")

dissTOMFive = 1-TOMsimilarityFromExpr(kFive, power =softPowerFive);

plotTOMFive = dissTOMFive^softPowerFive;
# Set diagonal to NA for a nicer plot
diag(plotTOMFive) = NA;
# Call the plot function
TOMplot(plotTOMFFive, geneTreeFive, moduleColorsFive, main = glue("Network heatmap plot, all sz genes, stage Five"))

par(cex = 1.0)
if(ncol(MEsFive) > 2) plotEigengeneNetworks(MEsFive, glue("Eigengene dendrogram for stage Five"), marDendro = c(0,4,2,0), plotHeatmaps = FALSE)

par(cex = 1.0)
if(ncol(MEsFive) > 2) plotEigengeneNetworks(MEsFive, glue("Eigengene adjacency heatmap for stage  Five"), marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

moduleEigengenes(kFive,
                 moduleColorsFive,
                 impute = TRUE,
                 nPC = 1,
                 align = "along average",
                 excludeGrey = FALSE,
                 grey = if (is.numeric(moduleColorsTwo)) 0 else "grey",
                 subHubs =TRUE,
                 softPower =softPowerFive,
                 scale = TRUE,
                 verbose = 0, indent = 0)

chooseTopHubInEachModule(
  kFive,
  moduleColorsFive,
  omitColors = "grey",
  power = softPowerFour,
  type = "unsigned")

geneModuleMembershipFive = as.data.frame(cor(kFive, MEsOne44, use = "p"))

## New

setwd("C:/Users/elena/Documents")

library(anRichment)
library(org.Hs.eg.db)
GOcollection = buildGOcollection(organism = "human")


allLLIDs<-read_csv("ABALocusLinked.csv")
GOenrFive = enrichmentAnalysis(
  classLabels = moduleColorsFive, identifiers = allLLIDs$LOCUSLINK_ID,
  refCollection = GOcollection,
  useBackground = "allOrgGenes",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")


GO_per_setFive <- cbind(GOenrFive$enrichmentTable$class, GOenrFive$enrichmentTable$dataSetID,  GOenrFive$enrichmentTable$dataSetName, GOenrFive$enrichmentTable$FDR)
GO_per_set_Fivecont <- cbind(GOenrFive$enrichmentTable$class, GOenrFive$enrichmentTable$dataSetID,  GOenrFive$enrichmentTable$dataSetName, GOenrFive$enrichmentTable$FDR,GOenrFive$enrichmentTable$overlapGenes)
tail(GO_per_setFive)

colnames(GO_per_setFive) <- c("Module", "GO term", "GO process", "FDR")
colnames(GO_per_set_Fivecont) <- c("Module", "GO Term", "GO Process", "FDR", "Genes")

write.csv(GO_per_set_Fivecont, file="GO_per_setFivecont.csv")
# Pull the top (most significant) for each module

top_GO_per_moduleFive <- GO_per_setFive[!duplicated(GO_per_setFive[,1]), ]
colnames(top_GO_per_moduleFive) <- c("Module", "GO term", "GO process", "FDR")

# Recalculate module eigengenes
MEsFive = moduleEigengenes(kFive, moduleColorsFive)$eigengenes

plotEigengeneNetworks(MEsFive, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 

write.csv(GO_per_setFive, file="GO_per_setFive.csv")

GO_per_setFive <- as.data.frame(GO_per_setFive)
GO_per_setFivecont <- as.data.frame(GO_per_set_Fivecont)
#Black
GO_per_setFivecont_Black <- GO_per_setFivecont[GO_per_setFivecont$Module == "black",]
write.csv(GO_per_setFivecont_Black, file="GO_per_setFive_blackcont.csv")
#Blue
GO_per_setFivecont_Blue <- GO_per_setFivecont[GO_per_setFivecont$Module == "blue",]
write.csv(GO_per_setFivecont_Blue, file="GO_per_setFive_Blue.csv")
#Greenyellow
GO_per_setFivecont_Greenyellow <- GO_per_setFivecont[GO_per_setFivecont$Module == "greenyellow",]
write.csv(GO_per_setFivecont_Greenyellow, file="GO_per_setFive_Greenyellow.csv")
#Brown
GO_per_setFivecont_Brown <- GO_per_setFivecont[GO_per_setFivecont$Module == "brown",]
write.csv(GO_per_setFivecont_Brown, file="GO_per_setFive_Brown.csv")
#Lightcyan
GO_per_setFivecont_Lightcyan <- GO_per_setFivecont[GO_per_setFivecont$Module == "lightcyan",]
write.csv(GO_per_setFivecont_Lightcyan, file="GO_per_setFive_Lightcyan.csv")
#Magenta
GO_per_setFivecont_Magenta <- GO_per_setFivecont[GO_per_setFivecont$Module == "magenta",]
write.csv(GO_per_setFivecont_Magenta, file="GO_per_setFive_Magenta.csv")
#Pink
GO_per_setFivecont_Pink <- GO_per_setFivecont[GO_per_setFivecont$Module == "pink",]
write.csv(GO_per_setFivecont_Pink, file="GO_per_setFive_Pink.csv")
#Purple
GO_per_setFivecont_Purple <- GO_per_setFivecont[GO_per_setFivecont$Module == "purple",]
write.csv(GO_per_setFivecont_Purple, file="GO_per_setFive_Purple.csv")
#Red
GO_per_setFivecont_Red <- GO_per_setFivecont[GO_per_setFivecont$Module == "red",]
write.csv(GO_per_setFivecont_Red, file="GO_per_setFive_Red.csv")


pdf("plotEigeneNetworks.pdf")
plotEigengeneNetworks(MEsFive, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabel = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(kFive), moduleColor=moduleColorsFive)

blueModule<-subset(geneInfoALL, moduleColor == "blue", select = c("geneIDs", "moduleColor"))
write.csv(blueModule, "blueModuleFive.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[blueModule$geneIDs,]->geneModuleMembershipFive_blue #get modulemembership for genes in blue module
geneModuleMembershipFive_blue$MEblue->geneModuleMembership45lue #correlation of each gene with the with module eigengene for genes in blue module

kFive[,blueModule$geneIDs] #get expression values for genes in blue module only

cor((kFive[,blueModule$geneIDs]))->cor3 #correlation of each gene expression value with every other gene in the blue module


rownames(cor3)->vec3 # list of gene names in blue module
exportNetworkToCytoscape(cor3, altNodeNames =vec3, nodeAttr=geneModuleMembershipFive_blue)->cytoblue5
write.csv(cytoblue5$edgeData, "cyto_blueEDGE5.csv",quote=FALSE)
write.csv(cytoblue5$nodeData, "cyto_blueNODE5.csv",quote=FALSE)
save(cytoblue5, file="cytoblue5.Rdata")

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


datablue5 <- read.csv("cyto_blueEDGE5.csv")

datablue5$fromNode2 = datablue5$fromNode

names(datablue5)[8] <- "genes"

cytoblue5 <- datablue5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoblue55 = cytoblue5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)
write.csv(cytoblue5, "cyto_blueEDGE5.csv",quote=FALSE)
write.csv(cytoblue55, "cyto_blueEDGE55.csv",quote=FALSE)


#Greenyellow
GreenyellowModule5<-subset(geneInfoALL, moduleColor == "greenyellow", select = c("geneIDs", "moduleColor"))
write.csv(GreenyellowModule5, "GreenyellowModule5.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[GreenyellowModule5$geneIDs,]->geneModuleMembership5_Greenyellow#get modulemembership for genes in blue module
geneModuleMembership5_Greenyellow$MEGreenyellow->geneModuleMembershipGreenyellow#correlation of each gene with the with module eigengene for genes in blue module

kFive[,GreenyellowModule5$geneIDs] #get expression values for genes in turquoise module only

cor((kFive[,GreenyellowModule5$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership5_Greenyellow)->cytoGreenyellow5


write.csv(cytoGreenyellow5$edgeData, "cyto_GreenyellowEDGE5.csv",quote=FALSE)
write.csv(cytoGreenyellow5$nodeData, "cyto_GreenyellowNODE5.csv",quote=FALSE)
save(cytoGreenyellow5, file="cytoGreenyellow5.Rdata")

GreenyellowModule<-subset(geneInfoALL, moduleColor == "greenyellow", select = c("geneIDs", "moduleColor"))
write.csv(GreenyellowModule, "GreenyellowModule5.csv")

dataGreenyellow5 <- read.csv("cyto_GreenyellowEDGE5.csv")

dataGreenyellow5$fromNode2 = dataGreenyellow5$fromNode

names(dataGreenyellow5)[8] <- "genes"

cytoGreenyellow5 <- dataGreenyellow5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoGreenyellow5, "cyto_Greenyellow5.csv",quote=FALSE)

cytoGreenyellow55 = cytoGreenyellow5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoGreenyellow55, "cyto_GreenyellowEDGE55.csv",quote=FALSE)
#Brown

BrownModule5<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(BrownModule4, "BrownModule5.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[BrownModule5$geneIDs,]->geneModuleMembership5_Brown#get modulemembership for genes in blue module
geneModuleMembership5_Brown$MEBrown->geneModuleMembershipBrown #correlation of each gene with the with module eigengene for genes in blue module

kFive[,BrownModule5$geneIDs] #get expression values for genes in turquoise module only


cor((kFive[,BrownModule5$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership5_Brown)->cytoBrown5


write.csv(cytoBrown5$edgeData, "cyto_BrownEDGE5.csv",quote=FALSE)
write.csv(cytoBrown5$nodeData, "cyto_BrownNODE5.csv",quote=FALSE)
save(cytoBrown5, file="cytoBrown5.Rdata")

BrownModule<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(BrownModule, "BrownModule5.csv")

dataBrown5 <- read.csv("cyto_BrownEDGE5.csv")

dataBrown5$fromNode2 = dataBrown5$fromNode

names(dataBrown5)[8] <- "genes"

cytoBrown5 <- dataBrown5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoBrown5, "cyto_Brown5.csv",quote=FALSE)

cytoBrown55 = cytoBrown5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoBrown55, "cyto_BrownDGE55.csv",quote=FALSE)

#Black

BlackModule5<-subset(geneInfoALL, moduleColor == "black", select = c("geneIDs", "moduleColor"))
write.csv(BlackModule4, "BlackModule5.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[BlackModule5$geneIDs,]->geneModuleMembership5_Black#get modulemembership for genes in blue module
geneModuleMembership5_Black$MEBlack->geneModuleMembershipBlack #correlation of each gene with the with module eigengene for genes in blue module

kFive[,BlackModule5$geneIDs] #get expression values for genes in turquoise module only


cor((kFive[,BlackModule5$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership5_Black)->cytoBlack5


write.csv(cytoBlack5$edgeData, "cyto_BlackEDGE5.csv",quote=FALSE)
write.csv(cytoBlack5$nodeData, "cyto_BlackNODE5.csv",quote=FALSE)
save(cytoBlack5, file="cytoBlack5.Rdata")

BlackModule<-subset(geneInfoALL, moduleColor == "black", select = c("geneIDs", "moduleColor"))
write.csv(BlackModule, "BlackModule5.csv")

dataBlack5 <- read.csv("cyto_BlackEDGE5.csv")

dataBlack5$fromNode2 = dataBlack5$fromNode

names(dataBlack5)[8] <- "genes"

cytoBlack5 <- dataBlack5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoBlack5, "cyto_Black5.csv",quote=FALSE)

cytoBlack55 = cytoBlack5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoBlack55, "cyto_BlackEDGE55.csv",quote=FALSE)

#Lightcyan

LightcyanModule5<-subset(geneInfoALL, moduleColor == "lightcyan", select = c("geneIDs", "moduleColor"))
write.csv(LightcyanModule4, "LightcyanModule5.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[LightcyanModule5$geneIDs,]->geneModuleMembership5_Lightcyan#get modulemembership for genes in blue module
geneModuleMembership5_Lightcyan$MELightcyan->geneModuleMembershipLightcyan #correlation of each gene with the with module eigengene for genes in blue module

kFive[,LightcyanModule5$geneIDs] #get expression values for genes in turquoise module only


cor((kFive[,LightcyanModule5$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership5_Lightcyan)->cytoLightcyan5


write.csv(cytoLightcyan5$edgeData, "cyto_LightcyanEDGE5.csv",quote=FALSE)
write.csv(cytoLightcyan5$nodeData, "cyto_LightcyanNODE5.csv",quote=FALSE)
save(cytoLightcyan5, file="cytoLightcyan5.Rdata")

LightcyanModule<-subset(geneInfoALL, moduleColor == "lightcyan", select = c("geneIDs", "moduleColor"))
write.csv(LightcyanModule, "LightcyanModule5.csv")

dataLightcyan5 <- read.csv("cyto_LightcyanEDGE5.csv")

dataLightcyan5$fromNode2 = dataLightcyan5$fromNode

names(dataLightcyan5)[8] <- "genes"

cytoLightcyan5 <- dataLightcyan5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoLightcyan5, "cyto_Lightcyan5.csv",quote=FALSE)

cytoLightcyan55 = cytoLightcyan5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoLightcyan55, "cyto_LightcyanEDGE55.csv",quote=FALSE)

#Magenta

MagentaModule5<-subset(geneInfoALL, moduleColor == "magenta", select = c("geneIDs", "moduleColor"))
write.csv(MagentaModule4, "MagentaModule5.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[MagentaModule5$geneIDs,]->geneModuleMembership5_Magenta#get modulemembership for genes in blue module
geneModuleMembership5_Magenta$MEMagenta->geneModuleMembershipMagenta #correlation of each gene with the with module eigengene for genes in blue module

kFive[,MagentaModule5$geneIDs] #get expression values for genes in turquoise module only

cor((kFive[,MagentaModule5$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module

rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership5_Magenta)->cytoMagenta5

write.csv(cytoMagenta5$edgeData, "cyto_MagentaEDGE5.csv",quote=FALSE)
write.csv(cytoMagenta5$nodeData, "cyto_MagentaNODE5.csv",quote=FALSE)
save(cytoMagenta5, file="cytoMagenta5.Rdata")

dataMagenta5 <- read.csv("cyto_MagentaEDGE5.csv")

dataMagenta5$fromNode2 = dataMagenta5$fromNode

names(dataMagenta5)[8] <- "genes"

cytoMagenta5 <- dataMagenta5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytoMagenta5, "cyto_Magenta5.csv",quote=FALSE)

cytoMagenta55 = cytoMagenta5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoMagenta55, "cyto_MagentaEDGE55.csv",quote=FALSE)

#Pink

PinkModule5<-subset(geneInfoALL, moduleColor == "pink", select = c("geneIDs", "moduleColor"))
write.csv(PinkModule4, "PinkModule5.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[PinkModule5$geneIDs,]->geneModuleMembership5_Pink#get modulemembership for genes in blue module
geneModuleMembership5_Pink$MEPink->geneModuleMembershipPink #correlation of each gene with the with module eigengene for genes in blue module

kFive[,PinkModule5$geneIDs] #get expression values for genes in turquoise module only

cor((kFive[,PinkModule5$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module

rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership5_Pink)->cytoPink5

write.csv(cytoPink5$edgeData, "cyto_PinkEDGE5.csv",quote=FALSE)
write.csv(cytoPink5$nodeData, "cyto_PinkNODE5.csv",quote=FALSE)
save(cytoPink5, file="cytoPink5.Rdata")

dataPink5 <- read.csv("cyto_PinkEDGE5.csv")

dataPink5$fromNode2 = dataPink5$fromNode

names(dataPink5)[8] <- "genes"

cytoPink5 <- dataPink5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytoPink5, "cyto_Pink5.csv",quote=FALSE)

cytoPink55 = cytoPink5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoPink55, "cyto_PinkEDGE55.csv",quote=FALSE)

#Purple
PurpleModule5<-subset(geneInfoALL, moduleColor == "purple", select = c("geneIDs", "moduleColor"))
write.csv(PurpleModule4, "PurpleModule5.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[PurpleModule5$geneIDs,]->geneModuleMembership5_Purple#get modulemembership for genes in blue module
geneModuleMembership5_Purple$MEPurple->geneModuleMembershipPurple #correlation of each gene with the with module eigengene for genes in blue module

kFive[,PurpleModule5$geneIDs] #get expression values for genes in turquoise module only

cor((kFive[,PurpleModule5$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module

rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership5_Purple)->cytoPurple5

write.csv(cytoPurple5$edgeData, "cyto_PurpleEDGE5.csv",quote=FALSE)
write.csv(cytoPurple5$nodeData, "cyto_PurpleNODE5.csv",quote=FALSE)
save(cytoPurple5, file="cytoPurple5.Rdata")

dataPurple5 <- read.csv("cyto_PurpleEDGE5.csv")

dataPurple5$fromNode2 = dataPurple5$fromNode

names(dataPurple5)[8] <- "genes"

cytoPurple5 <- dataPurple5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytoPurple5, "cyto_Purple5.csv",quote=FALSE)

cytoPurple55 = cytoPurplek5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoPurple55, "cyto_PurpleEDGE55.csv",quote=FALSE)

#Red
RedModule5<-subset(geneInfoALL, moduleColor == "red", select = c("geneIDs", "moduleColor"))
write.csv(RedModule4, "RedModule5.csv")

geneModuleMembership5 = as.data.frame(cor(kFive, MEsOne55, use = "p")) #module membership for all genes all modules
geneModuleMembership5[RedModule5$geneIDs,]->geneModuleMembership5_Red#get modulemembership for genes in blue module
geneModuleMembership5_Red$MEPurple->geneModuleMembershipRed #correlation of each gene with the with module eigengene for genes in blue module

kFive[,RedModule5$geneIDs] #get expression values for genes in turquoise module only

cor((kFive[,RedModule5$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module

rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership5_Red)->cytoRed5

write.csv(cytoRed5$edgeData, "cyto_RedEDGE5.csv",quote=FALSE)
write.csv(cytoRed5$nodeData, "cyto_RedNODE5.csv",quote=FALSE)
save(cytoRed5, file="cytoRed5.Rdata")

dataRed5 <- read.csv("cyto_RedEDGE5.csv")

dataRed5$fromNode2 = dataRed5$fromNode

names(dataRed5)[8] <- "genes"

cytoRed5 <- dataRed5 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytoRed5, "cyto_Red5.csv",quote=FALSE)

cytoRed55 = cytoRed5 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoRed55, "cyto_RedEDGE55.csv",quote=FALSE)
