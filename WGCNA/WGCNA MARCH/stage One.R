
if (!require("pacman")) install.packages("pacman")
pacman::p_load(WGCNA, ABAData, httr, readxl, tidyverse, patchwork, ggplotify, janitor, ggrepel, ggtext, glue, dplyr)
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
  geom_hline(yintercept = 0.77, col = "
red" )

p2 = sftOne$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = mean_k, label = powersOne)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  geom_point() +
  geom_text_repel()

pOne + p2

softPowerOne = 8
adjacencyOne = adjacency(kOne, power = softPowerOne);

# Turn adjacency into topological overlap 
TOMOne = TOMsimilarity(adjacencyOne); 

dissTOMOne = 1-TOMOne

# Call the hierarchical clustering function 
geneTreeOne = hclust(as.dist(dissTOMOne), method = "average")
# Plot the resulting clustering tree (dendrogram) 
plot(geneTreeOne, xlab="", sub="", main = glue("Gene clustering on TOM-based dissimilarity\nfor schizophrenia genes during stage One"), labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = 40; 
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
if(ncol(MEsOne) > 2) plotEigengeneNetworks(MEsOne, glue("Eigengene adjacency heatmap for stage  One"), marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

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


allLLIDs<-read_csv("ABALocusLinked.csv")
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
write.csv(GO_per_setOnecont_Blue, file="GO_per_setOne_blackcont.csv")
#Brown
GO_per_setOnecont_Brown <- GO_per_setOnecont[GO_per_setOnecont$Module == "brown",]
write.csv(GO_per_setOnecont_Brown, file="GO_per_setOne_brown.csv")
#Cyan
GO_per_setOnecont_Cyan <- GO_per_setOnecont[GO_per_setOnecont$Module == "cyan",]
write.csv(GO_per_setOnecont_Cyan, file="GO_per_setOne_Cyan.csv")
#GreenYellow
GO_per_setOnecont_GreenYellow <- GO_per_setOnecont[GO_per_setOnecont$Module == "greenyellow",]
write.csv(GO_per_setOnecont_GreenYellow, file="GO_per_setOne_GreenYellow.csv")
#Magenta
GO_per_setOnecont_Magenta <- GO_per_setOnecont[GO_per_setOnecont$Module == "magenta",]
write.csv(GO_per_setOnecont_Magenta, file="GO_per_setOne_Magenta.csv")
#Purple
GO_per_setOnecont_Purple <- GO_per_setOnecont[GO_per_setOnecont$Module == "purple",]
write.csv(GO_per_setOnecont_Purple, file="GO_per_setOne_Purple.csv")
#Red
GO_per_setOnecont_Red <- GO_per_setOnecont[GO_per_setOnecont$Module == "red",]
write.csv(GO_per_setOnecont_Red, file="GO_per_setOne_Red.csv")

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


datablack1 <- read.csv("cyto_blackEDGE1.csv")
datablacl1$fromNode2 = datablack1$fromNode

names(datablack1)[8] <- "genes"

cytoblack1 <- datablack1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoblack11 = cytoblack1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)
write.csv(cytoblack1, "cyto_blackEDGE1.csv",quote=FALSE)
write.csv(cytoblack11, "cyto_blackEDGE11.csv",quote=FALSE)


#Cyan
CyanModule1<-subset(geneInfoALL, moduleColor == "cyan", select = c("geneIDs", "moduleColor"))
write.csv(CyanModule1, "CyanModule1.csv")

geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[CyanModule1$geneIDs,]->geneModuleMembership1_cyan#get modulemembership for genes in blue module
geneModuleMembership1_cyan$MEturquoise->geneModuleMembershipcyan #correlation of each gene with the with module eigengene for genes in blue module

kOne[,cyanModule1$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,cyanModule1$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_cyan)->cytocyan1


write.csv(cytocyan1$edgeData, "cyto_cyanEDGE1.csv",quote=FALSE)
write.csv(cytocyan1$nodeData, "cyto_cyanNODE1.csv",quote=FALSE)
save(cytocyan1, file="cytocyan1.Rdata")

cyanModule<-subset(geneInfoALL, moduleColor == "cyan", select = c("geneIDs", "moduleColor"))
write.csv(cyanModule, "cyanModule.csv")

datacyan1 <- read.csv("cyto_cyanEDGE1.csv")

datacyan1$fromNode2 = datacyan1$fromNode

names(datacyan1)[8] <- "genes"

cytocyan1 <- datacyan1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 

write.csv(cytocyan1, "cyto_cyanEDGE1.csv",quote=FALSE)

cytocyan11 = cytocyan1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytocyan11, "cyto_cyanEDGE11.csv",quote=FALSE)
#Brown
geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

brownModule<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(brownModule, "brownModule.csv")
geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[brownModule$geneIDs,]->geneModuleMembership1_brown#get modulemembership for genes in blue module
geneModuleMembership1_brown$MEturquoise->geneModuleMembershipbrown #correlation of each gene with the with module eigengene for genes in blue module

kOne[,brownModule$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,brownModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_brown)->cytobrown1

#GreenYellow
geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

greenyellowModule<-subset(geneInfoALL, moduleColor == "greenyellow", select = c("geneIDs", "moduleColor"))
write.csv(greenyellowModule, "greenyellowModule.csv")
geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[greenyellowModule$geneIDs,]->geneModuleMembership1_greenyellow#get modulemembership for genes in blue module
geneModuleMembership1_greenyellow$MEgreenyellow->geneModuleMembershipgreenyellow #correlation of each gene with the with module eigengene for genes in blue module

kOne[,greenyellowModule$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,greenyellowModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_greenyellow)->cytogreenyellow1

write.csv(cytogreenyellow1$edgeData, "cyto_greenyellowEDGE1.csv",quote=FALSE)
write.csv(cytogreenyellow1$nodeData, "cyto_greenyellowNODE1.csv",quote=FALSE)
save(cytogreenyellow1, file="cytogreenyellow1.Rdata")

greenyellow1 <- read.csv("cyto_greenyellowEDGE1.csv")

datagreenyellow1$fromNode2 = datagreenyellow1$fromNode

names(datagreenyellow1)[8] <- "genes"

cytogreenyellow1 <- datagreenyellow1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytogreenyellow1, "cyto_greenyellowEDGE1.csv",quote=FALSE)

cytogreenyellow11 = cytogreenyellow1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytogreenyellow11, "cyto_greenyellowEDGE11.csv",quote=FALSE)

#Magenta
geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

MagentaModule<-subset(geneInfoALL, moduleColor == "magenta", select = c("geneIDs", "moduleColor"))
write.csv(MagentaModule, "MagentaModule.csv")
geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[MagentaModule$geneIDs,]->geneModuleMembership1_Magenta#get modulemembership for genes in blue module
geneModuleMembership1_Magenta$MEMagenta->geneModuleMembershipMagenta #correlation of each gene with the with module eigengene for genes in blue module

kOne[,MagentaModule$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,MagentaModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_Magenta)->cytoMagenta1

write.csv(cytoMagenta1$edgeData, "cyto_MagentaEDGE1.csv",quote=FALSE)
write.csv(cytoMagenta1$nodeData, "cyto_MagentaNODE1.csv",quote=FALSE)
save(cytoMagenta1, file="cytoMagenta1.Rdata")

Magenta1 <- read.csv("cyto_MagentaEDGE1.csv")

dataMagenta1$fromNode2 = dataMagenta1$fromNode

names(dataMagenta1)[8] <- "genes"

cytoMagenta1 <- dataMagenta1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoMagenta1, "cyto_MagentaEDGE1.csv",quote=FALSE)

Magenta11 = Magenta1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoMagenta11, "cyto_MagentaEDGE11.csv",quote=FALSE)

#Purple
geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

PurpleModule<-subset(geneInfoALL, moduleColor == "purple", select = c("geneIDs", "moduleColor"))
write.csv(PurpleModule, "PurpleModule.csv")
geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[PurpleModule$geneIDs,]->geneModuleMembership1_Purple#get modulemembership for genes in blue module
geneModuleMembership1_Purple$MEPurple->geneModuleMembershipPurple #correlation of each gene with the with module eigengene for genes in blue module

kOne[,PurpleModule$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,PurpleModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_Purple)->cytoPurple1

write.csv(cytoPurple1$edgeData, "cyto_PurpleEDGE1.csv",quote=FALSE)
write.csv(cytoPurple1$nodeData, "cyto_PurpleNODE1.csv",quote=FALSE)
save(cytoPurple1, file="cytoPurple1.Rdata")

Purple1 <- read.csv("cyto_PurpleEDGE1.csv")

dataPurple1$fromNode2 = dataPurple1$fromNode

names(dataPurple1)[8] <- "genes"

cytoPurple1 <- dataPurple1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoPurple1, "cyto_PurpleEDGE1.csv",quote=FALSE)

Purple11 = Purple1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoPurple11, "cyto_PurpleEDGE11.csv",quote=FALSE)

#Red
geneInfoALL<-data.frame(geneIDs = colnames(kOne), moduleColor=moduleColorsOne)

RedModule<-subset(geneInfoALL, moduleColor == "red", select = c("geneIDs", "moduleColor"))
write.csv(RedModule, "RedModule.csv")
geneModuleMembership1 = as.data.frame(cor(kOne, MEsOne11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[RedModule$geneIDs,]->geneModuleMembership1_Purple#get modulemembership for genes in blue module
geneModuleMembership1_Red$MERed->geneModuleMembershipRed #correlation of each gene with the with module eigengene for genes in blue module

kOne[,RedModule$geneIDs] #get expression values for genes in turquoise module only

cor((kOne[,RedModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_Red)->cytoRed1

write.csv(cytoRed1$edgeData, "cyto_RedEDGE1.csv",quote=FALSE)
write.csv(cytoRed1$nodeData, "cyto_RedNODE1.csv",quote=FALSE)
save(cytoRed1, file="cytoRed1.Rdata")

Purple1 <- read.csv("cyto_PurpleEDGE1.csv")

dataRed1$fromNode2 = dataRed1$fromNode

names(dataRed1)[8] <- "genes"

cytoRed1 <- dataRed1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoRed1, "cyto_RedEDGE1.csv",quote=FALSE)

Red11 = Red1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoRed11, "cyto_RedEDGE11.csv",quote=FALSE)
