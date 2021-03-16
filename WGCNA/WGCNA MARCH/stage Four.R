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
#Blue

GO_per_setFourcont_Blue <- GO_per_setFourcont[GO_per_setFourcont$Module == "blue",]
write.csv(GO_per_setFourcont_Blue, file="GO_per_setFour_bluecont.csv")
#Brown
GO_per_setFourcont_Brown <- GO_per_setFourcont[GO_per_setFourcont$Module == "brown",]
write.csv(GO_per_setFourcont_Brown, file="GO_per_setFour_brown.csv")
#Yellow
GO_per_setFourcont_Yellow <- GO_per_setFourcont[GO_per_setFourcont$Module == "yellow",]
write.csv(GO_per_setFourcont_Yellow, file="GO_per_setFour_Yellow.csv")


pdf("plotEigeneNetworks.pdf")
plotEigengeneNetworks(MEsFour, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(kFour), moduleColor=moduleColorsFour)

blueModule<-subset(geneInfoALL, moduleColor == "blue", select = c("geneIDs", "moduleColor"))
write.csv(blueModule, "blueModuleFour.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[blueModule$geneIDs,]->geneModuleMembershipFour_blue #get modulemembership for genes in blue module
geneModuleMembershipFour_blue$MEblue->geneModuleMembership4blue #correlation of each gene with the with module eigengene for genes in blue module

kFour[,blueModule$geneIDs] #get expression values for genes in blue module only

cor((kOne[,blueModule$geneIDs]))->cor3 #correlation of each gene expression value with every other gene in the blue module


rownames(cor3)->vec3 # list of gene names in blue module
exportNetworkToCytoscape(cor3, altNodeNames =vec3, nodeAttr=geneModuleMembershipFour_blue)->cytoblue4
write.csv(cytoblue4$edgeData, "cyto_blueEDGE4.csv",quote=FALSE)
write.csv(cytoblue4$nodeData, "cyto_blueNODE4.csv",quote=FALSE)
save(cytoblue4, file="cytoblue4.Rdata")

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


datablue4 <- read.csv("cyto_blueEDGE4.csv")

datablue4$fromNode2 = datablue4$fromNode

names(datablue4)[8] <- "genes"

cytoblue4 <- datablue4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoblue44 = cytoblue4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)
write.csv(cytoblue4, "cyto_blueEDGE4.csv",quote=FALSE)
write.csv(cytoblue44, "cyto_blueEDGE44.csv",quote=FALSE)


#yellow
YellowModule4<-subset(geneInfoALL, moduleColor == "yellow", select = c("geneIDs", "moduleColor"))
write.csv(YellowModule4, "yellowModule4.csv")

geneModuleMembership4 = as.data.frame(cor(kFour, MEsOne44, use = "p")) #module membership for all genes all modules
geneModuleMembership4[YellowModule4$geneIDs,]->geneModuleMembership4_yellow#get modulemembership for genes in blue module
geneModuleMembership4_yellow$MEyellow->geneModuleMembershipyellow #correlation of each gene with the with module eigengene for genes in blue module

kFour[,YellowModule4$geneIDs] #get expression values for genes in turquoise module only

cor((kFour[,YellowModule4$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership4_yellow)->cytoyellow4


write.csv(cytoyellow4$edgeData, "cyto_yellowEDGE4.csv",quote=FALSE)
write.csv(cytoyellow4$nodeData, "cyto_yellowNODE4.csv",quote=FALSE)
save(cytoyellow4, file="cytoyellow4.Rdata")

yellowModule<-subset(geneInfoALL, moduleColor == "yellow", select = c("geneIDs", "moduleColor"))
write.csv(yellowModule, "yellowModule4.csv")

datayellow4 <- read.csv("cyto_yellowEDGE4.csv")

datayellow4$fromNode2 = datayellow4$fromNode

names(datayellow4)[8] <- "genes"

cytoyellow4 <- datayellow4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoyellow4, "cyto_yellow4.csv",quote=FALSE)

cytoyell44 = cytoyellow4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoyell44, "cyto_yellowEDGE44.csv",quote=FALSE)
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


write.csv(cytoyellow4$edgeData, "cyto_yellowEDGE4.csv",quote=FALSE)
write.csv(cytoyellow4$nodeData, "cyto_yellowNODE4.csv",quote=FALSE)
save(cytoyellow4, file="cytoyellow4.Rdata")

yellowModule<-subset(geneInfoALL, moduleColor == "yellow", select = c("geneIDs", "moduleColor"))
write.csv(yellowModule, "yellowModule4.csv")

datayellow4 <- read.csv("cyto_yellowEDGE4.csv")

datayellow4$fromNode2 = datayellow4$fromNode

names(datayellow4)[8] <- "genes"

cytoyellow4 <- datayellow4 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoyellow4, "cyto_yellow4.csv",quote=FALSE)

cytoyell44 = cytoyellow4 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoyell44, "cyto_yellowEDGE44.csv",quote=FALSE)