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



sz_stage_2 <-
  dataset_5_stages %>%
  dplyr::filter(hgnc_symbol %in% all_sz_genes$genes) %>% 
  dplyr::filter(age_category == 2) %>% 
  dplyr::select(-c(entrezgene, ensembl_gene_id, age_category)) %>% 
  pivot_wider(names_from = structure, values_from = signal) %>% 
  column_to_rownames(var = "hgnc_symbol")%>%
  t() %>% 
  scale() %>% 
  t()


kTwo <- sz_stage_2 %>% t()

powersTwo = c(c(1:10), seq(from = 12, to=20, by=2))

sftTwo = pickSoftThreshold(kTwo, powerVector = powersTwo, verbose = 5)

cexOne = 0.9;

pTwo = sftTwo$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = -sign(slope)*sft_r_sq, label = powersTwo)) +
  labs(x = "Soft Threshold (power)",
       y = "Scale Free Topology Model Fit, signed R<sup>2</sup>",
       title = "Scale independence") +
  geom_point() +
  geom_text_repel() +
  theme(axis.title.y = element_textbox_simple(
    orientation = "left-rotated",
    width = NULL)) +
  geom_hline(yintercept = 0.77, col = "red" )

pTwo2 = sftTwo$fitIndices %>% 
  clean_names() %>% 
  ggplot(aes(power, y = mean_k, label = powersTwo)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  geom_point() +
  geom_text_repel()

pTwo + pTwo2

softPowerTwo = 7
adjacencyTwo = adjacency(kTwo, power = softPowerTwo);

# Turn adjacency into topological overlap 
TOMTwo = TOMsimilarity(adjacencyTwo); 

dissTOMTwo = 1-TOMTwo

# Call the hierarchical clustering function 
geneTreeTwo = hclust(as.dist(dissTOMTwo), method = "average")
# Plot the resulting clustering tree (dendrogram) 
plot(geneTreeTwo, xlab="", sub="", main = glue("Gene clustering on TOM-based dissimilarity\nfor schizophrenia genes during stage Two"), labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = 10; 
# Module identification using dynamic tree cut: 
dynamicModsTwo = cutreeDynamic(dendro = geneTreeTwo, 
                               distM = dissTOMTwo, 
                               deepSplit = 2, 
                               pamRespectsDendro = FALSE, 
                               minClusterSize = minModuleSize,
                               cutHeight = 0.95); 

table(dynamicModsTwo)

dynamicColorsTwo = labels2colors(dynamicModsTwo) 
table(dynamicColorsTwo) 

plotDendroAndColors(geneTreeTwo, 
                    dynamicColorsTwo, 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = glue("Age category Two colour dendrogram"))

MEListTwo = moduleEigengenes(sz_stage_2 %>% t(), colors = dynamicModsTwo)


MEsTwo = MEListTwo$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissTwo = 1-cor(MEsTwo); 
# Cluster module eigengenes 
METreeTwo = hclust(as.dist(MEDissTwo), method = "average"); 
# Plot the result 
plot(METreeTwo, main = glue("Clustering of module eigengenes for stage Two"), xlab = "", sub = "")
MEDissThresTwo = 0.9
# Plot the cut line into the dendrogram
abline(h=MEDissThresTwo, col = "red")

mergeTwo = mergeCloseModules(kTwo, dynamicColorsTwo, cutHeight = MEDissThresTwo, verbose = 3)

mergedColorsTwo = mergeTwo$colors;
# Eigengenes of the new merged modules:
mergedMEsTwo = mergeTwo$newMEs

plot(geneTreeTwo, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 0.5, 
     cex.axis = 1.5, 
     cex.main = 2, 
     cex = 0.6)

plotDendroAndColors(geneTreeTwo, 
                    cbind(dynamicColorsTwo, mergedColorsTwo),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



MEListTwo = moduleEigengenes(sz_stage_2 %>% t(), colors = dynamicModsTwo)


#MEDissThresOne = 1
# Plot the cut line into the dendrogram
#abline(h=MEDissThresOne, col = "red")

mergeTwo = mergeCloseModules(kTwo, dynamicColorsTwo, cutHeight = MEDissThresTwo, verbose = 3)

mergedColorsTwo = mergeTwo$colors;
geneInfoALLTwo <- data.frame(geneIDs = colnames(kTwo), 
                             moduleColor=mergedColorsTwo)

mergedMEsTwo = mergeTwo$newMEs
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreeTwo, cbind(dynamicColorsTwo, mergedColorsTwo),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Rename to moduleColors
moduleColorsTwo = mergedColorsTwo
# Construct numerical labels corresponding to the colors
colorOrderTwo = c("grey", standardColors(50));
moduleLabelsTwo = match(moduleColorsTwo, colorOrderTwo)-1;
MEsTwo = mergedMEsTwo;

MEsOne2 = moduleEigengenes(kTwo, moduleColorsTwo)$eigengenes
MEsOne22 = orderMEs(MEsOne2)
#View(MEs50)
MEs2Two = MEListTwo$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDissTwo = 1-cor(MEs2Two); 
# Cluster module eigengenes 
METreeTwo = hclust(as.dist(MEDissTwo), method = "average"); 
# Plot the result 
plot(METreeTwo, main = glue("Clustering of module eigengenes for age category  Two"), xlab = "", sub = "")

dissTOMTwo = 1-TOMsimilarityFromExpr(kTwo, power =softPowerTwo);

plotTOMTwo = dissTOMTwo^softPowerTwo;
# Set diagonal to NA for a nicer plot
diag(plotTOMTwo) = NA;
# Call the plot function
TOMplot(plotTOMTwo, geneTreeTwo, moduleColorsTwo, main = glue("Network heatmap plot, all sz genes, stage Two"))

par(cex = 1.0)
if(ncol(MEsTwo) > 2) plotEigengeneNetworks(MEsTwo, glue("Eigengene dendrogram for stage Two"), marDendro = c(0,4,2,0), plotHeatmaps = FALSE)

par(cex = 1.0)
if(ncol(MEsTwo) > 2) plotEigengeneNetworks(MEsTwo, glue("Eigengene adjacency heatmap for stage  Two"), marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

moduleEigengenes(kTwo,
                 moduleColorsTwo,
                 impute = TRUE,
                 nPC = 1,
                 align = "along average",
                 excludeGrey = FALSE,
                 grey = if (is.numeric(moduleColorsTwo)) 0 else "grey",
                 subHubs =TRUE,
                 softPower =softPowerTwo,
                 scale = TRUE,
                 verbose = 0, indent = 0)

chooseTopHubInEachModule(
  kTwo,
  moduleColorsTwo,
  omitColors = "grey",
  power = softPowerTwo,
  type = "unsigned")

geneModuleMembershipTwo = as.data.frame(cor(kTwo, MEsOne22, use = "p"))

## New

setwd("C:/Users/elena/Documents")

library(anRichment)
library(org.Hs.eg.db)
GOcollection = buildGOcollection(organism = "human")


allLLIDs<-read_csv("ABALocusLinked.csv")
GOenrTwo = enrichmentAnalysis(
  classLabels = moduleColorsTwo, identifiers = allLLIDs$LOCUSLINK_ID,
  refCollection = GOcollection,
  useBackground = "allOrgGenes",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")


GO_per_setTwo <- cbind(GOenrTwo$enrichmentTable$class, GOenrTwo$enrichmentTable$dataSetID,  GOenrTwo$enrichmentTable$dataSetName, GOenrTwo$enrichmentTable$FDR)
GO_per_set_Twocont <- cbind(GOenrTwo$enrichmentTable$class, GOenrTwo$enrichmentTable$dataSetID,  GOenrTwo$enrichmentTable$dataSetName, GOenrTwo$enrichmentTable$FDR,GOenrTwo$enrichmentTable$overlapGenes)
tail(GO_per_setTwo)

colnames(GO_per_setTwo) <- c("Module", "GO term", "GO process", "FDR")
colnames(GO_per_set_Twocont) <- c("Module", "GO Term", "GO Process", "FDR", "Genes")

write.csv(GO_per_set_Twocont, file="GO_per_setTwocont.csv")
# Pull the top (most significant) for each module

top_GO_per_moduleTwo <- GO_per_setTwo[!duplicated(GO_per_setTwo[,1]), ]
colnames(top_GO_per_moduleTwo) <- c("Module", "GO term", "GO process", "FDR")

# Recalculate module eigengenes
MEsTwo = moduleEigengenes(kTwo, moduleColorsTwo)$eigengenes

plotEigengeneNetworks(MEsTwo, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 

write.csv(GO_per_setTwo, file="GO_per_setTwo.csv")

GO_per_setTwo <- as.data.frame(GO_per_setTwo)
GO_per_setTwocont <- as.data.frame(GO_per_set_Twocont)
#Blue

GO_per_setTwocont_Blue <- GO_per_setTwocont[GO_per_setTwocont$Module == "blue",]
write.csv(GO_per_setTwocont_Blue, file="GO_per_setTwo_bluecont.csv")
#Brown
GO_per_setTwocont_Brown <- GO_per_setTwocont[GO_per_setTwocont$Module == "brown",]
write.csv(GO_per_setTwocont_Brown, file="GO_per_setTwo_brown.csv")
#Turquoise
GO_per_setTwocont_Turquoise <- GO_per_setTwocont[GO_per_setTwocont$Module == "turquoise",]
write.csv(GO_per_setTwocont_Turquoise, file="GO_per_setTwo_Turquoise.csv")


pdf("plotEigeneNetworks.pdf")
plotEigengeneNetworks(MEsTwo, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(kTwo), moduleColor=moduleColorsTwo)

blueModule<-subset(geneInfoALL, moduleColor == "blue", select = c("geneIDs", "moduleColor"))
write.csv(blueModule, "blueModuleTwo.csv")

geneModuleMembership2 = as.data.frame(cor(kTwo, MEsOne22, use = "p")) #module membership for all genes all modules
geneModuleMembership2[blueModule$geneIDs,]->geneModuleMembershipTwo_blue #get modulemembership for genes in blue module
geneModuleMembershipTwo_blue$MEblue->geneModuleMembership2blue #correlation of each gene with the with module eigengene for genes in blue module

kTwo[,blueModule$geneIDs] #get expression values for genes in blue module only

cor((kTwo[,blueModule$geneIDs]))->cor3 #correlation of each gene expression value with every other gene in the blue module


rownames(cor3)->vec3 # list of gene names in blue module
exportNetworkToCytoscape(cor3, altNodeNames =vec3, nodeAttr=geneModuleMembershipTwo_blue)->cytoblue2
write.csv(cytoblue2$edgeData, "cyto_blueEDGE2.csv",quote=FALSE)
write.csv(cytoblue2$nodeData, "cyto_blueNODE2.csv",quote=FALSE)
save(cytoblue2, file="cytoblue2.Rdata")

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


datablue2 <- read.csv("cyto_blueEDGE2.csv")

datablue2$fromNode2 = datablue2$fromNode

names(datablue2)[8] <- "genes"

cytoblue2 <- datablue2 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoblue22 = cytoblue2 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)
write.csv(cytoblue2, "cyto_blueEDGE2.csv",quote=FALSE)
write.csv(cytoblue22, "cyto_blueEDGE22.csv",quote=FALSE)


#Turq
turquoiseModule2<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(turquoiseModule2, "turquoiseModule2.csv")

geneModuleMembership2 = as.data.frame(cor(kTwo, MEsOne22, use = "p")) #module membership for all genes all modules
geneModuleMembership2[turquoiseModule1$geneIDs,]->geneModuleMembership2_turquoise#get modulemembership for genes in blue module
geneModuleMembership2_turquoise$MEturquoise->geneModuleMembershipturquoise #correlation of each gene with the with module eigengene for genes in blue module

kTwo[,turquoiseModule2$geneIDs] #get expression values for genes in turquoise module only

cor((kTwo[,turquoiseModule2$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership2_turquoise)->cytoturquoise2


write.csv(cytoturquoise2$edgeData, "cyto_turquoiseEDGE2.csv",quote=FALSE)
write.csv(cytoturquoise2$nodeData, "cyto_turquoiseNODE2.csv",quote=FALSE)
save(cytoturquoise2, file="cytoturquoise2.Rdata")

turquoiseModule<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(turquoiseModule, "turquoiseModule2.csv")

dataturquoise2 <- read.csv("cyto_turquoiseEDGE2.csv")

dataturquoise2$fromNode2 = dataturquoise2$fromNode

names(dataturquoise2)[8] <- "genes"

cytoturquoise2 <- dataturquoise2 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoturquoise2, "cyto_turquoiseEDGE2.csv",quote=FALSE)

cytoturq22 = cytoturquoise2 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoturq22, "cyto_turqEDGE22.csv",quote=FALSE)
#Brown
geneInfoALL<-data.frame(geneIDs = colnames(kTwo), moduleColor=moduleColorsTwo)

brownModule<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(brownModule, "brownModule2.csv")
geneModuleMembership2 = as.data.frame(cor(kTwo, MEsOne22, use = "p")) #module membership for all genes all modules
geneModuleMembership2[brownModule$geneIDs,]->geneModuleMembership2_brown#get modulemembership for genes in blue module
geneModuleMembership2_brown$MEturquoise->geneModuleMembershipbrown #correlation of each gene with the with module eigengene for genes in blue module

kTwo[,brownModule$geneIDs] #get expression values for genes in turquoise module only

cor((kTwo[,brownModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership2_brown)->cytobrown2


write.csv(cytobrown2$edgeData, "cyto_brownEDGE2.csv",quote=FALSE)
write.csv(cytobrown2$nodeData, "cyto_brownNODE2.csv",quote=FALSE)
save(cytobrown2, file="cytobrown2.Rdata")

databrown2 <- read.csv("cyto_brownEDGE2.csv")

databrown2$fromNode2 = databrown2$fromNode

names(databrown2)[8] <- "genes"

cytobrown2 <- databrown2 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytobrown2, "cyto_brownEDGE2.csv",quote=FALSE)

cytobrown22 = cytobrown2 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytobrown22, "cyto_brownEDGE22.csv",quote=FALSE)
