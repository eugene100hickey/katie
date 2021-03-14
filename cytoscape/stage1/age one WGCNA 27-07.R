

# Tue May 12 11:15:45 2020 ------------------------------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(ABAData, httr, readxl, tidyverse, WGCNA, igraph, pins) 
library(anRichment)
library(org.Hs.eg.db)
GOcollection = buildGOcollection(organism = "human")


allLLIDs<-read_csv("locuslinkedIDs.csv")


url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
df <- read_excel(temp_file, sheet = 4, skip = 3) # reads into a dataframe. First six rows of the excel file are just header
unlink(temp_file)     # deletes the temporary file
pin(x = df, name = "pardinas_genes")


 ##################################################################
# makes sz_genes, a dataframe with a single column of the CLOZUK genes
#######################################################################

all_sz_genes <- df$`Gene(s) tagged` %>% 
  str_split(",") %>% 
  unlist() %>% 
  str_trim() %>% 
  as.data.frame() %>% 
  distinct()
names(all_sz_genes) <- "genes"
#View(all_sz_genes)


# Tue May 12 11:15:54 2020 ------------------------------


data(dataset_5_stages)

sz_stage_1 <- dataset_5_stages %>%
  filter(hgnc_symbol %in% all_sz_genes$genes) %>% 
  filter(age_category == 1) %>% 
  dplyr::select(-c(entrezgene, ensembl_gene_id, age_category)) %>% 
  pivot_wider(names_from = structure, values_from = signal) %>% 
  column_to_rownames(var = "hgnc_symbol")%>%
  t() %>% 
  scale() %>% 
  t()

sz_stage_1 %>% t()

# Tue May 12 11:25:55 2020 ------------------------------
#sz_stage_4<-k4 <- t(k4)
k1 <- sz_stage_1 %>% t()
# Tue May 12 11:26:03 2020 ------------------------------


# Part 2:  Step-by-step network construction and module detection

# Choose a set of soft-thresholding powers
powers1 = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft1 = pickSoftThreshold(k1, powerVector = powers1, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft1$fitIndices[,1], -sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft1$fitIndices[,1], -sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],
     labels=powers1,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.55,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft1$fitIndices[,1], sft1$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft1$fitIndices[,1], sft1$fitIndices[,5], labels=powers1, cex=cex1,col="red")

#Calculating adjacincies with soft threshold of 7
softPower = 6
adjacency = adjacency(k1, power = softPower);

# Turn adjacency into topological overlap 
TOM1 = TOMsimilarity(adjacency); 
dissTOM1 = 1-TOM1

# Call the hierarchical cluste ring function 
geneTree1 = hclust(as.dist(dissTOM1), method = "average"); 
# Plot the resulting clustering tree (dendrogram) 
sizeGrWindow(12,9) 
plot(geneTree1, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = 10; 
# Module identification using dynamic tree cut: 
dynamicMods1 = cutreeDynamic(dendro = geneTree1, distM = dissTOM1, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize); 
table(dynamicMods1)

# Convert numeric lables into colors 
dynamicColors1 = labels2colors(dynamicMods1) 
table(dynamicColors1) 
# Plot the dendrogram and colors underneath 
sizeGrWindow(8,6) 
q1<-plotDendroAndColors(geneTree1, dynamicColors1, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Age one gene dendrogram")

# Calculate eigengenes 

# Tue May 12 11:35:46 2020 ------------------------------
MEList1 = moduleEigengenes(sz_stage_1 %>% t(), colors = dynamicMods1)
# Tue May 12 11:35:52 2020 ------------------------------

MEs = MEList1$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDiss1 = 1-cor(MEs); 
# Cluster module eigengenes 
METree1 = hclust(as.dist(MEDiss1), method = "average"); 
# Plot the result 
sizeGrWindow(7, 6) 
plot(METree1, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres1 = 1
# Plot the cut line into the dendrogram
abline(h=MEDissThres1, col = "red")
# Call an automatic merging function
merge1 = mergeCloseModules(k1, dynamicColors1, cutHeight = MEDissThres1, verbose = 3) # missing the 1 here
# The merged module colors
mergedColors1 = merge1$colors;
# Eigengenes of the new merged modules:
mergedMEs1 = merge1$newMEs
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree1, cbind(dynamicColors1, mergedColors1), #mergedColors1 was called mergedColours so I have corrected
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors1 = mergedColors1
# Construct numerical labels corresponding to the colors
colorOrder1 = c("grey", standardColors(50));
moduleLabels1 = match(moduleColors1, colorOrder1)-1;
MEs1 = mergedMEs1;
# Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree1, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")
MEs10 = moduleEigengenes(k1, moduleColors1)$eigengenes
MEs11 = orderMEs(MEs10)
#View(MEs11)
MEs1 = MEList1$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDiss1 = 1-cor(MEs1); 
# Cluster module eigengenes 
METree1 = hclust(as.dist(MEDiss1), method = "average"); 
# Plot the result 
sizeGrWindow(7, 6) 
plot(METree1, main = "Clustering of module eigengenes", xlab = "", sub = "")

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM1 = 1-TOMsimilarityFromExpr(k1, power = 7);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM1 = dissTOM1^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM1) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM1, geneTree1, moduleColors1, main = "Network heatmap plot for Sz genes for stage One")

MEs1 = moduleEigengenes(k1, moduleColors1)$eigengenes
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MEs1, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEs1, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MEs1, "Eigengene adjacency heatmap for developmental stage One", marHeatmap = c(3,4,2,2),
                     plotDendrograms = FALSE, xLabelsAngle = 9)

moduleEigengenes(k1, 
                moduleColors1, 
                impute = TRUE, 
                nPC = 1, 
                align = "along average", 
                excludeGrey = FALSE, 
                grey = if (is.numeric(moduleColors1)) 0 else "grey",
                subHubs =TRUE,
                softPower = 7,
                scale = TRUE,
                verbose = 0, indent = 0)

chooseTopHubInEachModule(
  k1, 
  moduleColors1, 
  omitColors = "grey", 
  power = 7, 
  type = "unsigned")

geneModuleMembership1 = as.data.frame(cor(k1, MEs11, use = "p"))

# Sun Mar 14 10:40:49 2021 ------------------------------
# this is the anRichment part

GOenr1 = enrichmentAnalysis(
  classLabels = moduleColors1, identifiers = allLLIDs$LOCUSLINK_ID,#moduleColors should be moduleColors1
  refCollection = GOcollection,
  useBackground = "allOrgGenes",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")


GO_per_set1 <- cbind(GOenr1$enrichmentTable$class, 
                     GOenr1$enrichmentTable$dataSetID,  
                     GOenr1$enrichmentTable$dataSetName, 
                     GOenr1$enrichmentTable$FDR)

tail(GO_per_set1)

colnames(GO_per_set1) <- c("Module", "GO term", "GO process", "FDR")

# Pull the top (most significant) for each module

top_GO_per_module1 <- GO_per_set1[!duplicated(GO_per_set1[,1]), ]
colnames(top_GO_per_module1) <- c("Module", "GO term", "GO process", "FDR")

# Recalculate module eigengenes
MEs1 = moduleEigengenes(k1, moduleColors1)$eigengenes

plotEigengeneNetworks(MEs1, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 

write.csv(GO_per_set1, file="GO_per_set1.csv")

GO_per_set1 <- as.data.frame(GO_per_set1)
#Blue
GO_per_set1_Blue <- GO_per_set1[GO_per_set1$Module == "blue",]
write.csv(GO_per_set1_Blue, file="GO_per_set1_blue.csv")
#Brown
GO_per_set1_Brown <- GO_per_set1[GO_per_set1$Module == "brown",]
write.csv(GO_per_set1_Brown, file="GO_per_set1_brown.csv")
#Turquoise
GO_per_set1_Turquoise <- GO_per_set1[GO_per_set1$Module == "turquoise",]
write.csv(GO_per_set1_Turquoise, file="GO_per_set1_Turquoise.csv")


pdf("plotEigeneNetworks.pdf")
plotEigengeneNetworks(MEs1, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90) 
dev.off()				  
#Cytoscape

geneInfoALL<-data.frame(geneIDs = colnames(k1), moduleColor=moduleColors1)

blueModule<-subset(geneInfoALL, moduleColor == "blue", select = c("geneIDs", "moduleColor"))
write.csv(blueModule, "blueModule.csv")

geneModuleMembership1 = as.data.frame(cor(k1, MEs11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[blueModule$geneIDs,]->geneModuleMembership1_blue #get modulemembership for genes in blue module
geneModuleMembership1_blue$MEblue->geneModuleMembershipblue #correlation of each gene with the with module eigengene for genes in blue module

k1[,blueModule$geneIDs] #get expression values for genes in blue module only

cor((k1[,blueModule$geneIDs]))->cor3 #correlation of each gene expression value with every other gene in the blue module


rownames(cor3)->vec3 # list of gene names in blue module
exportNetworkToCytoscape(cor3, altNodeNames =vec3, nodeAttr=geneModuleMembershipblue)->cytoblue1
write.csv(cytoblue1$edgeData, "cyto_blueEDGE1.csv",quote=FALSE)
write.csv(cytoblue1$nodeData, "cyto_blueNODE1.csv",quote=FALSE)
save(cytoblue1, file="cytoblue1.Rdata")

pacman::p_load(tidyverse, httr, readxl, janitor,dplyr)


url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
df <- read_excel(temp_file, sheet = 4, skip = 3) # reads into a dataframe. First six rows of the excel file are just header
unlink(temp_file)     # deletes the temporary file

pf <- df %>%
  clean_names() %>%
  separate_rows(gene_s_tagged) %>%
  dplyr::select(genes = gene_s_tagged, p_value) %>%
  distinct()

datablue1 <- read.csv("cyto_blueEDGE1.csv")

datablue1$fromNode2 = datablue1$fromNode

names(datablue1)[8] <- "genes"

cytoblue1 <- datablue1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


cytoblue11 = cytoblue1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)
write.csv(cytoblue1, "cyto_blueEDGE1.csv",quote=FALSE)
write.csv(cytoblue11, "cyto_blueEDGE11.csv",quote=FALSE)

#Turq
turquoiseModule<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(turquoiseModule, "turquoiseModule.csv")

geneModuleMembership1 = as.data.frame(cor(k1, MEs11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[turquoiseModule$geneIDs,]->geneModuleMembership1_turquoise#get modulemembership for genes in blue module
geneModuleMembership1_turquoise$MEturquoise->geneModuleMembershipturquoise #correlation of each gene with the with module eigengene for genes in blue module

k1[,turquoiseModule$geneIDs] #get expression values for genes in turquoise module only

cor((k1[,turquoiseModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_turquoise)->cytoturquoise1


write.csv(cytoturquoise1$edgeData, "cyto_turquoiseEDGE1.csv",quote=FALSE)
write.csv(cytoturquoise1$nodeData, "cyto_turquoiseNODE1.csv",quote=FALSE)
save(cytoturquoise1, file="cytoturquoise1.Rdata")

turquoiseModule<-subset(geneInfoALL, moduleColor == "turquoise", select = c("geneIDs", "moduleColor"))
write.csv(turquoiseModule, "turquoiseModule.csv")

dataturquoise1 <- read.csv("cyto_turquoiseEDGE1.csv")

dataturquoise1$fromNode2 = dataturquoise1$fromNode

names(dataturquoise1)[8] <- "genes"

cytoturquoise1 <- dataturquoise1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytoturquoise1, "cyto_turquoiseEDGE1.csv",quote=FALSE)

cytoturq11 = cytoturquoise1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytoturq11, "cyto_turqEDGE11.csv",quote=FALSE)
#Brown
geneInfoALL<-data.frame(geneIDs = colnames(k1), moduleColor=moduleColors1)

brownModule<-subset(geneInfoALL, moduleColor == "brown", select = c("geneIDs", "moduleColor"))
write.csv(brownModule, "brownModule.csv")
geneModuleMembership1 = as.data.frame(cor(k1, MEs11, use = "p")) #module membership for all genes all modules
geneModuleMembership1[brownModule$geneIDs,]->geneModuleMembership1_brown#get modulemembership for genes in blue module
geneModuleMembership1_brown$MEturquoise->geneModuleMembershipbrown #correlation of each gene with the with module eigengene for genes in blue module

k1[,brownModule$geneIDs] #get expression values for genes in turquoise module only

cor((k1[,brownModule$geneIDs]))->cor4 #correlation of each gene expression value with every other gene in the turquoise module


rownames(cor4)->vec4 # list of gene names in blue module
exportNetworkToCytoscape(cor4, altNodeNames =vec4, nodeAttr=geneModuleMembership1_brown)->cytobrown1


write.csv(cytobrown1$edgeData, "cyto_brownEDGE1.csv",quote=FALSE)
write.csv(cytobrown1$nodeData, "cyto_brownNODE1.csv",quote=FALSE)
save(cytobrown1, file="cytobrown1.Rdata")

databrown1 <- read.csv("cyto_brownEDGE1.csv")

databrown1$fromNode2 = databrown1$fromNode

names(databrown1)[8] <- "genes"

cytobrown1 <- databrown1 %>% 
  left_join(pf) %>% 
  dplyr::select(genes, p_value, fromNode, toNode, weight, direction, fromAltName, toAltName) 


write.csv(cytobrown1, "cyto_brownEDGE1.csv",quote=FALSE)

cytobrown11 = cytobrown1 %>% mutate(weight = abs(weight)) %>%
  dplyr::filter(weight > 0.7)

write.csv(cytobrown11, "cyto_brownEDGE11.csv",quote=FALSE)

