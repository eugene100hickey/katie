

# Tue May 12 11:15:45 2020 ------------------------------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(ABAData, httr, readxl, tidyverse, WGCNA, igraph, pins) 

stage <- 5 
softPower <- 7

# library(anRichment)
# library(org.Hs.eg.db)
# GOcollection = buildGOcollection(organism = "human")


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
  filter(age_category == stage) %>% 
  dplyr::select(-c(entrezgene, ensembl_gene_id, age_category)) %>% 
  pivot_wider(names_from = structure, values_from = signal) %>% 
  column_to_rownames(var = "hgnc_symbol")%>%
  t() %>% 
  scale() %>% 
  t()

# Tue May 12 11:25:55 2020 ------------------------------
#sz_stage_4<-k4 <- t(k4)
k1 <- sz_stage_1 %>% t()
# Tue May 12 11:26:03 2020 ------------------------------


# Part 2:  Step-by-step network construction and module detection


#Calculating adjacencies with soft threshold of 6

adjacency = adjacency(k1, power = softPower);

# Turn adjacency into topological overlap 
TOM1 = TOMsimilarity(adjacency); 
dissTOM1 = 1-TOM1

# Call the hierarchical cluste ring function 
geneTree1 = hclust(as.dist(dissTOM1), method = "average"); 

# Calculate eigengenes 
minModuleSize = 10; 
# Module identification using dynamic tree cut: 
dynamicMods1 = cutreeDynamic(dendro = geneTree1, distM = dissTOM1, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize); 
table(dynamicMods1)

# Convert numeric lables into colors 
dynamicColors1 = labels2colors(dynamicMods1) 
table(dynamicColors1) 

# Tue May 12 11:35:46 2020 ------------------------------
MEList1 = moduleEigengenes(sz_stage_1 %>% t(), colors = dynamicMods1)
# Tue May 12 11:35:52 2020 ------------------------------

MEs = MEList1$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDiss1 = 1-cor(MEs); 
# Cluster module eigengenes 
METree1 = hclust(as.dist(MEDiss1), method = "average"); 
MEDissThres1 = 1

merge1 = mergeCloseModules(k1, dynamicColors1, cutHeight = MEDissThres1, verbose = 3) # missing the 1 here
# The merged module colors
mergedColors1 = merge1$colors;
# Eigengenes of the new merged modules:
mergedMEs1 = merge1$newMEs

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
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM1 = 1-TOMsimilarityFromExpr(k1, power = softPower);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM1 = dissTOM1^softPower;
# Set diagonal to NA for a nicer plot
diag(plotTOM1) = NA;
# Call the plot function

MEs1 = moduleEigengenes(k1, moduleColors1)$eigengenes


moduleEigengenes(k1, 
                moduleColors1, 
                impute = TRUE, 
                nPC = 1, 
                align = "along average", 
                excludeGrey = FALSE, 
                grey = if (is.numeric(moduleColors1)) 0 else "grey",
                subHubs =TRUE,
                softPower = softPower,
                scale = TRUE,
                verbose = 0, indent = 0)

hub_genes <- chooseTopHubInEachModule(
  k1, 
  moduleColors1, 
  omitColors = "grey", 
  power = softPower, 
  type = "unsigned")

geneModuleMembership1 = as.data.frame(cor(k1, MEs11, use = "p"))

# Sun Mar 14 10:40:49 2021 ------------------------------
# this is the anRichment part

# GOenr1 = enrichmentAnalysis(
#   classLabels = moduleColors1, identifiers = allLLIDs$LOCUSLINK_ID,#moduleColors should be moduleColors1
#   refCollection = GOcollection,
#   useBackground = "allOrgGenes",
#   threshold = 1e-4,
#   thresholdType = "Bonferroni",
#   getOverlapEntrez = TRUE,
#   getOverlapSymbols = TRUE,
#   ignoreLabels = "grey")
# 
# 
# GO_per_set1 <- cbind(GOenr1$enrichmentTable$class, 
#                      GOenr1$enrichmentTable$dataSetID,  
#                      GOenr1$enrichmentTable$dataSetName, 
#                      GOenr1$enrichmentTable$FDR)
# 
# 
# colnames(GO_per_set1) <- c("Module", "GO term", "GO process", "FDR")
# 
# # Pull the top (most significant) for each module
# 
# top_GO_per_module1 <- GO_per_set1[!duplicated(GO_per_set1[,1]), ]
# colnames(top_GO_per_module1) <- c("Module", "GO term", "GO process", "FDR")

# Recalculate module eigengenes
MEs1 = moduleEigengenes(k1, moduleColors1)$eigengenes


# write.csv(GO_per_set1, file="GO_per_set1.csv")
# GO_per_set1 <- as.data.frame(GO_per_set1)
# Sun Mar 14 13:52:11 2021 ------------------------------
module_colours <- names(hub_genes)

# module_csv <- function(colour){
#   filename <- glue::glue("cytoscape/stage{stage}/GO_per_set1_{colour}.csv")
#   write_csv(GO_per_set1[GO_per_set1$Module == colour,], file = filename)
# }
# 
# map(module_colours, module_csv)

geneInfoALL<-data.frame(geneIDs = colnames(k1), moduleColor=moduleColors1)
geneModuleMembership1 <- as.data.frame(cor(k1, MEs11, use = "p")) #module membership for all genes all modules

colourModules <- function(colour){
  z <- subset(geneInfoALL, moduleColor == colour, select = c("geneIDs", "moduleColor"))
  write.csv(z, glue::glue("cytoscape/stage{stage}/{colour}Module.csv"))
  geneModuleMembership1_blue <- geneModuleMembership1[z$geneIDs,] #get module membership for genes in "colour" module
  names(geneModuleMembership1_blue) <- names(geneModuleMembership1_blue) %>% str_remove("ME")
  geneModuleMembershipblue <- geneModuleMembership1_blue %>%
    dplyr::select({{colour}}) #correlation of each gene with the with module eigengene for genes in blue module
  k1[,z$geneIDs] #get expression values for genes in blue module only
  cor3 <- cor((k1[,z$geneIDs])) #correlation of each gene expression value with every other gene in the blue module
  vec3 <- rownames(cor3) # list of gene names in blue module
  cytoblue1 <- exportNetworkToCytoscape(cor3,
                                        altNodeNames =vec3,
                                        nodeAttr=geneModuleMembershipblue)
  cytoblue1$edgeData <- cytoblue1$edgeData %>% 
    left_join(pf)
  cytoblue11 = cytoblue1$edgeData %>% mutate(weight = abs(weight)) %>%
    dplyr::filter(weight > 0.7)
  write.csv(cytoblue1$edgeData, glue::glue("cytoscape/stage{stage}/cyto_{colour}EDGE{stage}.csv"), quote=FALSE)
  write.csv(cytoblue1$nodeData, glue::glue("cytoscape/stage{stage}/cyto_{colour}NODE{1}.csv"), quote=FALSE)
  save(cytoblue1, file=glue::glue("cytoscape/stage{stage}/cyto{colour}{stage}.Rdata"))
  write.csv(cytoblue11, 
            glue::glue("cytoscape/stage{stage}/cyto_{colour}EDGE{stage}1.csv"),
                       quote=FALSE)
}

pf <- df %>%
  janitor::clean_names() %>%
  separate_rows(gene_s_tagged) %>%
  dplyr::select(fromNode = gene_s_tagged, p_value) %>%
  distinct()

map(module_colours, colourModules)

write_csv(tibble(colours = moduleColors1), 
          glue::glue("anRichment/stage{stage}/module_colours.csv"))
