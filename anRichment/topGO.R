pacman::p_load(tidyverse, httr, readxl, ABAData, topGO)

url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
df <- read_excel(temp_file, sheet = 4, skip = 3) # reads into a dataframe. First six rows of the excel file are just header
unlink(temp_file)     # deletes the temporary file

##################################################################
# makes sz_genes, a dataframe with a single column of the CLOZUK genes
#######################################################################

sz_genes <- df$`Gene(s) tagged` %>% 
  str_split(",") %>% 
  unlist() %>% 
  str_trim() %>% 
  as.data.frame() %>% 
  distinct()
names(sz_genes) <- "genes"

data("dataset_5_stages")

# CLOZUK_GENE <- read_tsv("/home/eugene/Desktop/Academic/Research/SZ_eugene/Pardinas_Genes/GWAS-Genes/CLOZUK_GENE.txt")
# PGC2_GENE <- read_tsv("/home/eugene/Desktop/Academic/Research/SZ_eugene/Pardinas_Genes/GWAS-Genes/PGC2_GENE.txt")
# sz_genes <- tibble(genes = c(CLOZUK_GENE$Gene, PGC2_GENE$Gene) %>% unique())

sz_genes <- read_csv("cytoscape/stage1/cyto_greenyellowNODE1.csv") %>% 
  dplyr::select(genes = nodeName)

genes <- sz_genes %>% 
  semi_join(dataset_5_stages, by = c("genes" = "hgnc_symbol")) %>% 
  pull(genes)

geneList <- rep(0.04, length.out = length(genes))
names(geneList) <- genes

sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")


allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

#Defined colMap, ref. https://github.com/Bioconductor-mirror/topGO/blob/master/vignettes/topGO.Rnw
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
      elim = pValue.elim[sel.go],
      classic = pValue.classic[sel.go])
par(cex = 0.4)
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
