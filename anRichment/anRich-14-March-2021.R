

# Tue May 12 11:15:45 2020 ------------------------------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse, anRichment, org.Hs.eg.db) 

stage <- 5 

if(!exists("GOcollection")){
  GOcollection = buildGOcollection(organism = "human")
}


allLLIDs <- read_csv("locuslinkedIDs.csv")

moduleColors1 <- read_csv(glue::glue("anRichment/stage{stage}/module_colours.csv"))

GOenr1 <- enrichmentAnalysis(
  classLabels = moduleColors1, 
  identifiers = allLLIDs$LOCUSLINK_ID,#moduleColors should be moduleColors1
  refCollection = GOcollection,
  useBackground = "allOrgGenes",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")


GO_per_set1 <- tibble(Module = GOenr1$enrichmentTable$class,
                      `GO term` = GOenr1$enrichmentTable$dataSetID,
                     `GO process` = GOenr1$enrichmentTable$dataSetName,
                     FDR = GOenr1$enrichmentTable$FDR)


# Pull the top (most significant) for each module


module_csv <- function(colour){
  filename <- glue::glue("anRichment/stage{stage}/GO_per_set1_{colour}.csv")
  write_csv(x = GO_per_set1[GO_per_set1$Module == colour,], file = filename)
}

map(table(moduleColors1) %>% names(), module_csv)
