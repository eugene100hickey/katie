library(tidyverse)
library(ABAData)
library(factoextra)
library(NbClust)
library(cluster)
library(genes) # from devtools::install_github("eugene100hickey/genes")
library(ggtext)

stage <- 4

if(!exists("dataset_5_stages")){
  data("dataset_5_stages")
}

NbClust_by_stage <- function(stage = 1) {
  wide_allen_sz <-
    dataset_5_stages %>% 
    dplyr::filter(age_category == stage) %>% 
    dplyr::select(hgnc_symbol, structure, signal) %>% 
    dplyr::filter(hgnc_symbol %in% pardinas()$genes) %>% 
    distinct() %>% 
    pivot_wider(names_from = "structure", 
                values_from = "signal") %>% 
    remove_rownames() %>%
    column_to_rownames(var = "hgnc_symbol")
  
  diss_matrix<- dist(wide_allen_sz, method = "euclidean", diag=FALSE)
  
  z <- NbClust(wide_allen_sz, 
               diss=diss_matrix, 
               distance = NULL, 
               min.nc=2, 
               max.nc=20,
               method = "ward.D2", 
               index = "all")   
  fviz_nbclust(z) + 
    labs(subtitle = glue::glue("NbClust Analysis for Stage <i style = 'color: darkred;'>{stage}</i>"),
         caption = " Malika Charrad, Nadia Ghazzali, Veronique Boiteau, Azam Niknafs (2014). NbClust: An
    R Package for Determining the Relevant Number of Clusters in a Data Set. Journal of
    Statistical Software, 61(6), 1-36. URL http://www.jstatsoft.org/v61/i06/.") +
    theme_minimal() +
    theme(plot.subtitle = element_markdown())
}


NbClust_by_stage(stage = stage)
