library(tidyverse)
library(ABAData)
library(factoextra)
library(NbClust)
library(cluster)
library(genes) # from devtools::install_github("eugene100hickey/genes")
library(ggtext)
source("kmeans/NbClust_without_Beale.R")

stage <- 5


if(!exists("dataset_5_stages")){
  data("dataset_5_stages")
}

fviz_nbclust_new <- function(x, print.summary = TRUE,
                         barfill = "steelblue", barcolor = "steelblue")
{
  best_nc <- x$Best.nc
  if(class(best_nc) == "numeric") print(best_nc)
  else if(class(best_nc) == "matrix"){
    best_nc <- as.data.frame(t(best_nc), stringsAsFactors = TRUE)
    best_nc$Number_clusters <- as.factor(best_nc$Number_clusters)
    
    # Summary
    if(print.summary){
      ss <- summary(best_nc$Number_clusters)
      cat ("Among all indices: \n===================\n")
      for(i in 1 :length(ss)){
        cat("*", ss[i], "proposed ", names(ss)[i], "as the best number of clusters\n" )
      }
      cat("\nConclusion\n=========================\n")
      cat("* According to the majority rule, the best number of clusters is ",
          names(which.max(ss)),  ".\n\n")
    }
    df <- data.frame(Number_clusters = names(ss), freq = ss, stringsAsFactors = TRUE )
    max_cluster <- max(as.numeric(df$Number_clusters))
    suppressWarnings(df <- df %>% 
      mutate(Number_clusters = fct_relevel(Number_clusters, 
                                           0:max_cluster %>% as.character)))
    p <- ggpubr::ggbarplot(df,  x = "Number_clusters", y = "freq", fill = barfill, color = barcolor)+
      labs(x = "Number of clusters k", y = "Frequency among all indices",
           title = paste0("Optimal number of clusters - k = ", names(which.max(ss)) ))
    
    return(p)
  }
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
  
  z <- NbClust_without_Beale(wide_allen_sz, 
               diss=diss_matrix, 
               distance = NULL, 
               min.nc=2, 
               max.nc=6,
               method = "ward.D2", 
               index = "all")   
  fviz_nbclust_new(z) + 
    labs(subtitle = glue::glue("NbClust Analysis for Stage <i style = 'color: darkred;'>{stage}</i>"),
         caption = " Malika Charrad, Nadia Ghazzali, Veronique Boiteau, Azam Niknafs (2014). NbClust: An
    R Package for Determining the Relevant Number of Clusters in a Data Set. Journal of
    Statistical Software, 61(6), 1-36. URL http://www.jstatsoft.org/v61/i06/.") +
    theme_minimal() +
    theme(plot.subtitle = element_markdown())
}


NbClust_by_stage(stage = stage)

