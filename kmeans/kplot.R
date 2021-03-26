library(tidyverse)
library(ABAData)
library(factoextra)
library(cluster)
library(genes) # from devtools::install_github("eugene100hickey/genes")
library(ggtext)
library(showtext)

font_add(family = "Kalam", regular = "fonts/Kalam/Kalam-Regular.ttf")
showtext::showtext_auto()


stage <- 3
centres <- 3
repel <- TRUE


if(!exists("dataset_5_stages")){
  data("dataset_5_stages")
}

fviz_clust_by_stage <- function(stage = 1, centres = 3, axes = c(1, 2), repel = T) {
  wide_allen_sz <-
    dataset_5_stages %>% 
    dplyr::filter(age_category == stage) %>% 
    dplyr::select(hgnc_symbol, structure, signal) %>% 
    dplyr::filter(hgnc_symbol %in% pardinas()$genes) %>% 
    distinct() %>% 
    pivot_wider(names_from = "structure", 
                values_from = "signal") %>% 
    remove_rownames() %>%
    column_to_rownames(var = "hgnc_symbol") %>%
    t() %>% 
    scale() %>% 
    t() %>% 
    as.data.frame()
  
  clusters <- kmeans(wide_allen_sz, centers = centres, nstart = 25)
  fviz_cluster(clusters, wide_allen_sz, 
               repel = repel, 
               axes = axes, 
               max.overlaps = 50,
               labelsize = 8,
               pointsize = 0.8)
}

fviz_clust_by_stage(stage = stage, centres = centres, axes = c(1,2), repel = repel) + 
  labs(title = glue::glue("Cluster Plot for Stage {stage}"),
       caption = "<i style='font-family: Kalam; color:darkred; font-size: 15px;'>factoextra: Extract and Visualize the Results of Multivariate Data Analyses</i>
       <br>
       <i style='font-family: Kalam; color:darkred; font-size: 15px;'>Alboukadel Kassambara and Fabian Mundt, 2020</i>") +
  theme_minimal() +
  theme(plot.caption = element_markdown(),
        plot.title = element_text(size = 18, family = "Kalam"),
        legend.position = "none",
        )
ggsave(glue::glue("kmeans/stage{stage}/kplot{stage}_centres_{centres}_repel-{repel}.png"))
  