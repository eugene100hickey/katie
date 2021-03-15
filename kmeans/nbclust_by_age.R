library(tidyverse)
library(ABAData)
library(ABAEnrichment)
library(ggridges)
require(httr)
require(readxl)
library(NbClust)
library(factoextra)

url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
df <- read_excel(temp_file, sheet = 5, skip = 6) # reads into a dataframe. First six rows of the excel file are just header
unlink(temp_file)     # deletes the temporary file

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

data("dataset_5_stages")

z <- dataset_5_stages %>% 
  dplyr::select(hgnc_symbol, structure, signal, age_category) %>% 
  filter(hgnc_symbol %in% all_sz_genes$genes) %>% 
  mutate(ages = factor(age_category, levels = 5:1))

# stage corresponds to the age of subjects from ABA
# 1: prenatal
# 2: infant (0-2 yrs)
# 3: child (3-11)
# 4: adolescent (12-19 yrs)
# 5: adult (>19 yrs)
stage <- 3

age_profile <- z %>%
  filter(ages == stage) %>%
  select(-c(age_category, ages)) %>%
  pivot_wider(names_from = structure, values_from = signal) %>%
  column_to_rownames(var = "hgnc_symbol") %>%
  as.matrix() %>%
  t() %>%
  scale() %>% 
  t() %>%
  as.data.frame()

selected <- c( "kl", "ch", "hartigan",  "cindex", "db", 
               "silhouette", "duda", "pseudot2", "beale"#, 
               #"ratkowsky", "ball", "ptbiserial", "gap", 
               # "frey", 
               #"mcclain", "gamma", "gplus", "tau", 
              # "dunn", "hubert", "sdindex", "dindex", "sdbw"
               ) 

if("res.nbclust" %in% list.files(glue::glue("kmeans/stage{stage}"))){
  res.nbclust <- readRDS(glue::glue("kmeans/stage{stage}/res.nbclust"))
}

lista.methods = c("kl", "ch", "hartigan","mcclain", "gamma", "gplus",
                  "tau", "dunn", "sdindex", "sdbw", "cindex", "silhouette",
                  "ball","ptbiserial", "gap","frey")
lista.distance = c("metodo","euclidean", "maximum", "manhattan", "canberra")

tabla = as.data.frame(matrix(ncol = length(lista.distance), nrow = length(lista.methods)))
names(tabla) = lista.distance

for (j in 2:length(lista.distance)){
  for(i in 1:length(lista.methods)){
    
    nb = NbClust(age_profile, distance = lista.distance[j],
                 min.nc = 2, max.nc = 10, 
                 method = "complete", index =lista.methods[i])
    tabla[i,j] = nb$Best.nc[1]
    tabla[i,1] = lista.methods[i]
    
  }}

if(!exists("res.nbclust")){
  res.nbclust <- NbClust(age_profile, distance = "euclidean",
                         min.nc = 2, max.nc = 9,
                         method = "complete", index = selected)
  saveRDS(res.nbclust, glue::glue("kmeans/stage{stage}/nbclust"))
} 

fviz_nbclust(res.nbclust) +
  theme_minimal() +
  ggtitle("NbClust's optimal number of clusters")
  
