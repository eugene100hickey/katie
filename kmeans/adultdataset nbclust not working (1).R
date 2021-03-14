# Tue May 12 11:15:45 2020 ------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr, readxl, ABAData, tidyverse, ABAEnrichment, ggridges, NbClust, factoextra)  

url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
df <- read_excel(temp_file, sheet = 4, skip = 3) # reads into a dataframe. First six rows of the excel file are just header
unlink(temp_file)     # deletes the temporary file


##################################################################
# makes sz_genes, a dataframe with a single column of the CLOZUK genes
#######################################################################

all_sz_genes <- df$`Gene(s) tagged` %>% 
  str_split(",") %>% 
  unlist() %>% 
  as.data.frame() %>% 
  distinct()
names(all_sz_genes) <- "genes"

library(ABAData)
data("dataset_adult")

pacman::p_load(httr, readxl, ABAData, tidyverse, ABAEnrichment, ggridges, NbClust, factoextra)


sz_genes <- df$`Gene(s) tagged` %>% 
  str_split(",") %>% 
  unlist() %>% 
  as.data.frame() %>% 
  distinct()
names(sz_genes) <- "genes"

data("dataset_adult")


########################################################################################
# makes wide_allen_sz, a dataframe which has sz_genes as rownames, brain areas as columns
# it's wide format and scaled
########################################################################################

wide_allen_sz <-
  dataset_adult %>% 
  select(hgnc_symbol, structure, signal) %>% 
  filter(hgnc_symbol %in% sz_genes$genes) %>% 
  distinct() %>% 
  spread(key = structure, value = signal) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "hgnc_symbol") %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame()

stage<-wide_allen_sz

selected <- c( "kl", "ch", "hartigan",  "cindex", "db",
               "silhouette", "duda", "pseudot2", 
               "ratkowsky", "ball", "ptbiserial", "gap",
               # "frey",
               "mcclain", "gamma", "gplus", "tau",
               "dunn", "hubert", "sdindex", "dindex", "sdbw")
results <- vector("list",length(selected))
for (i in 1:length(selected)) {
  results[[i]] <- try(NbClust(wide_allen_sz,
                              min.nc=1, max.nc=10,
                              method="ward.D", index=selected[i]))
}
output_nbclust <- map(1:length(selected), function(x) results[[x]]$Best.nc) %>%
  unlist()
output_nbclust <- output_nbclust[names(output_nbclust) == "Number_clusters"] %>%
  as.data.frame()

names(output_nbclust) <- "cluster_number"
output_nbclust %>% ggplot(aes(cluster_number)) + 
  geom_bar(fill = "firebrick4") + 
  scale_x_continuous(breaks = seq(1,10)) +
  labs(title = "Optimal Cluster Number",
       subtitle = "Stage 5",
       caption = "@data from ABA"
  ) +
  xlab("Number of Clusters") + 
  ylab("") +
  theme_minimal() 
