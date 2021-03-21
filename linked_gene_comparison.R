if (!require("pacman")) install.packages("pacman")
pacman::p_load(ABAData, httr, readxl, tidyverse, patchwork, ggplotify, janitor, ggrepel, ggtext, glue, dplyr)
theme_set(theme_minimal())

url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx")))

df_3 <- xlsx::read.xlsx(temp_file, sheetName = "Supp Table 3", startRow = 4) %>% 
  clean_names()
df_4 <- xlsx::read.xlsx(temp_file, sheetName = "Supp Table 4", startRow = 8) %>% 
  clean_names()
df_5 <- xlsx::read.xlsx(temp_file, sheetName = "Supp Table 5", startRow = 5) %>% 
  clean_names()
unlink(temp_file)

d4_genes <- df_4$gene_s_tagged %>% 
  str_split(",") %>% 
  unlist() %>% 
  str_trim() %>% 
  as.data.frame() %>% 
  distinct()

d3_genes <- df_3$gene_s_tagged %>% 
  str_split(",") %>% 
  unlist() %>% 
  str_trim() %>% 
  as.data.frame() %>% 
  distinct()
