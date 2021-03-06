pacman::p_load(tidyverse, httr, readxl, ABAData, here, genes)

locusLinked <- read_csv(here("GPL21185-21174.csv"))

locusLinked <- locusLinked[,1:16]

# url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
# GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
# df <- read_excel(temp_file, sheet = 5, skip = 6) # reads into a dataframe. First six rows of the excel file are just header
# unlink(temp_file)     # deletes the temporary file
# 
# 
# ##################################################################
# # makes sz_genes, a dataframe with a single column of the CLOZUK genes
# #######################################################################
# 
# all_sz_genes <- df$`Gene(s) tagged` %>% 
#   str_split(",") %>% 
#   unlist() %>% 
#   str_trim() %>% 
#   as.data.frame() %>% 
#   distinct()
all_sz_genes <- pardinas()
names(all_sz_genes) <- "GENE_SYMBOL"

data("dataset_5_stages")
#all_sz_genes <- all_sz_genes %>% filter(GENE_SYMBOL %in% dataset_adult$hgnc_symbol)

ABA_genes <- dataset_5_stages %>% 
  dplyr::select(hgnc_symbol) %>% 
  distinct()

ABA_sz_genes <- all_sz_genes %>% 
  filter(GENE_SYMBOL %in% ABA_genes$hgnc_symbol)

ABA_locuslinked <- ABA_sz_genes %>% 
  left_join(locusLinked) %>% 
  select(GENE_SYMBOL, LOCUSLINK_ID) %>% 
  distinct()

write_csv(ABA_locuslinked, "locuslinkedIDs.csv")
