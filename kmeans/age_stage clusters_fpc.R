pacman::p_load(ABAEnrichment, tidyverse, httr,
               ABAData, readxl, fpc, factoextra)

url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
df <- read_excel(temp_file, sheet =5 , skip = 6) # reads into a dataframe. First six rows of the excel file are just header
unlink(temp_file)     # deletes the temporary file

all_sz_genes <- df$"Gene(s) tagged" %>%
  # note, there should be inverted commas around Gene(s) tagged
  str_split(",") %>%
  unlist() %>%
  str_trim() %>% 
  as.data.frame() %>%
  distinct()
names(all_sz_genes) <- "genes"


data(dataset_5_stages)
number_centres_vec <- c(2, 3, 3, 2, 2)

stage <- 5

number_centres <- number_centres_vec[stage]

# to see clozuk genes in 5 stages
sz_stage <- dataset_5_stages %>%
  filter(hgnc_symbol %in% all_sz_genes$genes) %>% 
  filter(age_category == stage) %>% 
  select(-c(entrezgene, ensembl_gene_id, age_category)) %>% 
  pivot_wider(names_from = structure, values_from = signal) %>%         
  column_to_rownames(var = "hgnc_symbol") %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame()


k <- kmeans(sz_stage, centers = number_centres, nstart = 25)
kplot <- fviz_cluster(k, sz_stage, repel = T) +
  theme_minimal() +
  scale_color_discrete(labels = c(glue::glue("{round(k$withinss[1] / k$totss *100, 1)}%"),
                                  glue::glue("{round(k$withinss[2] / k$totss *100, 1)}%"),
                                  glue::glue("{round(k$withinss[3] / k$totss *100, 1)}%"))) +
  scale_shape_discrete(guide = F) +
  scale_fill_discrete(guide = F) +
  labs(title = glue::glue("kmeans - stage {stage}"),
       subtitle = glue::glue("Total Variance = {round(k$betweenss / k$totss *100, 1)}%")) +
  theme(text = element_text(size = 18, family = "Ink Free"))

ggsave(filename = glue::glue("kmeans/stage{stage}/kplot{stage}.png"))
