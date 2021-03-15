pacman::p_load(tidyverse, 
               httr, 
               readxl, 
               ABAData, 
               anRichment,
               GO.db, 
               org.Hs.eg.db)


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

data("dataset_adult")
ABA_genes <- unique(dataset_adult$hgnc_symbol)
good_genes <- all_sz_genes %>% filter(genes %in% ABA_genes)

good_genes_GO_terms <- AnnotationDbi::select(org.Hs.eg.db, 
                           keytype = "SYMBOL", 
                           keys = (as.character(good_genes$genes)), 
                           columns = c("SYMBOL", "GO", "ONTOLOGY"))

ABA_genes_GO_terms <- AnnotationDbi::select(org.Hs.eg.db, 
                                             keytype = "SYMBOL", 
                                             keys = ABA_genes, 
                                             columns = c("SYMBOL", "GO", "ONTOLOGY"))


z <- table(good_genes_GO_terms$GO)# /length(good_genes$genes)
z1 <- as.data.frame(z) %>% arrange(desc(Freq))
z2 <- table(ABA_genes_GO_terms$GO)/length(ABA_genes)
z3 <- as.data.frame(z2) %>% arrange(desc(Freq))
names(z1) <- c("GO", "CLOZUK_freq")
names(z3) <- c("GO", "ABA_freq")
GO_table <- left_join(z3, z1)
z4 <- ABA_genes_GO_terms %>% dplyr::select(GO, ONTOLOGY) %>% distinct()
GO_table <- left_join(GO_table, z4)

z5 <- AnnotationDbi::select(GO.db, 
                            keys = GO_table$GO, 
                            keytype = "GOID", 
                            columns = c("GOID", "TERM"))
names(z5) <- c("GO", "TERM")
GO_table <- left_join(GO_table, z5)

GO_table %>% ggplot(aes(ABA_freq, CLOZUK_freq, col = ONTOLOGY)) + 
  geom_point() + 
  scale_x_log10() + 
  geom_jitter() +
  scale_y_log10() + 
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_minimal()

GO_table$prob <- pbinom(GO_table$CLOZUK_freq, 103, GO_table$ABA_freq)
(z <- GO_table %>% 
    filter(CLOZUK_freq > 10) %>% 
    arrange(desc(prob)) %>% 
    head(15))
