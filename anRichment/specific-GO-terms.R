pacman::p_load(tidyverse, 
               GO.db, 
               org.Hs.eg.db)

p_threshold <- 0.05

# Need to change this depending on what stage
directory <- "anRichment/GO_per_set_One"

# Need to change this depending on what module
filename <- "GO_per_setOne_brown.csv"

anrichment_module <- read_csv(file = glue::glue("{directory}/{filename}")) %>% 
  janitor::clean_names() %>% 
  dplyr::rename(GOID = go_term)

search_GO <- function(index = 1, df = anrichment_module){
  case_when(
    df$ONTOLOGY[index] == 'BP' ~ 
      ifelse(sum(GOBPOFFSPRING[[df$GOID[index]]] %in% df$GOID) == 0, 'yes', 'no'),
    df$ONTOLOGY[index] == 'MF' ~ 
      ifelse(sum(GOMFOFFSPRING[[df$GOID[index]]] %in% df$GOID) == 0, 'yes', 'no'),
    df$ONTOLOGY[index] == 'CC' ~ 
      ifelse(sum(GOCCOFFSPRING[[df$GOID[index]]] %in% df$GOID) == 0, 'yes', 'no')
  )
}

anRich_ontology <- function(filename) {
  anrichment_module <- read_csv(filename) %>% 
    janitor::clean_names() %>% 
    dplyr::rename(GOID = go_term)
  
  
  anrichment_module <- AnnotationDbi::select(GO.db, 
                                             keys = anrichment_module$GOID, 
                                             keytype = "GOID", 
                                             columns = c("GOID", "ONTOLOGY")) %>% 
    right_join(anrichment_module)
  
  anrichment_module <- anrichment_module %>% 
    mutate(end_node = map(1:nrow(anrichment_module), search_GO) %>% unlist(),
           terminal_node = ifelse(GOID %in% c(GO_CC_child, GO_MF_child, GO_BP_child),
                                  'yes', 'no'))
  
  anrichment_children <- anrichment_module %>% dplyr::filter(end_node == 'yes')
  
  # gets the definitions and terms for all our GO's
  GO_terms_background <- AnnotationDbi::select(GO.db, 
                                               keys = anrichment_children$GOID, 
                                               keytype = "GOID", 
                                               columns = c("DEFINITION", "GOID", "ONTOLOGY", "TERM"))
  
  anrichment_children <- GO_terms_background %>% 
    left_join(anrichment_children) 
  
  anrichment_children
}



# GO_CC_child is vector of all CC GO terms that have no
# children of their own, lowest level terms
GO_CC <- as.list(GOCCCHILDREN)
GO_CC_child <- GO_CC[is.na(GO_CC)] %>% names()

# GO_MF_child is vector of all MF GO terms that have no
# children of their own, lowest level terms
GO_MF <- as.list(GOMFCHILDREN)
GO_MF_child <- GO_MF[is.na(GO_MF)] %>% names()

# GO_BP_child is vector of all BP GO terms that have no
# children of their own, lowest level terms
GO_BP <- as.list(GOBPCHILDREN)
GO_BP_child <- GO_BP[is.na(GO_BP)] %>% names()


anrichment_module <- AnnotationDbi::select(GO.db, 
                                           keys = anrichment_module$GOID, 
                                           keytype = "GOID", 
                                           columns = c("GOID", "ONTOLOGY")) %>% 
  right_join(anrichment_module)

anrichment_module <- anrichment_module %>% 
  mutate(end_node = map(1:nrow(anrichment_module), search_GO) %>% unlist(),
         terminal_node = ifelse(GOID %in% c(GO_CC_child, GO_MF_child, 
                                            GO_BP_child),
                                'yes', 'no'))

anrichment_children <- anrichment_module %>% dplyr::filter(end_node == 'yes')

# gets the definitions and terms for all our GO's
GO_terms_background <- AnnotationDbi::select(GO.db, 
                                             keys = anrichment_children$GOID, 
                                             keytype = "GOID", 
                                             columns = c("DEFINITION", "GOID", "ONTOLOGY", "TERM"))

anrichment_children <- GO_terms_background %>% 
  left_join(anrichment_children) %>% 
  dplyr::select(GOID, DEFINITION, ONTOLOGY, TERM, module, fdr, genes)

#Creates new file with "specific-" pre-pending the original filename, in the same directory
write_csv(anrichment_children, file = glue::glue("{directory}/specific-{filename}"))
