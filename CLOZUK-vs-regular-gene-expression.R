library(tidyverse)
library(ABAData)
library(genes) # from devtools::install_github("eugene100hickey/genes")

stage <- 1


if(!exists("dataset_5_stages")){
  data("dataset_5_stages")
}

datset_5_stages <- dataset_5_stages %>% 
  mutate(pardinas = hgnc_symbol %in% pardinas()$genes,
         log_signal = log(signal))

datset_5_stages %>% 
  ggplot(aes(y = log_signal, x = pardinas, colour = pardinas)) + 
  geom_boxplot(notch = T,varwidth = T, show.legend = F) + 
  scale_color_manual(values = c("darkolivegreen", "firebrick4")) +
  theme_minimal()

t.test(signal ~ pardinas, data = datset_5_stages)
