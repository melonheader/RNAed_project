##
require(tidyverse)
require(FactoMineR)
require(factoextra)
##
#Data <- read_csv("marshfield_world.csv", col_names = TRUE) #For fun and knowledge
##
Data <- read_tsv("vcf_comp.tsv", col_names = TRUE)
Data.Pca <- Data %>% 
  mutate(ed.id = paste(`#CHROM`, POS, EDIT, sep = "_"), 
         samp.id = paste(`KD gene ID`, SampleID, `Cell type ID`, sep = "_")) %>%
  select(ed.id, samp.id) %>% 
  group_by(ed.id, samp.id) %>%
  summarise(counts = n()) %>%
  spread(key = ed.id, value = counts, fill = "0")
##
RowNames <- Data.Pca %>% pull(samp.id)
Data.Pca <- as.matrix(Data.Pca)
Data.Pca <- Data.Pca[, -1]
rownames(Data.Pca) <- RowNames  
##
Res.Mca <- MCA(Data.Pca, graph = FALSE)
save(Res.Mca, file = "Res.Mca.Rds")