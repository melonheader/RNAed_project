#Plotting Script

#Attach libraries
library(tidyverse)

#Establish filepaths
Path.Data <- "Data/"

#Get data
Data <- read_tsv(file.path(Path.Data, "vcf_comp.tsv"), col_names = TRUE)

#####Temporary workspace#####


#Plot editing events

#SM.Type wise
ggplot() +
  geom_bar(data = Data, aes(x = Data$EDIT, fill = Data$SM.Type)) + 
  theme_bw()
#Set wise
ggplot() +
  geom_bar(data = Data, aes(x = Data$EDIT, fill = Data$Set)) + 
  geom_vline(aes(xintercept = median(weight)), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07") +
  theme_bw()

#Check density
ggplot(data = Data) +
  geom_density(aes(x = EF, y = ..scaled.., fill = `KD gene ID`), alpha = 1/2) +
  geom_vline(aes(xintercept = median(EF)), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07") +
  facet_grid(`Cell type ID` ~ EDIT) +
  theme_bw()

ggplot(data = Data) +
  geom_density(aes(x = DP, y = ..scaled.., fill = Set), alpha = 1/2) +
  facet_grid(SM.Type ~ .) +
  theme_bw()


ggplot(data = Data) +
  geom_violin(aes(x = Set, y = EF, fill = Set), alpha = 1/2) +
  facet_grid(SM.Type ~ .) +
  theme_bw()

#Check editing distributions across chromosomes 
#Attach additional factorial variable
Data <- Data %>% 
  mutate(EF.Factor = EF)

Data[Data$EF < 0.25, ]$EF.Factor <- rep("LowEf", length(Data[Data$EF < 0.25, ]$EF.Factor))
Data[Data$EF >= 0.25 & Data$EF < 0.5, ]$EF.Factor <- rep("MidEf Left", length(Data[Data$EF >= 0.25 & Data$EF < 0.5, ]$EF.Factor))
Data[Data$EF >= 0.5 & Data$EF < 0.75, ]$EF.Factor <- rep("MidEf Right", length(Data[Data$EF >= 0.5 & Data$EF < 0.75, ]$EF.Factor))
Data[Data$EF >= 0.75, ]$EF.Factor <- rep("HighEf", length(Data[Data$EF >= 0.75, ]$EF.Factor))

ggplot(data = Data %>%
         filter(DP > 30)) + 
  geom_bar(aes(x = `#CHROM`, y = ..count.., fill = EF.Factor)) +
  facet_grid(SM.Type ~ EDIT) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))


##Look at the most efficiently edited sites
Data.EEE <- Data %>%
  filter(EF > 0.85 & DP > 30)



#Intersect with bed regions
write_tsv(Data.EEE %>% 
            mutate(POS1 = POS) %>%
            select('#CHROM', POS, POS1),
          "bed/EEE.pos.bed",
          col_names = FALSE)

Data.EEE %>% 
  mutate(POS1 = POS) %>%
  select('#CHROM', POS, POS1)


#Press to save
ggsave("Plots/PCA_SNHG5.png", 
       plot = gg, 
       device = "png", 
       units = "mm",
       width = 300,
       heigh = 160,
       dpi = 500)




#Check intersections 
Data.EEE.intersect <- read_tsv("bed/EEE.intersect.txt", col_names = FALSE)

