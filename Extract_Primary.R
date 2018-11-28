#Create a function to read&process all vcf KD/ref pairs
#And merge them into large tibble with new column SampleID

#Attach libraries
library(tidyverse)


#Functions
numextract.int <- function(string){ 
  str_extract(string, "\\-*\\d*\\d*")
} 
numextract.double <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

Load.Data <- function(Path, Save){
  #Establish filepaths
  Path.Data <- Path
  Files <- list.files(file.path(Path.Data), pattern = "*.vcf$")
  Meta <- read_csv("FANTOM6_samples.csv", col_names = TRUE)
  #Merging loop
  for (i in 1:length(Files)) {
    #Read vcf
    tmp.vcf <- readLines(file.path(Path.Data, Files[i]))
    tmp.vcf.data <- read.table(file.path(Path.Data, Files[i]))
    #Extract and atach colnames
    tmp.vcf <- tmp.vcf[-(grep("#CHROM", tmp.vcf) + 1):-(length(tmp.vcf))]
    vcf.names <- unlist(strsplit(tmp.vcf[length(tmp.vcf)], "\t"))
    names(tmp.vcf.data)<-vcf.names
    #Extract DP, EF and Metadata
    tmp.vcf.data <- tmp.vcf.data %>% 
      mutate(DP = INFO %>% 
               str_split_fixed(., pattern = ";", n = 2) %>% 
               as.tibble() %>%
               pull(V1) %>%
               numextract.double(.), 
             EF = INFO %>% 
               str_split_fixed(., pattern = ";", n = 2) %>% 
               as.tibble() %>%
               pull(V2) %>%
               numextract.double(.),
             EDIT = paste(REF, ALT, sep = " to "),
             File.Name = rep(paste(Files[i]), length(.$`#CHROM`))) %>%
      mutate(SampleID = File.Name %>% 
               str_sub(., 1, -5), 
             EDIT = paste(REF, ALT, sep = " to ")
      ) %>%
      select(-INFO, -FORMAT, -File.Name, -starts_with("RDhi", ignore.case = FALSE))
    if (i == 1) {
      Data <- tmp.vcf.data
      rm(tmp.vcf.data)
    } else {
      if (i > 1 & i < length(Files)) {
        Data <- rbind(Data, tmp.vcf.data) 
      } else {
        if (i == length(Files)) {
          Data <- rbind(Data, tmp.vcf.data)
          #Create Factorial EF variable
          Data <- Data %>% 
            mutate(EF.Factor = EF)
          Data[Data$EF <= 0.25, ]$EF.Factor <- rep("4.LowEf", length(Data[Data$EF <= 0.25, ]$EF.Factor))
          Data[Data$EF > 0.25 & Data$EF <= 0.5, ]$EF.Factor <- rep("3.MidEf Left", length(Data[Data$EF > 0.25 & Data$EF <= 0.5, ]$EF.Factor))
          Data[Data$EF > 0.5 & Data$EF <= 0.75, ]$EF.Factor <- rep("2.MidEf Right", length(Data[Data$EF > 0.5 & Data$EF <= 0.75, ]$EF.Factor))
          Data[Data$EF > 0.75, ]$EF.Factor <- rep("1.HighEf", length(Data[Data$EF > 0.75, ]$EF.Factor))
          #Attach metadata
          Data <- Data %>%
            inner_join(., Meta, by = c("SampleID" = "Sample ID"))
          #Substitute NA for 'Control'
          Data$`KD gene ID`[is.na(Data$`KD gene ID`)] <- "Control"
          return(Data)
          if(Save == TRUE){
            write_tsv(Data, file.path(Path.Data, "vcf_comp.tsv"), col_names = TRUE)
          }
          rm(tmp.vcf.data, vcf.names, tmp.vcf, i)
        }
      }
    }
  }
}




Data <- Data %>% 
  mutate(EF.Factor = EF)
Data[Data$EF <= 0.25, ]$EF.Factor <- rep("4.LowEf", length(Data[Data$EF <= 0.25, ]$EF.Factor))
Data[Data$EF > 0.25 & Data$EF <= 0.5, ]$EF.Factor <- rep("3.MidEf Left", length(Data[Data$EF > 0.25 & Data$EF <= 0.5, ]$EF.Factor))
Data[Data$EF > 0.5 & Data$EF <= 0.75, ]$EF.Factor <- rep("2.MidEf Right", length(Data[Data$EF > 0.5 & Data$EF <= 0.75, ]$EF.Factor))
Data[Data$EF > 0.75, ]$EF.Factor <- rep("1.HighEf", length(Data[Data$EF > 0.75, ]$EF.Factor))

#Get Metadata
Data.Meta <- read_csv("FANTOM6_samples.csv", col_names = TRUE)

#Attach meta
Data <- Data %>%
  inner_join(., Data.Meta, by = c("SampleID" = "Sample ID"))
Data$`KD gene ID`[is.na(Data$`KD gene ID`)] <- "Control"

#Establish filepaths
Path.Data <- "Data/"
Files <- list.files(file.path(Path.Data), pattern = "*.vcf$")

#Merging loop
for (i in 1:length(Files)) {
    #Read vcf
    tmp.vcf <- readLines(file.path(Path.Data, Files[i]))
    tmp.vcf.data <- read.table(file.path(Path.Data, Files[i]))
    #Extract and atach colnames
    tmp.vcf <- tmp.vcf[-(grep("#CHROM", tmp.vcf) + 1):-(length(tmp.vcf))]
    vcf.names <- unlist(strsplit(tmp.vcf[length(tmp.vcf)], "\t"))
    names(tmp.vcf.data)<-vcf.names
    #Extract DP, EF and Metadata
    tmp.vcf.data <- tmp.vcf.data %>% 
      mutate(DP = INFO %>% 
               str_split_fixed(., pattern = ";", n = 2) %>% 
               as.tibble() %>%
               pull(V1) %>%
               numextract.double(.), 
             EF = INFO %>% 
               str_split_fixed(., pattern = ";", n = 2) %>% 
               as.tibble() %>%
               pull(V2) %>%
               numextract.double(.),
             EDIT = paste(REF, ALT, sep = " to "),
             File.Name = rep(paste(Files[i]), length(.$`#CHROM`))) %>%
      mutate(SampleID = File.Name %>% 
               str_sub(., 1, -5), 
             EDIT = paste(REF, ALT, sep = " to ")
      ) %>%
      select(-INFO, -FORMAT, -File.Name, -starts_with("RDhi", ignore.case = FALSE))
    if (i == 1) {
      Data <- tmp.vcf.data
      rm(tmp.vcf.data)
    } else {
      if (i > 1 & i < length(Files)) {
        Data <- rbind(Data, tmp.vcf.data) 
      } else {
        if (i == length(Files)) {
          Data <- rbind(Data, tmp.vcf.data)
          #Write for later
          write_tsv(Data, file.path(Path.Data, "vcf_comp.tsv"), col_names = TRUE)
          return(Data)
          rm(tmp.vcf.data, vcf.names, tmp.vcf, i)
        }
      }
    }
}


