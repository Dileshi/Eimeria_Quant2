## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Quant2
## Script adapted from Victor's 1_Data_preparation

library(tidyverse)
library(dplyr)
library(tidyr)
library(janitor)
library(readr)


## Quant2_E57 includes both genotype and flotation data merged unlike in Victor's case
#Hence, Victor's R script is adapted accordingly 

##Load data
#Quant2_E57 data has flotation and genotype data
#Nanadrop_2.0 has the DNA concentration data of our extractions
##E88_fecweight has the weight of fecal matter used in the DNA extraction
sample.data <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Quant2_E57.csv")
Nanodrop <- read.csv ("Raw_Data/E88_primary_Nanodrop.csv")
fecweight <- read.csv ("Raw_Data/E88_primary_fecweight.csv", header=T, na.strings=c("","NA"))

## Filter out the E64 (E. ferrisi) mice in Quant2_E57 file
sample.data <- sample.data %>% 
  dplyr::filter(!primary_infection == "E64")
#remove first row
sample.data <- subset( sample.data, select = -1 )

##rename column 'DNA_conc' of Nanadrop_2.0 to Conc_Data to allow cohesion with Victor's data
Nanodrop <- Nanodrop %>% 
  rename(
    Conc_DNA = DNA_Conc,
    R260.280 = X260.280,
    R260.230 = X260.230
    )

##rename the entries under the column 'Sample_ID' in file Nanodrop to 
##E57aXXX to allow merge with file Quant2_E57
Nanodrop$labels = paste('E57a',Nanodrop$Sample_ID, sep = '')
##remove Nanodrop duplicates based on 260/280 ratio
Nanodrop$Sample_ID[duplicated(Nanodrop$Sample_ID)]
Nanodrop <- Nanodrop %>% distinct(Sample_ID, .keep_all = TRUE)
#remove column Sample_ID
Nanodrop = subset(Nanodrop, select = -c(Sample_ID, Unit) )

##rename column 'feces_weight' of Quant2_E57 to 'fecweight_flot' to allow cohesion with Victor's data
sample.data <- sample.data %>% 
  rename(fecweight_flot = feces_weight)

##merge files Quant2_E57, Nanadrop and fecweight
sample.data = left_join(sample.data, Nanodrop) %>%
  left_join(., fecweight)

##remove column 'Sample_ID' from merged file
#sample.data = subset(sample.data, select = -c(Sample_ID) )


##R identifies the entries under column fecweight_DNA as 'NULL'
#Hence, class changed to numeric
sample.data$fecweight_DNA<-sapply(sample.data$fecweight_DNA, as.numeric)

###Use function from Alice Balard to calculate OPG
    #volume of each large square (1mm x 1mm x 0.1mm height = 0.1mm3)
    #0.1mm3 = 0.1µl (x10000 when converted into ml)
    #https://www.emsdiasum.com/microscopy/technical/datasheet/68052-14.aspx <- explanation 
calculateOPG <- function(sample.data){
  sample.data$mean_Neubauer <- 
    (sample.data$oocyst_sq1 + sample.data$oocyst_sq2 + sample.data$oocyst_sq3 + sample.data$oocyst_sq4) / 4
                                           # NB! Limit of detection = 1 oocysts
  sample.data$mean_Neubauer[sample.data$oocyst_sq1 + sample.data$oocyst_sq2 + sample.data$oocyst_sq3 + sample.data$oocyst_sq4 == 1] <- 0
  sample.data$oocysts.per.tube <- sample.data$mean_Neubauer * 10000 * sample.data$dilution
  sample.data$OPG <- sample.data$oocysts.per.tube / sample.data$fecweight_flot
  ## If we don't have the fecal weight BUT we counted in Neubauer chamber 0, then OPG = 0
  sample.data$oocysts.per.tube[sample.data$fecweight_flot == 0 & sample.data$mean_Neubauer == 0] <- 0
  sample.data$OPG[sample.data$fecweight_flot == 0 & sample.data$mean_Neubauer == 0] <- 0
  return(sample.data)
}

sample.data <- calculateOPG(sample.data = sample.data)

##Check for spaces
sample.data$EH_ID <- gsub(pattern = " ", replacement = "", x = sample.data$EH_ID)



##remove?
###Estimate microbiota density (MD) from Contijoch et al. 2019
### MD = total DNA per sample (µg)/ mg of fresh feces
##considering 40µL of elution volume 
#Conc_DNA is in ng/µl, hence, x0.001 (fecal mass measured in grams)
#sample.data %>%
  #mutate(Total_DNA = (sample.data$Conc_DNA*40)*0.001) -> sample.data ### Add a new variable that will contain total DNA extracted per sample in µg

#sample.data %>%
  #mutate(Microbial_density = sample.data$Total_DNA/(sample.data$fecweight_DNA*1000)) -> sample.data ### Total DNA extracted per sample in µg by feces weight in mg

#sample.data$labels<- as.vector(sample.data$labels)
#rownames(sample.data) <- make.unique(sample.data$labels)

rm(Nanodrop,fecweight, calculateOPG)

##Merge qPCR data 
#Ms. Fay Webster's function was adpated to merge all qPCR files into one
#setwd("Raw_Data/qPCR/qPCR_Primary/")

#all individual qPCR.csv files form a list
#list_faeces <- as.list(list.files())

#to perform a function to a list, it is transformed into a vector
#list_names <- as.vector(unlist(list_faeces))

#function set by Fay Webster to deleted rows 1 to 24 and set a new column that 
#includes the file name i.e. date of qPCR experiment
#read_qPCR_file <- function(x) {
  df1 = read_csv(x)
  filename <- colnames(df1[2])
  df1 <- df1 %>%
    filter(!row_number() %in% c(1:24))
  colnames(df1) <- df1[1, ]
  df1 <- df1 %>% filter(!row_number() %in% 1)
  df1 <- df1 %>% dplyr::mutate(plate = filename)
#}

#list_results <- lapply(list_names, read_qPCR_file)

#bind/merge all indivual csv files into one
#df_results <- Reduce(rbind, list_results)

#removes duplicates
#df_results <- unique(df_results)  

#rm(list_results, list_faeces,list_names, calculateOPG, read_qPCR_file)

#setwd("/Users/vinuri/Documents/GitHub/Eimeria_Quant2/Output_Data/")
#write.csv(df_results, "qPCR_E88_primary_merged.csv", row.names=FALSE)
#rm(df_results)


