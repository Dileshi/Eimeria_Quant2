###script for secondary infection (E64)
##Run scripts on Root of Repo (Eimeria_Quant2)

library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)

data = read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")

challenge <- data[data$EH_ID %in% c("LM0258", "LM0259","LM0260","LM0262",
                                    "LM0265","LM0269","LM0270","LM0279","LM0280", "LM0284"),]
challenge <- challenge %>% 
  dplyr::filter(infection == "challenge")
##write.csv(challenge,"E64_challenged.csv")

#load fecal_weight and Nanodrop data 
fecweight <- read.csv("Raw_Data/E64_challenged_fecweight.csv")
Nanodrop <- read.csv("Raw_Data/E64_challenged_Nanodrop.csv")

#remove and rename columns 
fecweight <- fecweight[-c(1)]
Nanodrop <- Nanodrop[-c(1,3)]
Nanodrop %>% rename(
    Conc_DNA = Nucleic.Acid.Conc...ng.µl.,
    labels = Sample.ID) -> Nanodrop
Nanodrop$labels = paste('E57by',Nanodrop$labels, sep = '')

#merge fecal data with main challenge 
challenge = left_join(challenge, fecweight)

#rename mislabelled samples in nanodrop
Nanodrop$labels[Nanodrop$labels == "E57byLKT"] <- "E57byLNT"
Nanodrop$labels[Nanodrop$labels == "E57byDSZ"] <- "E57byOSZ"
Nanodrop$labels[Nanodrop$labels == "E57byAVM"] <- "E57byAVW"

#remove 5 charters from column label and replace with E57by to allow merge with Nanodrop (alteration fixed later)
challenge$labels<-gsub("E57bx","",as.character(challenge$labels))
challenge$labels<-gsub("E57by","",as.character(challenge$labels))
challenge$labels = paste('E57by',challenge$labels, sep = '')
#merge Nanodrop data with main challenge 
challenge = left_join(challenge, Nanodrop)
##rename MICE LM0258 & LM0259 as E57bxAAA (the rest are E57byAAA)
challenge$labels[challenge$EH_ID == "LM0258" |challenge$EH_ID == "LM0259"] <- gsub("E57by","",as.character(challenge$labels))
challenge$labels[challenge$EH_ID == "LM0258"|challenge$EH_ID == "LM0259"] <- paste('E57bx',challenge$labels, sep = '')

rm(Nanodrop, fecweight,data)

####################################################
#OPG calculations
####################################################

challenge$oocyst_sq1 <- as.numeric(as.character(challenge$oocyst_sq1))
challenge$oocyst_sq2 <- as.numeric(as.character(challenge$oocyst_sq2)) 
challenge$oocyst_sq3 <- as.numeric(as.character(challenge$oocyst_sq3)) 
challenge$oocyst_sq4 <- as.numeric(as.character(challenge$oocyst_sq4)) 
challenge$fecweight_DNA <- as.numeric(as.character(challenge$fecweight_DNA)) 

challenge %>% rename(
  fecweight_flot = feces_weight,
  fecweight_DNA = fecweight_qPCR) -> challenge

###Use function from Alice Balard to calculate OPG
#volume of each large square (1mm x 1mm x 0.1mm height = 0.1mm3)
#0.1mm3 = 0.1µl (x10000 when converted into ml)
#https://www.emsdiasum.com/microscopy/technical/datasheet/68052-14.aspx <- explanation 
calculateOPG <- function(challenge){
  challenge$mean_Neubauer <- 
    (challenge$oocyst_sq1 + challenge$oocyst_sq2 + challenge$oocyst_sq3 + challenge$oocyst_sq4) / 4
  # NB! Limit of detection = 1 oocysts
  challenge$mean_Neubauer[challenge$oocyst_sq1 + challenge$oocyst_sq2 + challenge$oocyst_sq3 + challenge$oocyst_sq4 == 1] <- 0
  challenge$oocysts.per.tube <- challenge$mean_Neubauer * 10000 * challenge$dilution
  challenge$OPG <- challenge$oocysts.per.tube / challenge$fecweight_flot
  ## If we don't have the fecal weight BUT we counted in Neubauer chamber 0, then OPG = 0
  challenge$oocysts.per.tube[challenge$fecweight_flot == 0 & challenge$mean_Neubauer == 0] <- 0
  challenge$OPG[challenge$fecweight_flot == 0 & challenge$mean_Neubauer == 0] <- 0
  return(challenge)
}
challenge <- calculateOPG(challenge = challenge)
##Check for spaces
challenge$EH_ID <- gsub(pattern = " ", replacement = "", x = challenge$EH_ID)

###Estimate microbiota density (MD) from Contijoch et al. 2019
### MD = total DNA per sample (µg)/ mg of fresh feces
##considering 40µL of elution volume 
#Conc_DNA is in ng/µl, hence, x0.001 
challenge %>%
  mutate(Total_DNA = (challenge$Conc_DNA*40)*0.001) -> challenge ### Add a new variable that will contain total DNA extracted per sample in µg

challenge %>%
  mutate(Microbial_density = challenge$Total_DNA/(challenge$fecweight_DNA*1000)) -> challenge ### Total DNA extracted per sample in µg by feces weight in mg

challenge$labels<- as.vector(challenge$labels)
rownames(challenge) <- make.unique(challenge$labels)

##Transform to zero OPGs for DPI 1 and 2 and 3
#challenge$OPG[challenge$dpi==1] <- 0
#challenge$OPG[challenge$dpi==2] <- 0  
#challenge$OPG[challenge$dpi==3] <- 0

rm(calculateOPG)
#write.csv(challenge, "Output_Data/All_parameters_merged_challenge.csv", row.names = FALSE)


####################################################
#qPCR data calculations
####################################################



##Merge qPCR data 
#Ms. Fay Webster's function was adpated to merge all qPCR files into one
#setwd("Raw_Data/qPCR/qPCR_Secondary/")

#all individual qPCR.csv files form a list
#list_faeces <- as.list(list.files())

#to all a function to perform itself on the list, it is transformed into a vector
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

#rm(list_results, list_faeces,list_names, read_qPCR_file)

#write.csv(df_results, "qPCR_fecal_lab_E64_challenge_merged.csv", row.names=FALSE)
##moved to "\Eimeria_Quant2\Output_Data"
#rm(df_results)



















