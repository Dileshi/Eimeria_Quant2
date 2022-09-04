##Plot challenge infections 

library(tidyverse)
library(dplyr)
library(tidyr)

challenge <- read.csv("Output_Data/Quant2_challenge.csv")

##rename column 'feces_weight'  to 'fecweight_flot' to allow cohesion with Victor's data
challenge <- challenge %>% 
  rename(fecweight_flot = feces_weight)

##R identifies the entries under column ... as 'NULL'
#Hence, class changed to numeric
challenge$oocyst_sq1<-sapply(challenge$oocyst_sq1, as.numeric)
challenge$oocyst_sq2<-sapply(challenge$oocyst_sq2, as.numeric)
challenge$oocyst_sq3<-sapply(challenge$oocyst_sq3, as.numeric)
challenge$oocyst_sq4<-sapply(challenge$oocyst_sq4, as.numeric)
challenge$dilution<-sapply(challenge$dilution, as.numeric)


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
  challenge$OPG<- challenge$oocysts.per.tube / challenge$fecweight_flot
  ## If we don't have the fecal weight BUT we counted in Neubauer chamber 0, then OPG = 0
  challenge$oocysts.per.tube[challenge$fecweight_flot == 0 & challenge$mean_Neubauer == 0] <- 0
  challenge$OPG[challenge$fecweight_flot == 0 & challenge$mean_Neubauer == 0] <- 0
  return(challenge)
}

challenge <- calculateOPG(challenge = challenge)

##Check for spaces
challenge$EH_ID <- gsub(pattern = " ", replacement = "", x = challenge$EH_ID)

#NOT DONE ###############3###Estimate microbiota density (MD) from Contijoch et al. 2019
### MD = total DNA per sample (µg)/ mg of fresh feces
##considering 40µL of elution volume 
#Conc_DNA is in ng/µl, hence, x0.001 
#challenge %>%
#  mutate(Total_DNA = (challenge$Conc_DNA*40)*0.001) -> challenge ### Add a new variable that will contain total DNA extracted per sample in µg

#challenge %>%
 # mutate(Microbial_density = challenge$Total_DNA/(challenge$fecweight_DNA*1000)) -> challenge ### Total DNA extracted per sample in µg by feces weight in mg

#challenge$labels<- as.vector(challenge$labels)
#rownames(challenge) <- make.unique(challenge$labels)




