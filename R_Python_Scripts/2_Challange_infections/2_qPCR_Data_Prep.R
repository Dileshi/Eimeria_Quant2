##Run scripts on Root of Repo (Eimeria_Quant2)

##Plot challenge infections 
library(tidyverse)
library(dplyr)
library(tidyr)

##BEGIN WITH qPCR data -> plot DNA using primary infection standard curve
##load standard curve (Josi)
if(!exists("data.std")){
  source("R_Python_Scripts/1_Primary_Infections/2_qPCR_Data_Prep.R")
  }

#load qPCR data challenge infection and clean
DNA.data <- read.csv("Output_Data/qPCR_fecal_lab_challenge_merged.csv")
DNA.data <- rename (DNA.data, Cq_mean = Cq.Mean)

#create new column in dataframe in challenge to allow merge with DNA.data dataframe
#challenge dataframe new column named sample will contain labels without E57x or E57y to allow merge
challenge <- read.csv("Output_Data/All_parameters_merged_challenge.csv")
#subset required info
challenge %>%
  select(labels, EH_ID, dpi,Conc_DNA,fecweight_DNA,OPG,relative_weight) -> challenge.DNA
#introduce new column named Sample containing 3-letter code
challenge.DNA$Sample = challenge.DNA$labels
challenge.DNA$Sample<-gsub("E57bx","",as.character(challenge.DNA$Sample))
challenge.DNA$Sample<-gsub("E57by","",as.character(challenge.DNA$Sample))
#merge qPCR data and sample data 
challenge.DNA = left_join(DNA.data,challenge.DNA)
rm(DNA.data,challenge)

#rename: cohesive with other scripts
challenge.DNA <- rename (challenge.DNA, c(Ct = Cq, Tm = Tm1))
data.std <- rename (data.std, c(Tm = Tm1))

##PREDICT GENOME COPIES USING STANDARD
lm.data.std<- lm(log10(Genome_copies)~Cq_mean, data.std, na.action = na.exclude)
challenge.DNA$Genome_copies<- 10^predict(lm.data.std, challenge.DNA)  
rm(lm.data.std,data.std)

######### Challenge experiment data############
## Define real positive and negatives based on Tm 
### Estimate mean Eimeria Tm for positive controls
challenge.DNA%>%
  dplyr::select(Task,labels,Tm)%>%
  dplyr::filter(Task%in%c("Pos_Ctrl"))%>%
  dplyr::filter(complete.cases(.))%>%
  dplyr::mutate(Tm = as.numeric(Tm))%>%
  dplyr::summarise(mean = mean(Tm), sd= sd(Tm), n = n())%>%
  dplyr::mutate(se = sd / sqrt(n),
                lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
                upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)%>%
  dplyr::mutate(upper.ran = mean + (2*sd), 
                lower.ran = mean - (2*sd))

##Positive controls have a Tm at:
#mean        sd n  se lower.ci upper.ci upper.ran lower.ran
#74.1 0.8944272 5 0.4 72.98942 75.21058  75.88885  72.31115

##The mean Tm is at 74.1°C and sd of 0.89
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
challenge.DNA %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm)  ~ "Negative",
                                      (Tm >= 76 | Tm <= 72.2)  ~ "Negative",
                                      (Tm >=72.3 & Tm <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm == 3 ~ "Positive",
                                      Count.Tm != 3 ~ "Negative"))%>%
  dplyr::left_join(challenge.DNA, by= "labels")-> challenge.DNA

##Eliminate an unprocessed sample and controls
challenge.DNA%>%
  dplyr::mutate(Genome_copies = case_when((Infection=="Negative")~ 0,
                                          TRUE~ Genome_copies))%>% ## Make Negative zero
  dplyr::select(labels, Genome_copies, Tm, Infection,EH_ID,dpi,Cq_mean, Conc_DNA,fecweight_DNA,OPG,relative_weight)%>%
  dplyr::filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl"))%>% ## Replace NAs in real negative samples to 0 
  dplyr::mutate(Genome_copies= replace_na(Genome_copies, 0))%>%
  dplyr::mutate(Tm= replace_na(Tm, 0))%>%
  ##Get unique labels from qPCR data
  dplyr::distinct(labels, .keep_all = TRUE)-> challenge.DNA

challenge.DNA%>%
  dplyr::mutate(Genome_copies_ngDNA= Genome_copies/50, ## copies by ng of fecal DNA considering 1uL from 50 ng/uL DNA
                DNA_sample= Conc_DNA*40, ## Estimate total gDNA of sample considering 40uL of elution buffer
                DNA_g_feces= DNA_sample/fecweight_DNA,
                ## Transform it to ng fecal DNA by g of faeces
                Genome_copies_gFaeces= Genome_copies_ngDNA*DNA_g_feces) -> challenge.DNA ## Estimate genome copies by g of faeces

#####NOTE: Out 0f 73 samples 68 are calculated -> check 5 missing samples 











