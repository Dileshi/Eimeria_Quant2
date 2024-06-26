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
rm(sample.data)
#load qPCR data challenge infection and clean
DNA.data <- read.csv("Output_Data/qPCR_fecal_lab_E64_challenge_merged.csv")
DNA.data <- rename (DNA.data, Cq_mean = Cq.Mean)

DNA.data.brief <- DNA.data[c(2,4,13,22,26)]
#remove old replicates (use new plate 11.10 plate)
DNA.data<-DNA.data[-c(271:273,43:48,280:282,55:57,286:288,64:66,223:225,13:15,70:72,187:189,88:93,199:201,130:132,235:237), ]
#remove NTC and NA
DNA.data <- DNA.data[DNA.data$Sample != "NTC", ] 
DNA.data <- DNA.data[DNA.data$Sample != "NTC1", ] 
DNA.data <- DNA.data[DNA.data$Sample != "NTC2", ]
DNA.data <- DNA.data[DNA.data$Sample != "NTC3", ]
DNA.data <- DNA.data[!is.na(DNA.data$Sample),]
rm(DNA.data.brief)

#double check for dublicates
Freq_list <- data.frame(table(DNA.data$Sample))
Freq_list = Freq_list[Freq_list$Freq > 3,]
rm(Freq_list)

#create new column in dataframe in challenge to allow merge with DNA.data dataframe
#challenge dataframe new column named sample will contain labels without E57x or E57y to allow merge
    #challenge <- read.csv("Output_Data/All_parameters_merged_challenge.csv")
#remove "E57bx" & "E57by" in DNA.data dataframe
DNA.data$Sample<-gsub("E57bx","",as.character(DNA.data$Sample))
DNA.data$Sample<-gsub("E57by","",as.character(DNA.data$Sample))

#subset required info
challenge %>%
  select(labels, EH_ID, dpi,Conc_DNA,fecweight_DNA,OPG,relative_weight,weight_dpi0,weight) -> challenge.DNA
#introduce new column named Sample containing 3-letter code
challenge.DNA$Sample = challenge.DNA$labels
challenge.DNA$Sample<-gsub("E57bx","",as.character(challenge.DNA$Sample))
challenge.DNA$Sample<-gsub("E57by","",as.character(challenge.DNA$Sample))
#merge qPCR data and sample data 
challenge.DNA = left_join(challenge.DNA,DNA.data)
rm(DNA.data,challenge)

#rename: cohesive with other scripts
challenge.DNA <- rename (challenge.DNA, c(Ct = Cq, Tm = Tm1))
data.std <- rename (data.std, c(Tm = Tm1))

challenge.DNA$fecweight_DNA<- as.numeric(challenge.DNA$fecweight_DNA)

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
  dplyr::select(labels, Genome_copies, Tm, Infection,EH_ID,dpi,Cq_mean, Conc_DNA,fecweight_DNA,OPG,relative_weight,weight_dpi0,weight)%>%
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

#set NA to zero 
challenge.DNA$Genome_copies_gFaeces[challenge.DNA$labels=="E57byOSY"] <- 0
challenge.DNA$Genome_copies_gFaeces[challenge.DNA$labels=="E57bxOPZ"] <- 0  










