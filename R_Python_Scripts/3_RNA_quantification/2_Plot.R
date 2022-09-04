library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(reshape2)

##Load RNA-merged qPCR data
RNA.data <- read.csv("Output_Data/Sample_RNA_DEL.csv") 

##As we are only focusing on 

sample.data <- read.csv("Raw_Data/Quant2_E57.csv")
## Filter out the E64 (E. ferrisi) mice in Quant2_E57 file
sample.data <- sample.data %>% 
  dplyr::filter(!primary_infection == "E64")
##rename column 'feces_weight' of Quant2_E57 to 'fecweight_flot'
sample.data <- sample.data %>% 
  rename(fecweight_flot = feces_weight)

#change column name 'Sample' to 'labels' of qPCR file
RNA.data <- RNA.data %>% 
  rename(labels = Sample)
#keep essential data/colunms from sample.data dataframe
sample.data <- sample.data %>% select(EH_ID, experiment, labels, dpi)
#merge
RNA.data <- merge(RNA.data, sample.data, by = 'labels', all.x=TRUE)


#if Cq.Mean is NA/(empty:no amplification) add 40.  
RNA.data$Cq.Mean[is.na(RNA.data$Cq.Mean)] <- 40


##The mean Tm for Sporozoites is at XXX°C and sd of XXX
#filter out sporoziotes
RNA.data.Sporo <- RNA.data[RNA.data$Target == "Sporo", ]
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
RNA.data.Sporo %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm1)  ~ "Negative",
                                      (Tm1 >= 76 | Tm1 <= 72.2)  ~ "Negative",
                                      (Tm1 >=72.3 & Tm1 <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm1= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm1 == 3 ~ "Positive",
                                      Count.Tm1 != 3 ~ "Negative"))%>%
  dplyr::left_join(RNA.data.Sporo, by= "labels")-> RNA.data.Sporo 
#keep Infection == postive samples
#RNA.data.Sporo <- RNA.data.Sporo[RNA.data.Sporo$Infection == "Positive", ]

#repeat for all other developemnt stages (Merozoites=Mero, Oocyst=Oocy, 
#Early oocyts=Ear.Ooc, Microgamtes=Micro, Macrogamtes=Macro)

##The mean Tm for Merozoites is at XXX°C and sd of XXX
#filter out Merzoites
RNA.data.Mero <- RNA.data[RNA.data$Target == "Mero", ]
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
RNA.data.Mero %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm1)  ~ "Negative",
                                      (Tm1 >= 76 | Tm1 <= 72.2)  ~ "Negative",
                                      (Tm1 >=72.3 & Tm1 <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm1= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm1 == 3 ~ "Positive",
                                      Count.Tm1 != 3 ~ "Negative"))%>%
  dplyr::left_join(RNA.data.Mero, by= "labels")-> RNA.data.Mero 
#keep Infection == postive samples
#RNA.data.Mero <- RNA.data.Mero[RNA.data.Mero$Infection == "Positive", ]


##The mean Tm for Oocyst is at XXX°C and sd of XXX
#filter out Oocyst
RNA.data.Oocy <- RNA.data[RNA.data$Target == "Oocy", ]
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
RNA.data.Oocy %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm1)  ~ "Negative",
                                      (Tm1 >= 76 | Tm1 <= 72.2)  ~ "Negative",
                                      (Tm1 >=72.3 & Tm1 <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm1= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm1 == 3 ~ "Positive",
                                      Count.Tm1 != 3 ~ "Negative"))%>%
  dplyr::left_join(RNA.data.Oocy, by= "labels")-> RNA.data.Oocy 
#keep Infection == postive samples
#RNA.data.Oocy <- RNA.data.Oocy[RNA.data.Oocy$Infection == "Positive", ]


##The mean Tm for Microgametes is at XXX°C and sd of XXX
#filter out Microgametes
RNA.data.Micro <- RNA.data[RNA.data$Target == "Micro", ]
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
RNA.data.Micro %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm1)  ~ "Negative",
                                      (Tm1 >= 76 | Tm1 <= 72.2)  ~ "Negative",
                                      (Tm1 >=72.3 & Tm1 <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm1= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm1 == 3 ~ "Positive",
                                      Count.Tm1 != 3 ~ "Negative"))%>%
  dplyr::left_join(RNA.data.Micro, by= "labels")-> RNA.data.Micro 
#keep Infection == postive samples
#RNA.data.Micro <- RNA.data.Micro[RNA.data.Micro$Infection == "Positive", ]


##The mean Tm for Macrogametes is at XXX°C and sd of XXX
#filter out Macrogametes
RNA.data.Macro <- RNA.data[RNA.data$Target == "Macro", ]
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
RNA.data.Macro %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm1)  ~ "Negative",
                                      (Tm1 >= 76 | Tm1 <= 72.2)  ~ "Negative",
                                      (Tm1 >=72.3 & Tm1 <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm1= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm1 == 3 ~ "Positive",
                                      Count.Tm1 != 3 ~ "Negative"))%>%
  dplyr::left_join(RNA.data.Macro, by= "labels")-> RNA.data.Macro 
#keep Infection == postive samples
#RNA.data.Macro <- RNA.data.Macro[RNA.data.Macro$Infection == "Positive", ]


##The mean Tm for Early Oocyst is at XXX°C and sd of XXX
#filter out Early Oocyst
RNA.data.Ear.Ooc <- RNA.data[RNA.data$Target == "Ear-Ooc", ]
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
RNA.data.Ear.Ooc %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm1)  ~ "Negative",
                                      (Tm1 >= 76 | Tm1 <= 72.2)  ~ "Negative",
                                      (Tm1 >=72.3 & Tm1 <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm1= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm1 == 3 ~ "Positive",
                                      Count.Tm1 != 3 ~ "Negative"))%>%
  dplyr::left_join(RNA.data.Ear.Ooc, by= "labels")-> RNA.data.Ear.Ooc 
#keep Infection == postive samples
#RNA.data.Ear.Ooc <- RNA.data.Ear.Ooc[RNA.data.Ear.Ooc$Infection == "Positive", ]



#merge all development stage data frames 
RNA.data.all <- bind_rows(RNA.data.Mero,RNA.data.Oocy,RNA.data.Sporo,RNA.data.Micro,RNA.data.Macro,
                          RNA.data.Ear.Ooc) 
rm(RNA.data.Mero,RNA.data.Oocy,RNA.data.Sporo,RNA.data.Micro,RNA.data.Macro,RNA.data.Ear.Ooc) 



#subset containing single sample (remove dublicates)
RNA.single <- distinct(RNA.data.all, labels,Target, .keep_all= TRUE)
#display just important rows
RNA.single <- RNA.single[c("labels","Infection","Cq.Mean", "Target", "EH_ID","dpi","Cq.SD")]



#####How plot development stages 
#####Method 1: add negative sign in front of subtraction 
RNA.single %>% 
  group_by(labels)%>%
  summarise(Cq.max = max(Cq.Mean, na.rm=TRUE))%>%
  dplyr::left_join(RNA.single, by= "labels")-> RNA.single
RNA.single %>% 
  group_by(labels)%>%
  dplyr::mutate(Cq.Ratio = -1(Cq.Mean/Cq.max)) -> RNA.single



#####Method 2: calculate stage-scores by substracting lowest for each stage
RNA.single %>% 
  dcast(labels ~ Target, value.var = "Cq.Mean") %>% 
  dplyr::mutate(Mer_Gam = Mero-Micro)%>%
  dplyr::mutate(Mer_Ooc = Mero-Oocy)%>%
  dplyr::mutate(Ooc_Mer = Oocy-Mero)%>%
  dplyr::mutate(Ooc_Gam = Oocy-Micro)%>%
  dplyr::mutate(Gam_Mer = Micro-Mero)%>%
  dplyr::mutate(Gam_Ooc = Micro-Oocy)%>%
  dplyr::left_join(RNA.single, by= "labels")-> RNA.single
#remove addiotional rows
drop <- c("Macro","Mero","Oocy","Sporo","Micro","Ear-Ooc")
RNA.single = RNA.single[,!(names(RNA.single) %in% drop)]
rm(drop)




#####
#####RNA.data.all$dpi<-sapply(RNA.data.all$dpi, as.numeric)
###plot plot
RNA.single%>%
  filter(dpi%in%c("0","1","2","3","4","5","6", "7", "8", "9","10","11"))%>%
  dplyr::select(EH_ID, dpi, Cq.Mean, Cq.SD, Target, labels, Cq.Ratio)%>%
  dplyr::arrange(labels)%>%
  dplyr::arrange(dpi)%>% 
  ggplot(aes(x=dpi, y= Cq.Ratio, color=dpi))+
  scale_x_continuous(labels = 0:11, breaks = 0:11)+
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Target))-> A

##Figure X: Developemnt stages 
pdf(file = "Figures/RNA_Test_Ratio_3.pdf")
print(A)
dev.off()

#scale_color_manual(values=c("#9b5fe0","#A4E05F","#64E05F","#8bd346","#16a4d8",
                            #"#1643D8", "#16D8AB","#60dbe8","","#efdf48","#EF48AB","#f9a52c","#d64e12"))
  







