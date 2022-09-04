## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Quant2
## Script build for sample stage specific Eimeria quantification data set

library(tidyverse)
library(dplyr)
library(tidyr)
library(janitor)
library(readr)

###  First, data must be merged; use code from Ms. Fay Webster that was also applied
###   for primary infection data (script:1_Data_Prep)
setwd("Raw_Data/qPCR/qPCR_RNA/")
#all individual qPCR.csv files form a list
list_faeces <- as.list(list.files())
#to all a function to perform itself on the list, it is transformed into a vector
list_names <- as.vector(unlist(list_faeces))
#function set by Fay to deleted rows 1 to 24 and set a new column that includes the file
#name i.e. date of qPCR experiment
read_qPCR_file <- function(x) {
  df1 = read_csv(x)
  filename <- colnames(df1[2])
  df1 <- df1 %>%
    filter(!row_number() %in% c(1:24))
  colnames(df1) <- df1[1, ]
  df1 <- df1 %>% filter(!row_number() %in% 1)
  df1 <- df1 %>% dplyr::mutate(plate = filename)
  }
list_results <- lapply(list_names, read_qPCR_file)
#bind/merge all indivual csv files into one
df_results <- Reduce(rbind, list_results)
#removes duplicates
df_results <- unique(df_results)  
rm(list_results, list_faeces,list_names, read_qPCR_file)
setwd("/Users/vinuri/Documents/GitHub/Eimeria_Quant2")
write.csv(df_results, "Output_Data/qPCR_fecal_lab_RNA_merged.csv", row.names=FALSE)
rm(df_results)


###Merge qPCR RNA-data to dataframe containing sample data (dpi ect..)

##Load RNA-merged qPCR data
RNA.data <- read.csv("Output_Data/qPCR_fecal_lab_RNA_merge.csv")
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
RNA.data.Sporo <- RNA.data.Sporo[RNA.data.Sporo$Infection == "Positive", ]

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
RNA.data.Mero <- RNA.data.Mero[RNA.data.Mero$Infection == "Positive", ]


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
RNA.data.Oocy <- RNA.data.Oocy[RNA.data.Oocy$Infection == "Positive", ]


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
RNA.data.Micro <- RNA.data.Micro[RNA.data.Micro$Infection == "Positive", ]


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
RNA.data.Macro <- RNA.data.Macro[RNA.data.Macro$Infection == "Positive", ]


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
RNA.data.Ear.Ooc <- RNA.data.Ear.Ooc[RNA.data.Ear.Ooc$Infection == "Positive", ]



#merge all development stage data frames 
RNA.data.all <- bind_rows(RNA.data.Mero,RNA.data.Oocy,RNA.data.Sporo,RNA.data.Micro,RNA.data.Macro,
                          RNA.data.Ear.Ooc) 
rm(RNA.data.Mero,RNA.data.Oocy,RNA.data.Sporo,RNA.data.Micro,RNA.data.Macro,RNA.data.Ear.Ooc) 



##################PLOT
#RNA.data.all$dpi<-sapply(RNA.data.all$dpi, as.factor)

RNA.data.all%>%
  filter(dpi%in%c("2","3","4", "5", "7", "8", "10"))%>%
  dplyr::select(EH_ID, dpi, Cq.Mean, Cq.SD, Target, labels)%>%
  dplyr::arrange(labels)%>%
  dplyr::arrange(dpi)%>% 
  ggplot(aes(x= factor (dpi,levels=1:10), y= Cq.Mean, color=factor(dpi)))+
  scale_x_discrete('dpi', breaks=factor(1:10))+
  scale_color_manual(values=c("#9b5fe0", "#16a4d8", "#60dbe8","#8bd346","#efdf48","#f9a52c","#d64e12"))+
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Target, fill= factor(dpi)))+
  geom_errorbar(aes(ymin = Cq.Mean-Cq.SD, ymax=  Cq.Mean+Cq.SD), width=0.2, color="gray") -> A


##Figure X: Developemnt stages 
pdf(file = "Figures/RNA_Test2.pdf")
print(A)
dev.off()



####Jitter Cq.SD  whiskers from point (not align)
###®emove factor bar/legend









