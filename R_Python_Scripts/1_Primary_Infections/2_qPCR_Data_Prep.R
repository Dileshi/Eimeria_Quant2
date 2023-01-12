## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Quant2
## Script adapted from Victor's 2_qPCR_data_preparation


##Plot standard curve 
library(tidyverse)
library(dplyr)
library(tidyr)
library(janitor)
library(readr)
library(ggtext)
library(lmtest)
library(ggpubr)
library(rcompanion)
library(gridExtra)
library(rstatix)
library(lmtest)
library(ggtext) 

##load standard curve data (lab_fecal) and clean 
#qPCR executed by Josi 
data.std<-read.csv("Raw_Data/Standard_curve.csv")

#delete unwanted rows
data.std <- data.std %>%
  filter(!row_number() %in% c(1:24))
  colnames(data.std) <- data.std[1, ]
  data.std <- data.std %>% filter(!row_number() %in% 1)
  data.std <- data.std %>% filter(!row_number() %in% c(1:66))
  data.std <- data.std %>% filter(!row_number() %in% c(13:15))
  data.std <- data.std %>% filter(!row_number() %in% c(22:27))
  data.std <- data.std %>% filter(!row_number() %in% c(4:6))

#change column name 'Sample' to 'labels' of qPCR file
data.std <- data.std %>% 
  rename(labels = Sample)
data.std <- data.std %>% rename(Cq_mean = 'Cq Mean')

data.std <- data.std %>%
    dplyr::mutate(labels = case_when(
      labels == "v10_2" ~ "V10_2",
      labels == "v10_5" ~ "V10_5",
      TRUE ~ labels
    ))
data.std <- select (data.std, labels, Cq_mean, Tm1,Tm2,Tm3,Tm4)

y = str_extract(data.std$labels, "\\d$")  
x =as.numeric(y)
data.std%>%
  dplyr::mutate(Oocyst_count= 10**x)-> data.std
data.std$Cq_mean = as.numeric(data.std$Cq_mean)
                
data.std%>%
  dplyr::mutate(Genome_copies= Oocyst_count*8)-> data.std  ###Not multiplying by 8
###unknown if sporulated or unsporulated 

data.std%>%
  ggplot(aes(x = Oocyst_count, y = Cq_mean)) +
  geom_point(color = "tomato", size=3) +
  ylim(10, 35) +
  scale_x_log10("log 10 *Eimeria* genome copies", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_smooth(method="lm", se=FALSE, color="gray")+
  stat_cor(label.x =3, label.y = 32, size=3.5, 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+ 
  stat_regline_equation(label.x =3, label.y = 33)+
  theme_bw() +
  labs(title = "Standard Curve of *Eimeria falciformis*",
     y = "Cycle Threshold") + 
  theme(text = element_text(size=11)) +
  theme(axis.title.x = ggtext::element_markdown(color = "tomato", size = rel(1))) + 
  theme(axis.title = ggtext::element_markdown(color = "tomato", size = rel(1))) + 
  theme(plot.title = element_markdown(size = 12, hjust=0.5))

rm(x,y)

###################################################################################################

##Load data
if(!exists("sample.data")){
  source("R_Python_Scripts/1_Primary_Infections/1_Data_Prep.R")
}

##Infection experiment: load qPCR data (lab_fecal) and clean 
data.inf.exp<-read.csv("Output_Data/qPCR_fecal_lab_E88_primary_merged.csv")

data.inf.exp <- rename (data.inf.exp, Cq_mean = Cq.Mean)  
lm.data.std<- lm(log10(Genome_copies)~Cq_mean, data.std, na.action = na.exclude)
data.inf.exp$Genome_copies<- 10^predict(lm.data.std, data.inf.exp)  
rm(lm.data.std,data.std)

### adapt all data in files to allow merge into challenge_infection file:
# this includes changing column names of files to match that of challenge_infection file
# change column name 'Sample' to 'labels' of qPCR file
data.inf.exp <- data.inf.exp %>% 
  rename(labels = Sample)

#remove NA rows 
data.inf.exp <- data.inf.exp %>% drop_na(labels) 

#entries under column require the characters 'E57a' to allow cohesive merge
data.inf.exp$labels[223:399] <- paste0('E57a', data.inf.exp$labels[223:399])
data.inf.exp$labels[880:936] <- paste0('E57a', data.inf.exp$labels[880:936])
data.inf.exp$labels[1225:1263] <- paste0('E57a', data.inf.exp$labels[1225:1263])
data.inf.exp <- data.inf.exp %>%
  dplyr::mutate(labels = case_when(
    labels == "E57INR" ~ "E57aINR",
    labels == "E57CDE" ~ "E57aCDE", 
    labels == "E57CEW" ~ "E57aCEW", 
    labels == "E57EFU" ~ "E57aEFU",
    labels == "E57NTC" ~ "E57aNTC",
    labels == "E57aNTC1" ~ "E57aNTC",
    labels == "E57aNTC2" ~ "E57aNTC",
    labels == "E57aBHN2" ~ "E57aBHN",
    labels == "E57aFLV2" ~ "E57aFLV",
    labels == "E57aIJQ2" ~ "E57aIJQ",
    labels == "E57aRWY2" ~ "E57aRWY",
    labels == "NTC" ~ "E57aNTC",
    labels == "NTC1" ~ "E57aNTC",
    labels == "NTC2" ~ "E57aNTC",
    TRUE ~ labels
  ))

#unsuccessful qPCR replicates are deselected from data frame
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aABD" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aABL" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aACN" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aACP" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aACP" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aACY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aAFR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aAJT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aAJT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aALU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aATW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBCE" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBCE" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBCR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBCR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBDY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBHN" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBHN" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBLN" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBLX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBLX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBOR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBOR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aBTZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCEG" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_a.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCEG" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_11042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCEG" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCFG" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCIW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCIW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCLT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCLT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aCNS" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCQT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCQT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCUX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCUX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aDMV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aDTY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aDTY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aAKO" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aBGU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aBGU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCFW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCJQ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aCSX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aDEK" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEIP" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEIP" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEIY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEKO" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEKO" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEQU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEQU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aESV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aESV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aETU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aETU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEKR" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEMX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEUW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aEMX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFJQ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFKM" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFKW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFRW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aFRW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGHL" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGHU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGHU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGNZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aGRV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHIS" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHJT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHOX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHOX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aHQZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIJQ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIKZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIOX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIOY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIPT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_11042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aIPY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aJOS" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aKMX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aKOX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aKUV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aKUV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLMW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLNU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLNU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLPU" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLSZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aLXZ" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_26042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aMQT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_a.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aMRT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_11042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aMRT" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- filter(data.inf.exp, labels != "E57aNTC")
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aOPY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_01052022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aOPY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aQVY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aQVY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRTV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRTV" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_29042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aRWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aSUW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_25042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aTUY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aTUY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aUWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aUWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aUWY" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_b.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aVWX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_17042022.eds"),]
data.inf.exp <- data.inf.exp [!(data.inf.exp$labels == "E57aVWX" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_22042022.eds"),]

##E57aRSW is a peculiar case, check T1,T2,T3 values of all replicates to recall 
data.inf.exp <- data.inf.exp[!(data.inf.exp$labels == "E57aRSW" & data.inf.exp$plate == "qPCR_eimeria_lab_fecal_DB_30042022_c.eds"),]

#to unsure the absence of qPCR replicates, the following code should produce no results
#the number of qPCR values for each label (E57aXXX) should be at a frequency of 3
#code displays labels containing more than 3 replicates 
Freq_list <- data.frame(table(data.inf.exp$labels))
Freq_list = Freq_list[Freq_list$Freq > 3,]
rm(Freq_list)

##apply Victors script logically
data.inf.exp <- rename (data.inf.exp, c(Ct = Cq, Tm = Tm1))


######### Infection experiment data############
## Define real positive and negatives based on Tm 

### Estimate mean Eimeria Tm for positive controls
data.inf.exp%>%
  dplyr::select(Task,labels,Tm)%>%
  dplyr::filter(Task%in%c("Pos_Ctrl"))%>%
  dplyr::filter(complete.cases(.))%>%
  dplyr::mutate(Tm = as.numeric(Tm))%>%
  #dplyr::slice(1,4,7,11,13)%>%
  dplyr::summarise(mean = mean(Tm), sd= sd(Tm), n = n())%>%
  dplyr::mutate(se = sd / sqrt(n),
                lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
                upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)%>%
  dplyr::mutate(upper.ran = mean + (2*sd), 
                lower.ran = mean - (2*sd))

##Positive controls have a Tm at:
#mean        sd n  se lower.ci upper.ci upper.ran lower.ran
#74.1 0.8944272 5 0.4 72.98942 75.21058  75.88885  72.31115

#there is a second Tm higher than 80°C from an unspecific product

##The mean Tm is at 74.1°C and sd of 0.89
##Make all the things below/above Tm mean +/- 2sd negative 
##Check which samples have the three triplicates "positive"
data.inf.exp %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm)  ~ "Negative",
                                      (Tm >= 76 | Tm <= 72.2)  ~ "Negative",
                                      (Tm >=72.3 & Tm <= 75.9) ~ "Positive"))%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise(Count.Tm= sum(Infection=="Positive"))%>%
  dplyr::mutate(Infection = case_when(Count.Tm == 3 ~ "Positive",
                                      Count.Tm != 3 ~ "Negative"))%>%
  dplyr::left_join(data.inf.exp, by= "labels")-> data.inf.exp 


##Eliminate an unprocessed sample and controls
data.inf.exp%>%
  dplyr::mutate(Genome_copies = case_when((Infection=="Negative")~ 0,
                                          TRUE~ Genome_copies))%>% ## Make Negative zero
  dplyr::select(labels, Genome_copies, Tm, Infection,Cq_mean,)%>%
  dplyr::filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl"))%>% ## Replace NAs in real negative samples to 0 
  dplyr::mutate(Genome_copies= replace_na(Genome_copies, 0))%>%
  dplyr::mutate(Tm= replace_na(Tm, 0))%>%
  ##Get unique labels from qPCR data
  dplyr::distinct(labels, .keep_all = TRUE)-> data.inf.exp

####============

#remove column infection (shows:primary): could clash with infection +ive or _ive
sample.data <- sample.data[ , ! names(sample.data) %in% "infection"]

###Join all the data in the same dataframe
sdt<- left_join(sample.data, data.inf.exp, by="labels") ## Add qPCR data
### adjustment  
sdt$dpi<- as.factor(sdt$dpi)

rm(data.inf.exp, sample.data)

sdt%>%
  dplyr::mutate(Genome_copies_ngDNA= Genome_copies/50, ## copies by ng of fecal DNA considering 1uL from 50 ng/uL DNA
                DNA_sample= Conc_DNA*40, ## Estimate total gDNA of sample considering 40uL of elution buffer
                DNA_g_feces= DNA_sample/fecweight_DNA,
                ## Transform it to ng fecal DNA by g of faeces
                Genome_copies_gFaeces= Genome_copies_ngDNA*DNA_g_feces) -> sdt ## Estimate genome copies by g of faeces

##Transform to zero OPGs for DPI 1 and 2 and 3
sdt$OPG[sdt$dpi==1] <- 0
sdt$OPG[sdt$dpi==2] <- 0  
sdt$OPG[sdt$dpi==3] <- 0

#write.csv(sdt,"Output_Data/All_Parameters_Fecal_Lab_E88_primary.csv")



  
