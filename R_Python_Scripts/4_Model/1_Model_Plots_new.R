##Plot 3 
##Distribution plot of PRIMARY(E644,E88) & SECONDARY(E64)
#Compares effect of being primary/secondary for DNA, Weight, Oocyst

library(tidyverse)
library(ggeffects)

########IMPORT E88 PRIMARY DATA (E88_Prim)
if(!exists("sample.data")){
  source("R_Python_Scripts/1_Primary_Infections/1_Data_Prep.R")
}
if(!exists("sdt")){
  source("R_Python_Scripts/1_Primary_Infections/2_qPCR_Data_Prep.R")
}
##select the completed 18 mice of sdt dataframe 

E88_Prim <- sdt[sdt$EH_ID %in% c("LM0179","LM0180","LM0181", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229",
                                 "LM0238", "LM0191","LM0240","LM0244","LM0246","LM0247","LM0248","LM0194","LM0190","LM0199",
                                 "LM0197","LM0198","LM0188","LM0254","LM0255","LM0192"),]
E88_Prim = E88_Prim[,c("EH_ID","labels","dpi","weight","weight_dpi0","OPG","Tm","Infection","Genome_copies_ngDNA","Genome_copies_gFaeces","Cq_mean")]

########IMPORT E64 PRIMARY DATA (E64_Prim)
E64_Prim <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Quant/master/data/Eimeria_quantification_Inf_exp_data.csv")
E64_Prim <- E64_Prim %>% rename(Cq_mean = 'Cq.Mean')
E64_Prim <- E64_Prim[c(5,8)]
E64_Prim %>%   dplyr::distinct(Sample, .keep_all = TRUE)-> E64_Prim
E64_Prim <- E64_Prim %>% rename(labels = 'Sample')
E64_Prim_OPG <- read.csv("Output_Data/E64_Primary_Victor_all.csv") #OPG, Genome copy calulations + other data
E64_Prim<- left_join(E64_Prim_OPG, E64_Prim, by="labels") ## Add qPCR data
E64_Prim <- E64_Prim[c(1:6,9,23,55:57,60,61)]
E64_Prim <- E64_Prim[-c(6,7)]
rm(E64_Prim_OPG)
E64_Prim$labels<-paste0('E57a', E64_Prim$labels)
E64_Prim <- E64_Prim %>% rename(Tm = 'Tm_mean')

########IMPORT E64:E64 Secondary DATA (E64_Secon)
if(!exists("DNA.data")){
  source("R_Python_Scripts/2_Challange_infections/2_qPCR_Data_Prep.R")
}
if(!exists("challenge.DNA")){
  source("R_Python_Scripts/2_Challange_infections/2_qPCR_Data_Prep.R")
}
E64_E64 <- challenge.DNA
E64_E64 <- E64_E64[-c(2,8,9,11,15,16)]
rm(challenge.DNA)

####combine E64 primary and E64 challenge into one data frame 
##first add column primary
E64_Prim$Variable <- paste0('E64_Prim')
##first add column primary
E88_Prim$Variable <- paste0('E88_Prim')
##first add column secondary
E64_E64$Variable <- paste0('E64_Second')
#merge E64_Primary with E64:E64 secondary
E64_total <- rbind(E64_Prim,E64_E64)
E64_total <- rbind(E64_total,E88_Prim)

##contaminated 8 mice
E64_total%>%
  filter(dpi==0 & Variable== 'E88_Prim' & Genome_copies_gFaeces>0) -> contaminant
#make vector with EH_ID
EH_ID_Contaminant <- contaminant$EH_ID 

E64_total%>%
  filter(E64_total$EH_ID %in% EH_ID_Contaminant)%>%
  mutate(Variable = replace(Variable, Variable == 'E88_Prim', 'E88_Second'))->E64_total_conta

E64_total %>%
  filter(!E64_total$EH_ID %in% EH_ID_Contaminant) -> E64_total2
rbind(E64_total2,E64_total_conta) -> E64_total
             
rm(E64_E64,E64_Prim,E88_Prim,contaminant,EH_ID_Contaminant,E64_total2,E64_total_conta)



#########PLOT Distrubution plot for E64_Prim, E88_Prim, E64_Secon 
## filter(E64_total, dpi%in%c(3, 4, 5, 6, 7, 8)) %>%
E64_total$dpi <- as.factor(E64_total$dpi)
E64_total$dpi <-fct_inseq(E64_total$dpi)


lm(log10(Genome_copies_gFaeces+1)~dpi+Variable, data=E64_total) -> E64_Model 
##Distribution plot
#Genome
ggplot(E64_total,aes(y=log10(Genome_copies_gFaeces+1),x=Variable, color=dpi))+
  geom_violin()+
  geom_jitter(position=position_jitter(seed = 1, width = 0.2))+
  geom_smooth(method="lm")
#OPG
ggplot(E64_total,aes(y=log10(OPG+1),x=Variable, color=dpi))+
  geom_violin()+
  geom_jitter()+
  geom_smooth(method="lm")

summary(E64_Model)
##Distribution plot separated out into dpi
#Genome
ggplot(E64_total,aes(y=log10(Genome_copies_gFaeces+1),x=Variable, color=dpi))+
  geom_violin()+
  geom_jitter()+
  facet_wrap(~dpi)+
  geom_smooth(method="lm")
#OPG
ggplot(E64_total,aes(y=log10(OPG+1),x=Variable, color=dpi))+
  geom_violin()+
  geom_jitter()+
  facet_wrap(~dpi)+
  geom_smooth(method="lm")

################################################ 
#Residual plot

sdt <- E64_total

sdt%>%
  filter(!(OPG== 0))%>%
  filter(!(Genome_copies_gFaeces== 0))-> sdt.nozero

##Model 1: Genome copies/g faeces modeled by OPG
DNAbyOPG <- lm(log10(Genome_copies_gFaeces)~log10(OPG),
               data = sdt.nozero, na.action = na.exclude)
summary(DNAbyOPG)

### Plot and extract estimates
require(sjPlot)

sdt.nozero$predictedM1 <- predict(DNAbyOPG)   # Save the predicted values
sdt.nozero$residualsM1 <- residuals(DNAbyOPG) # Save the residual values

##Assign the colors for dpi and keep consistency with previous plots
colores<- c("0"="#e97269","4"="#00BD5C", "5"= "#00C1A7", "6"= "#00BADE", "7"= "#00A6FF", 
            "8" = "#B385FF", "9"= "#EF67EB", "10" = "#FF63B6","11"="#d557a3")
##my plot
sdt.nozero%>%
  ggplot(aes(OPG, Genome_copies_gFaeces))+
  geom_smooth(method = lm)+
  geom_jitter(aes(colour= dpi),position=position_jitter(0.2))+
  scale_fill_manual(values = colores, guide= "none")+
  labs(tag= "a", fill= "DPI")+
  facet_wrap(~Variable)+
  scale_x_log10(name = "log10 (Oocysts/g Faeces) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces)  \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  annotation_logticks()


##my plot residual
sdt.nozero%>%
  group_by(dpi) %>% 
  summarise(residualsM1_mean = mean(na.omit(residualsM1)))%>%
  inner_join(sdt.nozero, by= "dpi")%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(dpi, residualsM1, residualsM1_mean, Variable)%>%
  mutate(residualim= 0)%>%
  ggplot(aes(x= dpi, y= residualsM1, group=Variable))+
  geom_jitter(width = 0.5, shape=21, size=2.5, aes(fill= Variable), alpha= 0.75, color= "black")+
  geom_segment(aes(y= residualsM1_mean, yend= residualim, xend= dpi, colour=Variable, group=Variable), size= 1) +
  geom_point(aes(x = dpi, y = residualsM1_mean, group=Variable), size=4)+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  xlab("Day post infection")+
  scale_y_continuous(name = "Residuals\n (Genome copies/g Faeces)")+
  theme_bw()


##try 2
sdt.nozero%>%
  group_by(Variable,dpi) %>% 
  summarise(residualsM1_mean = mean(residualsM1,rm.na=TRUE))%>%
  inner_join(sdt.nozero, by=c( "dpi","Variable"))%>%
#  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
#  dplyr::select(dpi, residualsM1, residualsM1_mean, Variable)%>%
  mutate(residualim= 0)%>%
  ggplot(aes(x = dpi, y = residualsM1))+
  geom_point(aes(x = dpi, y = residualsM1_mean, colour = Variable),
             size = 5,
             #alpha = .25,
             position = position_dodge(width = 0.75)
  ) +
  geom_jitter(width = 0.1, shape=21, size=2.5, aes(fill= Variable), alpha= 0.75, color= "black")+
  #stat_summary(fun.y = mean(na.omit(residualsM1), geom = "point", size = 5, aes(colour = Variable), position = position_dodge(0.75)) +
  geom_linerange(aes(x = dpi, ymin = residualim, ymax = residualsM1_mean, colour = Variable), 
                 size = 1, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  theme_bw()



  
  
  

##Plot residuals
##Mean residuals for plot 
sdt.nozero%>%
  group_by(dpi) %>% 
  summarise(residualsM1_mean = mean(na.omit(residualsM1)))%>%
  inner_join(sdt.nozero, by= "dpi")%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(dpi, residualsM1, residualsM1_mean)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  mutate(residualim= 0)%>%
  ggplot(aes(x= dpi, y= residualsM1))+
  geom_jitter(width = 0.5, shape=21, size=2.5, aes(fill= dpi), alpha= 0.75, color= "black")+
  scale_fill_manual(values = colores, guide= "none")+
  geom_segment(aes(y= residualsM1_mean, yend= residualim, xend= dpi), color= "black", size= 1) +
  geom_point(aes(x = dpi, y = residualsM1_mean), size=4)+
  #geom_rect(aes(xmin=-0.1,xmax=4.5,ymin=-Inf,ymax=Inf),alpha=0.01,fill="grey")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  xlab("Day post infection")+
  scale_y_continuous(name = "Residuals\n (Genome copies/g Faeces)")+
  labs(tag= "b")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "none")-> B

##To visualize it externally 
#ggsave(filename = "Rplots.pdf", B)
















################################################
################################################
################################################EH WORK

E64_total %>%  ## dplyr::filter(Genome_copies_gFaeces>0) %>%
  transform(logGCopy=log10(Genome_copies_gFaeces+1))%>%
  lm(logGCopy ~ dpi * Variable, data=.) ->
  E64_Model

summary(E64_Model)



E64_total %>%   dplyr::filter(OPG>0) %>%
  transform(logOO=log10(OPG+1))%>%
  lm(logOO ~ dpi * Variable, data=.) ->
  E64_Model


ggpredict(E64_Model, terms = c("dpi", "Variable"), na.action = na.exclude) %>% as.data.frame()


plot(, add.data = TRUE)  + 
  scale_y_log10()

?ggpredict

lm(log10(OPG+1)~dpi+Variable, data=E64_total) -> E64_Model 



ggplot(E64_total,aes(y=log10(OPG),x=Variable, color=dpi))+
  geom_violin()+
  geom_jitter()+
  facet_wrap(~dpi)+
  geom_smooth(method="lm")

E64_total$dpiConti <- as.numeric(as.character(E64_total$dpi))

ggplot(E64_total,aes(y=log10(Genome_copies_gFaeces),x=log10(OPG), color=dpi))+
  geom_jitter()+
  facet_wrap(~Variable)+
  geom_smooth(method="lm", se = FALSE)

E64_total %>% filter(OPG>0 & Genome_copies_gFaeces>0 & !Variable%in%"E88_Primary")%>%
  lm(log10(Genome_copies_gFaeces) ~ log10(OPG)*dpi+Variable, data=.) %>%
  summary()


E64_total %>% filter(OPG>0 & Genome_copies_gFaeces>0 & !Variable%in%"E88_Primary")%>%
  lm(log10(Genome_copies_gFaeces) ~ log10(OPG)+Variable, data=.) %>% 
  residuals() ->
  OOqPCRResids 

E64_total %>% filter(OPG>0 & Genome_copies_gFaeces>0 & !Variable%in%"E88_Primary") %>% cbind(., OOqPCRResids) %>%
  ggplot(aes(x=dpiConti, OOqPCRResids, color=dpi, 
             shape=Variable)) +
  geom_jitter(size=4)










