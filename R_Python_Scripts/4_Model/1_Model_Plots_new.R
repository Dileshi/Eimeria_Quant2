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
##select the completed 26 mice of sdt dataframe 

E88_Prim <- sdt[sdt$EH_ID %in% c("LM0179","LM0180","LM0181", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229",
                                 "LM0238", "LM0191","LM0240","LM0244","LM0246","LM0247","LM0248","LM0194","LM0190","LM0199",
                                 "LM0197","LM0198","LM0188","LM0254","LM0255","LM0192"),]
E88_Prim$weightloss <-  ((E88_Prim$weight_dpi0 - E88_Prim$weight)/E88_Prim$weight_dpi0)*100

E88_Prim = E88_Prim[,c("EH_ID","labels","dpi","weightloss","OPG","Genome_copies_gFaeces")]

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
E64_Prim$weightloss <- ((E64_Prim$weight_dpi0 - E64_Prim$weight)/E64_Prim$weight_dpi0)*100

E64_Prim %>%select(labels, EH_ID, dpi,OPG,weightloss,Genome_copies_gFaeces) -> E64_Prim


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

challenge_DNA %>%select(labels, EH_ID, dpi,OPG,weightloss,Genome_copies_gFaeces) -> E64_E64


####combine E64 primary and E64 challenge into one data frame 
##first add column primary
E64_Prim$Variable <- paste0('E64_Prim')
##first add column primary
E88_Prim$Variable <- paste0('E88_Prim')
##first add column secondary
E64_E64$Variable <- paste0('E64_Sec')
#merge E64_Primary with E64:E64 secondary
E64_total <- rbind(E64_Prim,E64_E64)
E64_total <- rbind(E64_total,E88_Prim)

##contaminated 8 mice
E64_total%>%
  filter(dpi==0 & Variable== 'E88_Prim' & (OPG>0 | Genome_copies_gFaeces>0)) -> contaminant
#make vector with EH_ID
EH_ID_Contaminant <- contaminant$EH_ID 

E64_total%>%
  filter(E64_total$EH_ID %in% EH_ID_Contaminant)%>%
  mutate(Variable = replace(Variable, Variable == 'E88_Prim', 'E88_Sec'))->E64_total_conta

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
legend_title <- "dpi"
ggplot(E64_total,aes(y=log10(Genome_copies_gFaeces+1),x=Variable, color=factor(dpi)))+
  geom_violin()+
  scale_color_manual(legend_title,values=c("#6CB4EE", "#0066B2", "#89CFF0","#318CE7","#1E90FF",
                              "#4D4DFF","#3457D5","#00BFFF","#2A52BE","#0CAFFF","#6495ED","#1F75FE"))+
  stat_summary(fun.y =median, mult=1, 
               geom="pointrange", color="black")+
  geom_jitter()+
  facet_wrap(~dpi)+
  labs(x="Infection", y="log10 (Genome copies/g Faeces + 1) (qPCR)")+
  scale_fill_discrete(name = "dpi")+
  geom_smooth(method="lm")
#OPG
ggplot(E64_total,aes(y=log10(OPG+1),x=Variable, color=factor(dpi)))+
  geom_violin()+
  scale_color_manual(legend_title, values=c("#73F440", "#67DC38", "#439323","#76BA31","#A4DE3F",
                              "#D4FB41","#96B83D","#A7DB42","#77B430","#73F440","#00D2A9","#38D279"))+
  stat_summary(fun.y =median, mult=1, 
               geom="pointrange", color="black")+
  geom_jitter()+
  facet_wrap(~dpi)+
  labs(x="Infection", y="log10 (Oocysts/g Faeces + 1) (Flotations)")+
  geom_smooth(method="lm")

#weightloss
ggplot(E64_total,aes(y=weightloss,x=Variable, color=factor(dpi)))+
  geom_violin()+
  stat_summary(fun.y =median, mult=1, 
               geom="pointrange", color="black")+
  scale_color_manual(legend_title, values=c("#E63A5F", "#ED5487", "#E03B72","#EF6668","#F6C0D2","#ED67B3","#EC49B2","#F195C6","#F080AA","#EC4395","#EC47A8","#EC407E"))+
  geom_jitter()+
  facet_wrap(~dpi)+
  labs(x="Infection", y="Weight loss (%)")+
  geom_smooth(method="lm")
################################################ 
#Residual plot

sdt <- E64_total
sdt$OPG[sdt$labels == "E57aDMV"] <- NA
sdt$OPG[sdt$labels == "E57byFPV"] <- NA

sdt%>%
  filter(!(OPG== 0))%>%
  filter(!(Genome_copies_gFaeces== 0))-> sdt.nozero

##Model 1: Genome copies/g faeces modeled by OPG
DNAbyOPG <- lm(log10(Genome_copies_gFaeces)~(log10(OPG)),
               data = sdt.nozero, na.action = na.exclude)
summary(DNAbyOPG)

### Plot and extract estimates
require(sjPlot)

sdt.nozero$predictedM1 <- predict(DNAbyOPG)   # Save the predicted values
sdt.nozero$residualsM1 <- residuals(DNAbyOPG) # Save the residual values

##Assign the colors for dpi and keep consistency with previous plots
colores<- c("0"="#FFED4F","4"="#F3993E", "5"= "#E33530","6"="#EE6C9B", "7"= "#9561E2", "8"= "#6574CD", 
            "9" = "#3490DC", "10"= "#53C0B5", "11" = "#58C172")
#colores<- c("0"="#EE6C9B","4"="#9561E2", "5"= "#3490DC","6"="#6574CD", "7"= "#53C0B5", "8"= "#58C172", 
 #           "9" = "#FFED4F", "10"= "#F3993E", "11" = "#E33530")
#colores<- c("0"="#E25151","4"="#E161AA", "5"= "#D84B94","6"="#9C88F8", "7"= "#8D51E7", "8"= "#63DBFA", 
#           "9" = "#237EE0", "10"= "#9BD33B", "11" = "#489228")
##my plot
sdt.nozero%>%
  ggplot(aes(OPG, Genome_copies_gFaeces))+
  geom_smooth(method = lm,color= "black")+
  geom_point(aes(fill = factor(dpi), shape = Variable), size=4.0,stroke=0.35,color="#424242")+
  scale_fill_manual(values = colores)+
  scale_shape_manual(values = c(24, 25, 21,22))+
  labs(tag= "a", fill= "DPI")+
  scale_x_log10(name = "log10 (Oocysts/g Faeces) \n (Oocyst counts)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces)  \n (Eimeria DNA)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  annotation_logticks()


sdt.nozero <- subset(sdt.nozero, !(dpi %in% c(0))) 



  
legend_title2 <- "Infection\nCondition"
##my plot residual
sdt.nozero%>%
  group_by(Variable,dpi) %>% 
  summarise(residualsM1_mean = mean(residualsM1,rm.na=TRUE))%>%
  inner_join(sdt.nozero, by=c( "dpi","Variable"))%>%
  mutate(residualim= 0)%>%
  ggplot(aes(x = dpi, y = residualsM1, color=factor(Variable)))+
  scale_color_manual(legend_title2, values=c("#3490DC", "#53C0B5", "#F3993E","#E33530"))+
  scale_fill_manual(legend_title2, values=c("#3490DC", "#53C0B5", "#F3993E","#E33530"))+
  geom_point(aes(x = dpi, y = residualsM1_mean,fill=Variable),
             size = 5,
             shape=21,stroke=1,#alpha = .25,
             position = position_dodge(width = 0.75)) +
  geom_jitter(width = 0.1, shape=21, size=2.5, alpha= 0.75,aes(fill=Variable))+
  #stat_summary(fun.y = mean(na.omit(residualsM1), geom = "point", size = 5, aes(colour = Variable), position = position_dodge(0.75)) +
  geom_linerange(aes(x = dpi, ymin = residualim, ymax = residualsM1_mean, colour = Variable), 
                 size = 1, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  xlab("Days post infection (dpi)")+
  scale_y_continuous(name = "Residuals\n (Genome copies/g Faeces)")+
  theme_bw()+
  theme(legend.position="bottom")
 
  
  





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
  ggplot(aes(x=dpi, OOqPCRResids, color=dpi, 
             shape=Variable)) +
  geom_jitter(size=4)










