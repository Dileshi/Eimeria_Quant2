##Plot statistical model for primary (E64,E88) and secondary infection (E64)
#Compares effect of being primary/secondary E64:E64 & E88:E64 for DNA, Weight, Oocyst

library(dplyr)

########IMPORT E88 PRIMARY DATA (E88_Prim)
if(!exists("sample.data")){
  source("R_Python_Scripts/1_Primary_Infections/1_Data_Prep.R")
}
if(!exists("sdt")){
  source("R_Python_Scripts/1_Primary_Infections/2_qPCR_Data_Prep.R")
}
##select the completed 18 mice of sdt dataframe 
E88_Prim  <- sdt[sdt$EH_ID %in% c("LM0179","LM0180","LM0181", "LM0182","LM0184","LM0185","LM0186","LM0188","LM0190","LM0191",
                            "LM0192", "LM0194","LM0195","LM0197","LM0198","LM0199","LM0244","LM0246"),]
#sdt <- sdt[sdt$EH_ID %in% c("LM0179","LM0180","LM0181", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229",
                            #"LM0238", "LM0191","LM0240","LM0244","LM0246","LM0247","LM0248","LM0194","LM0190","LM0199",
                            #"LM0197","LM0198","LM0188","LM0254","LM0255","LM0192"),]
E88_Prim  <- E88_Prim [c(1,3:6,8,44,46,48,63,66:69,73)]


########IMPORT E64 PRIMARY DATA (E64_Prim)
E64_Prim <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Quant/master/data/Eimeria_quantification_Inf_exp_data.csv")
E64_Prim <- E64_Prim %>% rename(Cq_mean = 'Cq.Mean')
E64_Prim <- E64_Prim[c(5,8)]
E64_Prim %>%   dplyr::distinct(Sample, .keep_all = TRUE)-> E64_Prim
E64_Prim <- E64_Prim %>% rename(labels = 'Sample')
E64_Prim_OPG <- read.csv("Output_Data/E64_Primary_Victor_all.csv") #OPG, Genome copy calulations + other data
E64_Prim<- left_join(E64_Prim_OPG, E64_Prim, by="labels") ## Add qPCR data
E64_Prim <- E64_Prim[c(1:6,9,23,55:57,60,61)]
rm(E64_Prim_OPG)
E64_Prim$labels<-paste0('E57a', E64_Prim$labels)

########IMPORT E64 Secondary DATA (E64_Secon)
if(!exists("DNA.data")){
  source("R_Python_Scripts/2_Challange_infections/2_qPCR_Data_Prep.R")
}
if(!exists("challenge.DNA")){
  source("R_Python_Scripts/2_Challange_infections/2_qPCR_Data_Prep.R")
}
E64_E64 <- challenge.DNA
rm(challenge.DNA)


########
##E88 secondary samples in sdt dataframe 
E88_E88 <- sdt[sdt$challenge_infection == "E88", ]



#E88:E88 + E88:E64
rm(sdt)





#https://mspeekenbrink.github.io/sdam-r-companion/factorial-anova.html 
#








