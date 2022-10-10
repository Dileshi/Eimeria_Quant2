library (dplyr)
E64_prim = read.csv("https://raw.githubusercontent.com/derele/Eimeria_Quant/master/data/Eimeria_quantification_Inf_exp_data.csv")
E64_sam = select(E64_prim, c('Sample','Cq'))
##remove 'Neg_Ctrl' and 'Pos_Ctrl'
E64_sam<-E64_sam[!(E64_sam$Sample=="Neg_Ctrl" | E64_sam$Sample=="Pos_Ctrl"),]
##add E57a infront
E64_sam$Sample <- paste('E57a', E64_sam$Sample, sep='')
All = read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Quant2_E57.csv")
All_sel = select(All, c('EH_ID','labels'))
##rename and merge 
E64_sam %>% 
  rename(
    labels = Sample,
  )-> E64_sam
total <- merge(E64_sam, All_sel,by="labels")
##primary infected mice in Quant1 experiment 
##LM0204, LM0206 - LM0226
##list which mice were carried to secondary infection
challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")
##select LM0204, LM0206 - LM0226 primary infected mice
chall_total = merge(total, challenge, by ="EH_ID")
##primary infected mice in Quant1 project not avaible for secondary infection
##E88 primary mice I did are either re-infected with E64 or are uninfected
##Decided to continue with random E64 (E64:E64) or E88(E88:E88) secondary infected mice
#In challange infected dataframe select E64:E64 mice and E88:E88 mice
E64_challenge <- challenge %>% filter(primary_infection== 'E64' & challenge_infection == 'E64')
E88_challenge <- challenge %>% filter(primary_infection== 'E88' & challenge_infection == 'E88')
##select only E57 experiment
E64_chall57 <- E64_challenge %>% filter(experiment== 'E57'& infection == 'challenge')
## select E88 samples for secondary infection 
##write.csv(E88_challenge,"Desktop/E88_challenged.csv", row.names = FALSE)
## select E64 samples for secondary infection 
write.csv(E64_chall57,"Desktop/E64_challenged.csv", row.names = FALSE)