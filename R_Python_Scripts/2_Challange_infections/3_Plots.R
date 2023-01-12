##Run scripts on Root of Repo (Eimeria_Quant2)
##Plot challenge infection graph 

library(tidyverse)
library(ZIR)

##BEGIN with qPCR Data
if(!exists("DNA.data")){
  source("R_Python_Scripts/2_Challange_infections/2_qPCR_Data_Prep.R")
}
if(!exists("challenge.DNA")){
  source("R_Python_Scripts/2_Challange_infections/2_qPCR_Data_Prep.R")
}


challenge.DNA2 <- challenge.DNA[challenge.DNA$EH_ID %in% c("LM0258", "LM0259","LM0260","LM0262",
                                                           "LM0265","LM0269","LM0270","LM0279","LM0280", "LM0284"),]





##Secondar Infection plots 

##upload challnege infection file 


###Add to data cleaning script following 
challenge_DNA$OPG[challenge_DNA$labels == "E57byLQU"] <- NA
challenge_DNA$OPG[challenge_DNA$labels == "E57byNTU"] <- NA

## Oocysts
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
challenge_DNA%>%
  #filter(dpi%in%c("0","1", "2", "3","4", "5","6", "7", "8"))%>%
  #dplyr::select(EH_ID, dpi,OPG)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  #dplyr::mutate(OPG= OPG+1)%>% ##To check
  filter(!is.na(OPG))%>% 
  wilcox_test(OPG ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

#challenge.DNA2 <-challenge.DNA[!(challenge.DNA$dpi=="1" | challenge.DNA$dpi=="2" | challenge.DNA$dpi=="3"),]


challenge_DNA%>%
  filter(dpi%in%c("0","4", "5","6", "7", "8"))%>%
  dplyr::select(EH_ID, dpi,OPG)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  filter(!is.na(OPG))%>% 
  wilcox_test(OPG ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

ziw(challenge.DNA$OPG, challenge.DNA$dpi, perm = TRUE)
##$p.value[1] 0.0143
##$statistics[1] 2.477437



##Save statistical analysis
x <- stats.test
x$groups<- NULL
#write.csv(x, "Tables/Secondary_Infections/OPG_DPI_Challenged_Comparison_15mice.csv")

challenge_DNA%>%
  dplyr::select(EH_ID, dpi, OPG)%>%
  ggplot(aes(factor(dpi), OPG+1, fill = factor(dpi)))+ #remove colour=dpi to get black outline boxplot
  stat_boxplot(geom = "errorbar", width=0.4, size=1, color="#929292") +
  geom_boxplot(fill = "white",color="#929292")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  scale_fill_manual(values=c("#73F440", "#67DC38", "#439323","#76BA31","#A4DE3F","#D4FB41","#96B83D","#A7DB42","#77B430","#73F440","#00D2A9","#38D279")) +
  geom_point(aes(fill = factor(dpi)), position=position_jitter(0.05), size=3.0,shape=21,stroke=0.35)+
  scale_color_manual(values=c("#73F440", "#67DC38", "#439323","#76BA31","#A4DE3F","#D4FB41","#96B83D","#A7DB42","#77B430","#73F440","#00D2A9","#38D279"))+
  scale_shape_manual(values = c(24))+
  labs(x="Day Post Infection", y="log10 (Oocysts/g Faeces + 1) (OPG)")+
  labs(tag= "a", shape= "Flotation results of Primary Infection")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  guides(colour="none")+
  guides(fill="none")+
  theme(text = element_text(size=13), axis.title.x = element_blank(), legend.position = "top")+
  annotation_logticks(sides = "l")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> A
print(A)
#pdf(file = "E64_chall_test_2.pdf", width = 10, height = 15)
#grid.arrange(A,B)
#dev.off()




##Genome copies/g of faeces
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
challenge_DNA%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8"))%>%
  dplyr::select(EH_ID, dpi, Genome_copies_gFaeces)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  filter(!is.na(Genome_copies_gFaeces))%>%
  wilcox_test(Genome_copies_gFaeces ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Secondary_Infections/Genome_copies_gFaeces_DPI_Comparison_challanged.csv")

stats.test%>%
  dplyr::mutate(y.position = log10(y.position))%>%
  dplyr::mutate(dpi = c("1","2","3","4", "5","6", "7", "8"))-> stats.test

#challenge.DNA[,'dpi'] <- as.factor(as.character(challenge.DNA[,'dpi']))
#Let me try to plot 
challenge_DNA%>%
  dplyr::select(EH_ID, dpi, Genome_copies_gFaeces, Infection)%>%
  filter(!is.na(Genome_copies_gFaeces))%>%
  filter(!is.na(Infection))%>%
  ggplot(aes(factor(dpi), Genome_copies_gFaeces+1))+ #remove colour=dpi to get black outline boxplot #replace fill with colour for no balck outline and coloured boxplot
  stat_boxplot(geom = "errorbar", width=0.4, size=1, color="#929292") +
  geom_boxplot(fill = "white",color="#929292")+
  scale_fill_manual(values=c("#6CB4EE", "#0066B2", "#89CFF0","#318CE7","#1E90FF","#4D4DFF","#3457D5","#00BFFF","#2A52BE","#0CAFFF","#6495ED","#1F75FE")) +
  geom_point(aes(shape=Infection,fill = factor(dpi)), position=position_jitter(0.05), size=3.0,)+
  scale_color_manual(values=c("#6CB4EE", "#0066B2", "#89CFF0","#318CE7","#1E90FF","#4D4DFF","#3457D5","#00BFFF","#2A52BE","#0CAFFF","#6495ED","#1F75FE"))+
  scale_shape_manual(values = c(24, 21))+
  labs(x="Day Post Infection", y="log10 (Genome copies/g Faeces + 1) (qPCR)")+
  labs(tag= "a", shape= "qPCR results of E88 Secondary Infection")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  theme_bw()+
  guides(colour="none")+
  guides(fill="none")+
  theme(text = element_text(size=13), axis.title.x = element_blank(), legend.position = "top")+
  annotation_logticks(sides = "l")+
  stat_pvalue_manual(stats.test, label= "p.adj.signif", x= "dpi", y.position = 100000000000) ->B
print(B)
#pdf(file = "E64_Chall_test_1.pdf")
#print(A)
#dev.off()


##Weight loss 
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
#challenge.DNA$weightloss <- ((challenge.DNA$weight_dpi0 - challenge.DNA$weight)/challenge.DNA$weight_dpi0)*100

##Weight loss 
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
challenge_DNA%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8"))%>%
  dplyr::select(EH_ID, dpi, weightloss)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison
  filter(!is.na(weightloss))%>% 
  wilcox_test(weightloss ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Secondary_Infections/Weightloss_DPI_Challenge_Comparison.csv")

challenge_DNA%>%
  dplyr::select(EH_ID, dpi, weightloss)%>%
  filter(!is.na(weightloss))%>%
  ggplot(aes(factor(dpi), weightloss, fill =factor(dpi)))+ #remove colour=dpi to get black outline boxplot #replace fill with colour for no balck outline and coloured boxplot
  stat_boxplot(geom = "errorbar", width=0.4, size=1, color="#929292") +
  geom_boxplot(fill = "white",color="#929292")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  scale_fill_manual(values=c("#E63A5F", "#ED5487", "#E03B72","#EF6668","#F6C0D2","#ED67B3","#EC49B2","#F195C6","#F080AA","#EC4395","#EC47A8","#EC407E")) +
  geom_point(aes(fill = factor(dpi)), position=position_jitter(0.05), size=3.0,shape=21,stroke=0.35)+
  scale_color_manual(values=c("#E63A5F", "#ED5487", "#E03B72","#EF6668","#F6C0D2","#ED67B3","#EC49B2","#F195C6","#F080AA","#EC4395","#EC47A8","#EC407E"))+
  scale_shape_manual(values = c(24))+
  labs(x="Day Post Infection", y=" Weightloss (%)")+
  theme_bw()+
  guides(colour="none")+
  guides(fill="none")+
  theme(text = element_text(size=13), axis.title.x = element_blank(), legend.position = "top")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> C
print(C)

##Figure 2: Course of Eimeria Infection in genome copies, OPG, and weight loss
pdf(file = "Figure_challnege.pdf", width = 10, height = 15)
grid.arrange(A,B,C)
dev.off()
###rm(A,B,C, x, stats.test)















challenge.DNA%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8"))%>%
  dplyr::select(EH_ID, dpi, relative_weight)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  wilcox_test(relative_weight ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Secondary Infections/Weightloss_DPI_Comparison_challenge.csv")

challenge.DNA%>%
  dplyr::select(EH_ID, dpi, relative_weight)%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8"))%>%
  #filter(!is.na(relative_weight))%>%
  ggplot(aes(factor(dpi), relative_weight, colour = factor(dpi)))+
  geom_point(position=position_jitter(0.2), size=2.5,) + 
  geom_boxplot(fill = "transparent")+
  scale_color_manual(values=c("#f66d9b", "#9561e2", "#6574cd","#3490dc","#4dc0b5","#38c172","#ffed4a","#f6993f","#e3342f"))+
  labs(title="Weight loss of challenge infection", x="Day Post Infection", y="Weight loss (%)")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  #scale_x_continuous(labels = as.character(dpi),breaks = dpi)+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> C

##Figure 4: Course of Eimeria Infection in genome copies, OPG, and weight loss
pdf(file = "Figures/Figure_4_challenge2.pdf", width = 10, height = 15)
grid.arrange(A,B,C)
dev.off()
#rm(A,B,C, x, stats.test)







