## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Quant2
## Script adapted from Jarquín-Díaz et al. 2022

library(dplyr)
library(wesanderson)

if(!exists("sample.data")){
  source("R_Python_Scripts/1_Primary_Infections/1_Data_Prep.R")
}

if(!exists("sdt")){
  source("R_Python_Scripts/1_Primary_Infections/2_qPCR_Data_Prep.R")
}

##select the completed 25 mice
sdt <- sdt[sdt$EH_ID %in% c("LM0179","LM0180","LM0181", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229",
                            "LM0238", "LM0191","LM0240","LM0246","LM0247","LM0248","LM0194","LM0190","LM0199",
                            "LM0197","LM0198","LM0188","LM0254","LM0255","LM0192"),]

####First plot oocysts count of E. falciformis primary [OPG PLOT]
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","1", "2", "3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi,OPG)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  dplyr::mutate(OPG = OPG+1)%>% ##To check
  #filter(!is.na(OPG))->%>% 
  wilcox_test(OPG ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
#write.csv(x, "Tables/Primary_Infection/OPG_DPI_Primary_Comparison.csv")

sdt%>%
  dplyr::select(EH_ID, dpi, OPG)%>%
  filter(!is.na(OPG))%>%
  ggplot(aes(dpi, OPG+1, fill = dpi))+ #remove colour=dpi to get black outline boxplot
  stat_boxplot(geom = "errorbar", width=0.4, size=1, color="#929292") +
  geom_boxplot(fill = "white",color="#929292")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  scale_fill_manual(values=c("#73F440", "#67DC38", "#439323","#76BA31","#A4DE3F","#D4FB41","#96B83D","#A7DB42","#77B430","#73F440","#00D2A9","#38D279")) +
  geom_point(aes(fill = factor(dpi)), position=position_jitter(0.08), size=3.0,shape=21,stroke=0.35)+
  scale_color_manual(values=c("#73F440", "#67DC38", "#439323","#76BA31","#A4DE3F","#D4FB41","#96B83D","#A7DB42","#77B430","#73F440","#00D2A9","#38D279"))+
  scale_shape_manual(values = c(24))+
  labs(x="Day Post Infection (DPI)", y="log10 (Oocysts/g Faeces + 1) (Oocysts counts)")+
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


##Second plot Eimeria DNA in E falciformis primary infection [Genome copies/g of faeces PLOT]
##Remove contaminant samples: These were mice that had oocyst shed at 0 dpi that are 
#####later relabelled as challenged infection 
sdt_minus_contaminant <- sdt[!((sdt$OPG>0 | sdt$Genome_copies_gFaeces>0) & sdt$dpi ==0),]


#####ADD TO DATA CLEANING SCRPITTTTTTTTTT!!!!!!!!!!!!
##E57aJLX Genome copy value made NA due to error in qPCR experiment
#sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aJLX"] <- NA
#sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aRUV"] <- NA
#sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aHPU"] <- NA
#sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aEFG"] <- NA
#sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aDOU"] <- NA
sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aACJ"] <- NA
sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aLUW"] <- NA
sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aBLZ"] <- NA
sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$labels == "E57aEHR"] <- NA


##remove mouse 255dpi
sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$EH_ID == "LM0255"] <- NA
sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$EH_ID == "LM0248"] <- NA
sdt_minus_contaminant$Genome_copies_gFaeces[sdt_minus_contaminant$EH_ID == "LM0181"] <- NA


sdt_minus_contaminant%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi, Genome_copies_gFaeces)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% 
  #filter(!is.na(Genome_copies_gFaeces))%>%
  wilcox_test(Genome_copies_gFaeces ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
#write.csv(x, "Tables/Primary_Infection/Genome_copies_gFaeces_DPI_Primary_Comparison.csv")

stats.test%>%
  dplyr::mutate(y.position = log10(y.position))%>%
  dplyr::mutate(dpi = c("1","2","3","4", "5","6", "7", "8", "9", "10","11"))-> stats.test

sdt_minus_contaminant%>%
  dplyr::select(EH_ID, dpi, Genome_copies_gFaeces, Infection)%>%
  filter(!is.na(Genome_copies_gFaeces))%>%
  filter(!is.na(Infection))%>%
  ggplot(aes(dpi, Genome_copies_gFaeces+1))+ #remove colour=dpi to get black outline boxplot #replace fill with colour for no balck outline and coloured boxplot
  stat_boxplot(geom = "errorbar", width=0.4, size=1, color="#929292") +
  geom_boxplot(fill = "white",color="#929292")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  scale_fill_manual(values=c("#6CB4EE", "#0066B2", "#89CFF0","#318CE7","#1E90FF","#4D4DFF","#3457D5","#00BFFF","#2A52BE","#0CAFFF","#6495ED","#1F75FE")) +
  geom_point(aes(shape=Infection,fill = factor(dpi)), position=position_jitter(0.05), size=3.0,)+
  scale_color_manual(values=c("#6CB4EE", "#0066B2", "#89CFF0","#318CE7","#1E90FF","#4D4DFF","#3457D5","#00BFFF","#2A52BE","#0CAFFF","#6495ED","#1F75FE"))+
  scale_shape_manual(values = c(24, 21))+
  labs(x="Day Post Infection (DPI)", y="log10 (Genome copies/g Faeces + 1) (Eimeria DNA)")+
  labs(tag= "b", shape= "qPCR results of E88 Primary Infection")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  guides(colour="none")+
  guides(fill="none")+
  theme(text = element_text(size=13), axis.title.x = element_blank(), legend.position = "top")+
  annotation_logticks(sides = "l")+
  stat_pvalue_manual(stats.test, label= "p.adj.signif", x= "dpi", y.position = 100000000000) ->B

print(B)



##Third plot Weightloss in E falciformis primary infection 
####calculate relative weight loss using function
sdt_minus_contaminant$weightloss <- ((sdt_minus_contaminant$weight_dpi0 - sdt_minus_contaminant$weight)/sdt_minus_contaminant$weight_dpi0)*100

##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt_minus_contaminant%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi, weightloss)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  wilcox_test(weightloss ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Primary_Infection/Weightloss_DPI_Primary_Comparison.csv")

sdt_minus_contaminant%>%
  filter(!is.na(weightloss))%>%
  ggplot(aes(dpi, weightloss, fill =dpi))+ #remove colour=dpi to get black outline boxplot #replace fill with colour for no balck outline and coloured boxplot
  stat_boxplot(geom = "errorbar", width=0.4, size=1, color="#929292") +
  geom_boxplot(fill = "white",color="#929292")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  scale_fill_manual(values=c("#E63A5F", "#ED5487", "#E03B72","#EF6668","#F6C0D2","#ED67B3","#EC49B2","#F195C6","#F080AA","#EC4395","#EC47A8","#EC407E")) +
  geom_point(aes(fill = factor(dpi)), position=position_jitter(0.05), size=3.0,shape=21,stroke=0.35)+
  scale_color_manual(values=c("#E63A5F", "#ED5487", "#E03B72","#EF6668","#F6C0D2","#ED67B3","#EC49B2","#F195C6","#F080AA","#EC4395","#EC47A8","#EC407E"))+
  labs(x="Day Post Infection (DPI)", y=" Weightloss (%)")+
  labs(tag= "c", shape= "Percentage of Weightloss of E88 Primary Infection")+
  theme_bw()+
  guides(colour="none")+
  guides(fill="none")+
  theme(text = element_text(size=13), axis.title.x = element_blank(), legend.position = "top")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> C
print(C)

##Figure 2: Course of Eimeria Infection in genome copies, OPG, and weight loss
pdf(file = "Figures/Figure_33", width = 10, height = 15)
grid.arrange(A,B,C)
dev.off()
###rm(A,B,C, x, stats.test)



















