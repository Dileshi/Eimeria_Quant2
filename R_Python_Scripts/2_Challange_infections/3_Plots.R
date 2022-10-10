##Run scripts on Root of Repo (Eimeria_Quant2)
##Plot challenge infection graph 

library(tidyverse)

##BEGIN with qPCR Data
if(!exists("challenge.DNA")){
  source("R_Python_Scripts/2_Challange_infections/2_qPCR_Data_Prep.R")
}
if(!exists("challenge.DNA")){
  source("R_Python_Scripts/2_Challange_infections/2_qPCR_Data_Prep.R")
}

##Genome copies/g of faeces
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
challenge.DNA%>%
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
#write.csv(x, "Tables/Secondary Infections/Genome_copies_gFaeces_DPI_Comparison_challanged.csv")

stats.test%>%
  dplyr::mutate(y.position = log10(y.position))%>%
  dplyr::mutate(dpi = c("1","2","3","4", "5","6", "7", "8"))-> stats.test


#Let me try my own plot 
challenge.DNA%>%
  dplyr::select(EH_ID, dpi, Genome_copies_gFaeces, Infection)%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8"))%>%
  filter(!is.na(Genome_copies_gFaeces))%>%
  filter(!is.na(Infection))%>%
  ggplot(aes(factor(dpi), Genome_copies_gFaeces+1, colour = factor(dpi),group = Infection))+
  geom_point(aes(shape=factor(Infection)), position=position_jitter(0.05), size=2.5,) + 
  geom_boxplot(fill = "transparent")+
  scale_color_manual(values=c("#f66d9b", "#9561e2", "#6574cd","#3490dc","#4dc0b5","#38c172","#ffed4a","#f6993f","#e3342f"))+
  labs(title="qPCR results of challenge infection", x="Day Post Infection", y="log10 (Genome copies/g Faeces + 1) (qPCR)")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  #scale_x_discrete()+
  stat_pvalue_manual(stats.test, label= "p.adj.signif", x= "dpi", y.position = 100000000000) ->A

#stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> A


#pdf(file = "Figure_secondary_DNA_Test.pdf")
#print(A)
#dev.off()


## Oocysts
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
challenge.DNA%>%
  filter(dpi%in%c("0","1","2","3","4","5","6","7","8"))%>%
  dplyr::select(EH_ID, dpi, OPG)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  dplyr::mutate(OPG= OPG+1)%>% ##To check
  filter(!is.na(OPG))%>% 
  wilcox_test(OPG ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Secondary Infections//Challenge_OPG_DPI_Comparison.csv")

challenge.DNA%>%
  dplyr::select(EH_ID, dpi, OPG)%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8"))%>%
  filter(!is.na(OPG))%>%
  ggplot(aes(factor(dpi), OPG+1, colour = factor(dpi)))+
  geom_point(position=position_jitter(0.2), size=2.5,) + 
  geom_boxplot(fill = "transparent")+
  scale_color_manual(values=c("#f66d9b", "#9561e2", "#6574cd","#3490dc","#4dc0b5","#38c172","#ffed4a","#f6993f","#e3342f"))+
  labs(title="Flotation results of challenge infection", x="Day Post Infection", y="log10 (OPG/g Faeces + 1)")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  #scale_x_discrete()+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> B





##Weight loss 
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
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







