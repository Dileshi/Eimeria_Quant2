##Plot challenge infection graph 

library(RColorBrewer)

if(!exists("challenge")){
  source("R_Python_Scripts/2_Challange_infections/2_Data_Prep_chal.R")
}

challenge$dpi<-sapply(challenge$dpi, as.factor)
challenge$EH_ID<-sapply(challenge$EH_ID, as.factor)
#challenge <- challenge%>%dplyr::mutate(OPG= OPG+1)

## Oocysts
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
challenge%>%
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
write.csv(x, "Tables/Challenge_OPG_DPI_Comparison.csv")


challenge%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8"))%>%
  dplyr::select(EH_ID, dpi,OPG, challenge_infection)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= OPG+1, color=dpi))+
  scale_y_log10("log10 (Oocysts/g Faeces + 1) (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot()+
  geom_point(shape=21, position=position_jitter(0.2), size=2.5, aes(shape= challenge_infection, fill=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  labs(tag= "b")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "none")+
  annotation_logticks(sides = "l")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> B


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

challenge%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8"))%>%
  dplyr::select(EH_ID, dpi,OPG, challenge_infection)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  filter(!is.na(OPG))%>% 
  ggplot(aes(x= dpi, y= OPG+1, color=dpi))+
  scale_y_log10("log10 (Oocysts/g Faeces + 1) (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot()+
  geom_point(position=position_jitter(0.2), size=5.5, aes(shape= challenge_infection, fill=dpi))+
  scale_fill_manual(values=cbbPalette)+
  scale_shape_manual(values = c(21, 24))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  labs(tag= "a", shape= "infection status")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top")+
  annotation_logticks(sides = "l")+
  guides(scale="none")-> A


pdf(file = "Figures/Figure_secondary.pdf")
print(A)
dev.off()







