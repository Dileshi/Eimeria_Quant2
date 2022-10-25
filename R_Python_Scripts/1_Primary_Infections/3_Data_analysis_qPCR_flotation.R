## RUN SCRIPTS FROM THE ROOT OF THE REPOSITORY: Eimeria_Quant2
## Script adapted from Victor's 3_Data_analysis_qPCR_flotation


### Code to analyses
## 1) Correlation among qPCR and oocyst flotation quantification
### library(ggsci)

library(dplyr)
library(wesanderson)

if(!exists("sample.data")){
  source("R_Python_Scripts/1_Primary_Infections/1_Data_Prep.R")
  }

if(!exists("sdt")){
    source("R_Python_Scripts/1_Primary_Infections/2_qPCR_Data_Prep.R")
  }

##select the completed 26 mice <- correction 18 mice not 26 mice 
sdt <- sdt[sdt$EH_ID %in% c("LM0179","LM0180","LM0181", "LM0182","LM0184","LM0185","LM0186","LM0188","LM0190","LM0191",
                            "LM0192", "LM0194","LM0195","LM0197","LM0198","LM0199","LM0244","LM0246"),]
##sdt <- sdt[sdt$EH_ID %in% c("LM0179","LM0180","LM0181", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229",
                             ##"LM0238", "LM0191","LM0240","LM0244","LM0246","LM0247","LM0248","LM0194","LM0190","LM0199",
                             ##"LM0197","LM0198","LM0188","LM0254","LM0255","LM0192"),]

#sdt_breif <- sdt[c(1,3,7,8,9,10,50,57,63,64,65,66,67,68,69,70,71,72)]

###Let's start plotting and analysing the data!
### 1) Course of infection 
##Genome copies/g of faeces
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
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
#write.csv(x, "Tables/Primary_Infection/Genome_copies_gFaeces_DPI_Primary_Comparison.csv")

stats.test%>%
  dplyr::mutate(y.position = log10(y.position))%>%
  dplyr::mutate(dpi = c("1","2","3","4", "5","6", "7", "8", "9", "10","11"))-> stats.test

sdt%>%
  dplyr::select(EH_ID, dpi, Genome_copies_gFaeces, Infection)%>%
  filter(!is.na(Genome_copies_gFaeces))%>%
  filter(!is.na(Infection))%>%
  ggplot(aes(dpi, Genome_copies_gFaeces+1))+ #remove colour=dpi to get black outline boxplot #replace fill with colour for no balck outline and coloured boxplot
  geom_boxplot(fill = "white")+
  scale_fill_manual(values=c("#6CB4EE", "#0066B2", "#89CFF0","#318CE7","#1E90FF","#4D4DFF","#3457D5","#00BFFF","#2A52BE","#0CAFFF","#6495ED","#1F75FE")) +
  geom_point(aes(shape=Infection,fill = factor(dpi)), position=position_jitter(0.05), size=3.0,)+
  scale_color_manual(values=c("#6CB4EE", "#0066B2", "#89CFF0","#318CE7","#1E90FF","#4D4DFF","#3457D5","#00BFFF","#2A52BE","#0CAFFF","#6495ED","#1F75FE"))+
  scale_shape_manual(values = c(24, 21))+
  labs(x="Day Post Infection", y="log10 (Genome copies/g Faeces + 1) (qPCR)")+
  labs(tag= "a", shape= "qPCR results of Primary Infection")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  theme_bw()+
  guides(colour="none")+
  guides(fill="none")+
  theme(text = element_text(size=13), axis.title.x = element_blank(), legend.position = "top")+
  annotation_logticks(sides = "l")+
  stat_pvalue_manual(stats.test, label= "p.adj.signif", x= "dpi", y.position = 100000000000) ->A

#pdf(file = "E88_Prim_test_1.pdf")
#print(A)
#dev.off()



### Oocysts
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","1", "2", "3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi,OPG)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  dplyr::mutate(OPG= OPG+1)%>% ##To check
  filter(!is.na(OPG))->#%>% 
  wilcox_test(OPG ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")#-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
#write.csv(x, "Tables/Primary_Infection/OPG_DPI_Primary_Comparison.csv")

sdt%>%
  dplyr::select(EH_ID, dpi, OPG)%>%
  filter(!is.na(OPG))%>%
  ggplot(aes(dpi, OPG+1, colour = dpi))+ #remove colour=dpi to get black outline boxplot
  geom_boxplot(fill = "white")+
  scale_fill_manual(values=c("#73F440", "#67DC38", "#439323","#76BA31","#A4DE3F","#D4FB41","#96B83D","#A7DB42","#77B430","#73F440","#00D2A9","#38D279")) +
  geom_point(aes(fill = factor(dpi)), position=position_jitter(0.05), size=3.0,)+
  scale_color_manual(values=c("#73F440", "#67DC38", "#439323","#76BA31","#A4DE3F","#D4FB41","#96B83D","#A7DB42","#77B430","#73F440","#00D2A9","#38D279"))+
  scale_shape_manual(values = c(24))+
  labs(x="Day Post Infection", y="log10 (Oocysts/g Faeces + 1) (qPCR)")+
  labs(tag= "a", shape= "Flotation results of Primary Infection")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  theme_bw()+
  guides(colour="none")+
  guides(fill="none")+
  theme(text = element_text(size=13), axis.title.x = element_blank(), legend.position = "top")+
  annotation_logticks(sides = "l")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> B

#pdf(file = "E88_Prim_test_2.pdf", width = 10, height = 15)
#grid.arrange(A,B)
#dev.off()


####calculate relative weight loss using function
sdt$weightloss <- ((sdt$weight_dpi0 - sdt$weight)/sdt$weight_dpi0)*100

##Weight loss 
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10","11"))%>%
  dplyr::select(EH_ID, dpi, weightloss)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  wilcox_test(weightloss ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
#x <- stats.test
#x$groups<- NULL
#write.csv(x, "Tables/Primary_Infection/Weightloss_DPI_Primary_Comparison.csv")

sdt%>%
  dplyr::select(EH_ID, dpi, weightloss)%>%
  filter(!is.na(weightloss))%>%
  ggplot(aes(dpi, weightloss, colour =dpi))+ #remove colour=dpi to get black outline boxplot #replace fill with colour for no balck outline and coloured boxplot
  geom_boxplot(fill = "white")+
  scale_fill_manual(values=c("#E63A5F", "#ED5487", "#E03B72","#EF6668","#F6C0D2","#F6C0D2","#EC49B2","#F195C6","#F080AA","#EC4395","#EC47A8","#EC407E")) +
  geom_point(aes(colour = factor(dpi)), position=position_jitter(0.05), size=3.0,)+
  scale_color_manual(values=c("#E63A5F", "#ED5487", "#E03B72","#EF6668","#F6C0D2","#F6C0D2","#EC49B2","#F195C6","#F080AA","#EC4395","#EC47A8","#EC407E"))+
  scale_shape_manual(values = c(24))+
  labs(x="Day Post Infection", y=" Weightloss")+
  labs(tag= "a", shape= "Percentage of Weightloss of Primary Infections")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.6)+
  theme_bw()+
  guides(colour="none")+
  guides(fill="none")+
  theme(text = element_text(size=13), axis.title.x = element_blank(), legend.position = "top")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> C


##Figure 2: Course of Eimeria Infection in genome copies, OPG, and weight loss
pdf(file = "Figures/Figure_2A_test_3.pdf", width = 10, height = 15)
grid.arrange(A,B,C)
dev.off()
###rm(A,B,C, x, stats.test)






























### 2) Correlation among Eimeria quantification methods
###Aiming to state a quantitative measure of the relationship between both measurements 
###Generate a data set without samples with zero OPG counts and/or zero genome copies 
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

##Plot model
##Assign the colors for dpi and keep consistency with previous plots
colores<- c("0"="#B7A034", "1"= "#00BD5C", "6"= "#00C1A7", "7"= "#00BADE", 
            "8" = "#00A6FF", "9"= "#B385FF", "10" = "#EF67EB", "11" = "#FF63B6")

####Genome copies modeled by OPGs 
sdt.nozero%>%
  ggplot(aes(OPG, Genome_copies_gFaeces))+
  geom_smooth(method = lm, col= "black")+
  scale_x_log10(name = "log10 (Oocysts/g Faeces) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces)  \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  scale_fill_manual(values = colores, guide= "none")+
  labs(tag= "a", fill= "DPI")+
  theme_bw()+
  #stat_cor(label.y = 5,  label.x = 6, method = "spearman",
  #         aes(label= paste("rho","'='", ..r.., ..p.label.., sep= "~` `~")))+
  theme(text = element_text(size=16), legend.position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  annotation_logticks()-> A

##To visualize it 
#ggsave(filename = "Rplots.pdf", A)

##Plot residuals
##Mean residuals for plot 
sdt.nozero%>%
  group_by(dpi) %>% 
  summarise(residualsM1_mean = mean(na.omit(residualsM1)))%>%
  inner_join(sdt.nozero, by= "dpi")%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10", "11"))%>%
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

##Model 2: Genome copies/g faeces modeled by OPG without DPI interaction
DNAbyOPG_dpi <- lm(log10(Genome_copies_gFaeces)~log10(OPG)+dpi,
                   data = sdt.nozero, na.action = na.exclude)
summary(DNAbyOPG_dpi)

##Model 3: Genome copies/g faeces modeled by OPG with DPI interaction
DNAbyOPGxdpi <- lm(log10(Genome_copies_gFaeces)~log10(OPG)*dpi,
                   data = sdt.nozero, na.action = na.exclude)
summary(DNAbyOPGxdpi)

### Plot and extract estimates
sjPlot:: tab_model(DNAbyOPGxdpi, 
                   terms = c("log10(OPG):dpi0","log10(OPG):dpi1", "log10(OPG):dpi6", "log10(OPG):dpi7", 
                             "log10(OPG):dpi8", "log10(OPG):dpi9", "log10(OPG):dpi10", "log10(OPG):dpi11"))

plot_model(DNAbyOPGxdpi, terms = c("log10(OPG):dpi0","log10(OPG):dpi1", "log10(OPG):dpi6", "log10(OPG):dpi7", 
                                   "log10(OPG):dpi8", "log10(OPG):dpi9", "log10(OPG):dpi10", "log10(OPG):dpi11"))-> tmp.fig.1

#ggsave(filename = "Rplots.pdf", tmp.fig.1)

##Comparison of models
# test difference LRT or anova
lrtest(DNAbyOPG, DNAbyOPG_dpi) #--> Report this table in the results 
lrtest(DNAbyOPG, DNAbyOPGxdpi) #--> Report this table in the results 
lrtest(DNAbyOPG_dpi, DNAbyOPGxdpi) #--> Report this table in the results 

###Model 4: GLMM with dpi as random factor
require(lme4)      

DNAbyOPG_dpi_glmm <- lmer(log10(Genome_copies_gFaeces)~log10(OPG) + (1|dpi),
                          data = sdt.nozero, na.action = na.exclude, REML=TRUE)

summary(DNAbyOPG_dpi_glmm)
sjPlot:: tab_model(DNAbyOPG_dpi_glmm)

### Plot estimates
## Random effect estimates
plot_model(DNAbyOPG_dpi_glmm, type = "re", show.values = TRUE)-> tmp.fig.2
#ggsave(filename = "Rplots.pdf", tmp.fig.2)

##Plot model by DPI
sdt.nozero%>%
  mutate(dpi = fct_relevel(dpi, "0","1", "2", "3", "4", "5", 
                                   "6", "7", "8", "9", "10","11"))%>%
  ggplot(aes(OPG, Genome_copies_gFaeces, fill=dpi))+
  geom_point(shape=21, size=5) +
  geom_smooth(method = lm, se=FALSE, aes(OPG, Genome_copies_gFaeces, color=dpi))+
  scale_color_manual(values = colores, guide= "none")+
  scale_fill_manual(values = colores, guide= "none")+
  scale_x_log10(name = "log10 (Oocysts/g Faeces) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces) \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  labs(tag= "c")+
  theme(text = element_text(size=16), legend.position = "none")+
  annotation_logticks()-> C

##To visualize it externally 
#ggsave(filename = "Rplots.pdf", C)

##Extraction of coefficients by dpi
sdt.nozero%>% 
  nest(-dpi)%>% 
  mutate(cor=map(data,~lm(log10(Genome_copies_gFaeces) ~ log10(OPG), data = .))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()-> x

x$data<- NULL
x$cor<- NULL
x%>%
  arrange(dpi)%>%
  filter(term != "(Intercept)")-> fitted_models_dpi

#write.csv(fitted_models_dpi, "data/Experiment_results/Quant_Eimeria/Tables/Q1_OPG_DNA_estimates_DPI.csv",  row.names = F)

sdt.nozero%>% 
  nest(-dpi)%>% 
  mutate(cor=map(data,~cor.test(log10(.x$Genome_copies_gFaeces), log10(.x$OPG), method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()-> x

x$data<- NULL
x$cor<- NULL
x%>%
  arrange(dpi)-> corOPGbyDNA_DPI

#write.csv(corOPGbyDNA_DPI, "Tables/Q1_OPG_DNA_Correlation_DPI.csv",  row.names = F)
##Non significant correlation between measurements by DPI

rm(x, tmp.fig.1, tmp.fig.2, DNAbyOPG, DNAbyOPG_dpi, DNAbyOPG_dpi_glmm, DNAbyOPGxdpi)