## upload E64_total
##re-name as my Data 

myData = E64_total

# dpi as numbers
myData$dpi = as.numeric(as.character(myData$dpi))

#  -----------------

# 1. Calculate maximum value & at which dpi it happens
# 1.1 OPG
datMaxOPG = myData %>% dplyr::group_by(EH_ID) %>%
  filter(OPG == max(OPG, na.rm = T)) %>%
  dplyr::select(EH_ID, dpi, OPG,Variable)%>% data.frame()
names(datMaxOPG)[names(datMaxOPG)%in% "dpi"] <- "dpi_maxOPG"

table(datMaxOPG$dpi)
# 1.2 fecal DNA
datMaxDNA = myData %>% group_by(EH_ID) %>%
  dplyr::filter(Genome_copies_gFaeces == max(Genome_copies_gFaeces, na.rm = T)) %>%
  dplyr::select(EH_ID, dpi, Genome_copies_gFaeces,Variable)%>% data.frame()
names(datMaxDNA)[names(datMaxDNA)%in% "dpi"] <- "dpi_maxDNA"

myData %>% group_by(EH_ID) %>%
  dplyr::summarize(max_genome_copies =  max(Genome_copies_gFaeces, na.rm = T) )

myData %>% group_by(EH_ID) %>%
  dplyr::filter(Genome_copies_gFaeces == max(Genome_copies_gFaeces, na.rm = T)) -> Save

table(datMaxDNA$dpi)
# 1.3 relative weight loss
datMaxWL = myData %>% group_by(EH_ID) %>%
  filter(weightloss == max(weightloss, na.rm = T)) %>%
  dplyr::select(EH_ID, dpi, weightloss,Variable)%>% data.frame()
names(datMaxWL)[names(datMaxWL)%in% "dpi"] <- "dpi_maxWL"

table(datMaxWL$dpi)

datMaxALL = merge(merge(datMaxOPG[c("EH_ID", "Variable","OPG", "dpi_maxOPG")], 
                        datMaxDNA[c("EH_ID", "Variable","Genome_copies_gFaeces", "dpi_maxDNA")]),
                  datMaxWL[c("EH_ID", "Variable", "weightloss","dpi_maxWL")])

##remove trouble marker mice 
datMaxALL<- datMaxALL[datMaxALL$EH_ID != "LM0179" , ]
datMaxALL<-datMaxALL[datMaxALL$EH_ID != "LM0181" , ]
datMaxALL<-datMaxALL[datMaxALL$EH_ID != "LM0244" , ]




# ----------------------------
# Q2(a): Does DNA predict a health effect on the host “overall” (e.g. maximum)?
# ----------------------------

# Can we predict weight loss by a combination of OPG and fecDNA?
library(lmtest)

modNull = lm(weightloss ~ 1, data = datMaxALL)

modFullest = lm(weightloss ~ OPG * Genome_copies_gFaeces*Variable, data = datMaxALL)
summary(modFullest) 

modFull_red1 = lm(weightloss ~ OPG * Genome_copies_gFaeces+Variable, data = datMaxALL)
summary(modFull_red1)

modFull_red2 = lm(weightloss ~ OPG + Genome_copies_gFaeces*Variable, data = datMaxALL)
summary(modFull_red2)

modFull_red3 = lm(weightloss ~ OPG * Variable + Genome_copies_gFaeces, data = datMaxALL)
summary(modFull_red3)

anova(modFullest,modFull_red1)
anova(modFullest,modFull_red2)
anova(modFullest,modFull_red3)

AIC(modFullest, modFull_red1, modFull_red2, modFull_red3)

#We go on comapring against modFUll_red3

anova(modFull_red3, modred_red_1)


modOPG_Int = lm(weightloss ~ OPG * Variable, data = datMaxALL)
summary(modOPG_Int)


modGenome_Int = lm(weightloss ~  Genome_copies_gFaeces * Variable, data = datMaxALL)
summary(modGenome_Int)


modred_red_1 = lm(weightloss ~ OPG + Variable + Genome_copies_gFaeces, data = datMaxALL)
summary(modred_red_1)
modGenome = lm(weightloss ~  Genome_copies_gFaeces + Variable, data = datMaxALL)
summary(modGenome)
modOPG = lm(weightloss ~ OPG + Variable, data = datMaxALL)
summary(modOPG)

anova(modred_red_1, modOPG)
anova(modred_red_1, modGenome)

modred_red_1$AIC <- AIC(modred_red_1)
modGenome$AIC <- AIC(modGenome)
modOPG$AIC <- AIC(modOPG)

stargazer(modred_red_1, modGenome, modOPG, type = "html", out = "weightloss_model.html",
          title="Linear Models Predicting Weight Loss",
          keep.stat=c("aic"), omit.table.layout="n")

AIC(modOPG, modGenome, modOPG_Int,modGenome_Int,modred_red_1)

#modred_red_1 is the best 

library(ggeffects)

pred <- ggpredict(modred_red_1,terms=c("Genome_copies_gFaeces","OPG","Variable"), ci.lvl = NA,add.data = TRUE)

ggplot(pred, aes(x = x, y = predicted, colour = group)) +
  stat_smooth(method = "lm", se = FALSE) +
  facet_wrap(~facet) +
  #geom_point(aes())+ #add raw data argument 
  labs(
    y = get_y_title(pred),
    x = get_x_title(pred),
    colour = get_legend_title(pred))+
  theme_bw()

plot(pred, rawdata = TRUE, dot.alpha = 0.9)
methods(plot)

help (plot.ggeffects)













