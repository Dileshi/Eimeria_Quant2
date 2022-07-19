##Script to idetify stage specific genes and a house keeping gene for RNA primer design
library(hacksaw)
library(pheatmap)
library(tidyverse)
library(dplyr)

##Load processed sequncing data fom project Ehret et al., 2017 and Prep
data <- read.csv("https://raw.githubusercontent.com/derele/Ef_RNAseq/master/output_data/Ef_norm_counts.csv", header=FALSE,na.strings=c("NA","NaN", " ", "?"))
#shift fits row to the right                 
data <- data %>% shift_row_values(at = 1, .dir = "right")
#Assign first row as column headers
colnames(data)<-data[c(1),]
data<-data[-c(1),]
##Give column the heading 'Gene.ID'
names(data)[1] <- paste("Gene.ID")



#####Firstly, Begin with designing 2 primers for oocyst 
##Load  E. falciformis oocyst wall genes from ToxoDB
OW_genes <- read.csv("Raw_Data/Primer_Design/Genes_Oocyst_Wall.csv")

##create subset in data' that list oocyst wall genes (OW-genes) 
OW_genes_vec <- OW_genes$Gene.ID
data_OW <-data[data$Gene.ID %in% OW_genes_vec,]

##add COX1 and 60s genes to compare expressions
add <- data[data$Gene.ID %in% c("EfaB_PLUS_15576.g1276", "EfaB_MINUS_39892.g2639","EfaB_PLUS_25458.g2106"), ] 
data_OW <-bind_rows(data_OW,anti_join(add,data_OW,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in OW_genes must be converted back to rows names
data_OW <- data_OW %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
#rename COX1 genes 
row.names(data_OW)[51] <- 'Cytochrome_c_oxidase_subunit_2'
row.names(data_OW)[52] <- 'Cytochrome_c_oxidase_subunit'
row.names(data_OW)[53] <- '60S ribosomal protein L16'

##make columns numeric 
data_OW <- data_OW %>% mutate_if(is.character,as.numeric)

##As Oocyst shedding begins after dpi5 (see Ehret et al., 2017  figure 1), a few mice pre-dpi5 were excluded to better aid visualisation 
data_OW <- data_OW %>% select(-c('C57BL6_1stInf_5dpi_rep2','NMRI_1stInf_5dpi_rep3',
                                 'NMRI_1stInf_3dpi_rep2', 'NMRI_1stInf_5dpi_rep2','NMRI_1stInf_5dpi_rep1',
                                 'NMRI_2ndInf_3dpi_rep2','NMRI_2ndInf_5dpi_rep1','Rag_2ndInf_5dpi_rep1','Rag_1stInf_5dpi_rep2', 
                                 'Rag_1stInf_5dpi_rep1'))
##plot pheatmap
pheatmap(data_OW[1:53,],
         display_numbers = TRUE,
         cutree_rows = 3,)
rm(data_OW,OW_genes,OW_genes_vec,add)

##Note expression level to help distinguish between house keeping gene
##compare OW expression with genes house keeping genes and C57BL/6 control 
###DECIDE which gene to use based on expression level on map:EfaB_MINUS_1001.g115. , EfaB_PLUS_2387.g293 , EfaB_MINUS_40637.g2673
 ##exclude gene  EfaB_MINUS_6743.g610 with 39000 expression (see Excel sheet 2 for reasoning and details)

###Plot with over-expressed genes removed###
#OW_genes <- read.csv("Raw_Data/Primer_Design/Genes_Oocyst_Wall.csv")
#OW_genes <- OW_genes[!(OW_genes$Gene.ID=="EfaB_MINUS_6743.g610" | OW_genes$Gene.ID=="EfaB_MINUS_1001.g115" | OW_genes$Gene.ID== 'EfaB_PLUS_2387.g293' | OW_genes$Gene.ID== 'EfaB_PLUS_2074.g242'),]
#OW_genes_vec <- OW_genes$Gene.ID
#data_OW <-data[data$Gene.ID %in% OW_genes_vec,]
#add <- data[data$Gene.ID %in% c("EfaB_PLUS_15576.g1276", "EfaB_MINUS_39892.g2639","EfaB_PLUS_25458.g2106"), ] 
#data_OW <-bind_rows(data_OW,anti_join(add,data_OW,by="Gene.ID"))
##Now, to plot a pheatmap, the Gene.ID column in OW_genes must be converted back to rows names
#data_OW <- data_OW %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
#rename COX1 genes 
#row.names(data_OW)[47] <- 'Cytochrome_c_oxidase_subunit_2'
#row.names(data_OW)[48] <- 'Cytochrome_c_oxidase_subunit'
#row.names(data_OW)[49] <- '60S ribosomal protein L16'
##make columns numeric 
#data_OW <- data_OW %>% mutate_if(is.character,as.numeric)
##plot pheatmap
#pheatmap(data_OW[1:49,],
         #display_numbers = TRUE,
         #cutree_rows = 2,)
#rm(data_OW,OW_genes,OW_genes_vec,add)




#####Secondly, design 2 primers for microgametes 
##Load  E. falciformis 'flagella' and 'flagellar' genes from ToxoDB
Flagella_genes <- read.csv("Raw_Data/Primer_Design/Genes_Flagella_Microgamete.csv")
Flagellar_genes <- read.csv("Raw_Data/Primer_Design/Genes_Flagellar_Microgamete.csv") 

##create subset in 'data' that list flagella genes  
Flagella_genes_vec <- Flagella_genes$Gene.ID
Flagellar_genes_vec <- Flagellar_genes$Gene.ID
data_flagella <-data[data$Gene.ID %in% c(Flagella_genes_vec, Flagellar_genes_vec),]

##add actin and 18s ribosomal RNA genes to compare expressions
add <- data[data$Gene.ID %in% c("EfaB_PLUS_4114.g414", "EfaB_PLUS_31498.g2326"), ] 
data_flagella <-bind_rows(data_flagella,anti_join(add,data_flagella,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in flagella_genes must be converted back to rows names
data_flagella <- data_flagella %>% remove_rownames %>% column_to_rownames(var="Gene.ID")

#rename COX1 genes 
row.names(data_flagella)[57] <- 'actin'
row.names(data_flagella)[56] <- '18S_ribosomal_RNA'

##make columns numeric 
data_flagella <- data_flagella %>% mutate_if(is.character,as.numeric)

##oocyst shed peak after dpi7. Gametes must form after dpi3 -"sexual gametocytes can be observed from 4 to 6 dpi" (Jarquín-Díaz et al., 2022)
##hide a few dpi3 and oocyst expressions for better visulaisation 
#data_flagella <- data_flagella %>% select(-c('NMRI_1stInf_3dpi_rep1',
                                 #'NMRI_2ndInf_3dpi_rep2','NMRI_oocysts_rep2'))
##plot pheatmap
pheatmap(data_flagella[1:57,],
         display_numbers = TRUE,
         cutree_rows = 3,
         cutree_cols = 3,
         treeheight_col=0,)
rm(data_flagella, Flagella_genes,Flagella_genes_vec,Flagellar_genes,Flagellar_genes_vec,add)

##next microgamete candidate search term "armadillo". According to Walker et al., 2015 this genes was 
 #highly expressed in microgametes and is homolog to the plasmodium falciparum PF16 - a component of the 
 #axoneme central apparatus and is essential for microgamete motility and fertilisation. 
##Load  E. falciformis 'armadillo' genes from ToxoDB
armadillo_genes <- read.csv("Raw_Data/Primer_Design/Genes_microgametes_armadillo.csv")

##create subset in 'data_armadillo' that list  genes  
armadillo_genes_vec <- armadillo_genes$Gene.ID
data_armadillo <-data[data$Gene.ID %in% armadillo_genes_vec,]

##add COX gene to compare expressions
#add <- data[data$Gene.ID %in% "EfaB_PLUS_15576.g1276", ] 
#data_armadillo <-bind_rows(data_armadillo,anti_join(add,data_armadillo,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in armadillo_genes must be converted back to rows names
data_armadillo <- data_armadillo %>% remove_rownames %>% column_to_rownames(var="Gene.ID")

##remove genes with no change in expression avoid pheatmap overcrowing 
remove_rownames <- c("EfaB_PLUS_25177.g2048", "EfaB_PLUS_4849.g448", "EfaB_PLUS_3162.g328", "EfaB_PLUS_7048.g689", 
                     "EfaB_PLUS_20371.g1649", "EfaB_PLUS_24701.g2018", "EfaB_PLUS_5767.g515", "EfaB_PLUS_7048.g696", 
                     "EfaB_PLUS_3469.g356", "EfaB_PLUS_3469.g372", "EfaB_MINUS_28105.g2313", "EfaB_PLUS_5670.g138", 
                     "EfaB_MINUS_1068.g138", "EfaB_MINUS_17400.g1575", "EfaB_MINUS_33184.g2507", "EfaB_MINUS_43533.g2753", 
                     "EfaB_MINUS_31498.g2393", "EfaB_MINUS_32574.g2352", "EfaB_MINUS_25052.g2042","EfaB_7742.g762",
                     "EfaB_MINUS_22450.g1934", "EfaB_MINUS_5504.g463","EfaB_MINUS_7048.g842","EfaB_PLUS_5670.g495",
                     "Efab_PLUS_4849.g445","EfaB_PLUS_15697.g1349","EfaB_PLUS_14881.g1217","EfaB_PLUS_720.g48",
                     "EfaB_MINUS_31746.g2425","EfaB_MINUS_7742.g762")
data_armadillo <- data_armadillo[!(row.names(data_armadillo) %in% remove_rownames), ]

#rename COX1 genes 
#row.names(data_armadillo)[110] <- 'Cytochrome_c_oxidase_subunit'

##make columns numeric 
data_armadillo <- data_armadillo %>% mutate_if(is.character,as.numeric)

##oocyst shed peak after dpi7. Gametes must form after dpi3 -"sexual gametocytes can be observed from 4 to 6 dpi" (Jarquín-Díaz et al., 2022)
##hide a few dpi3 and oocyst expressions for better visulaisation 
#data_armadillo <- data_armadillo %>% select(-c('NMRI_1stInf_3dpi_rep1','NMRI_1stInf_3dpi_rep2', 'NMRI_oocysts_rep1',
                                             #'NMRI_2ndInf_3dpi_rep2','NMRI_oocysts_rep2','NMRI_1stInf_3dpi_rep1'))
##plot pheatmap
pheatmap(data_armadillo[1:85,],
         display_numbers = TRUE,
         cutree_rows = 2,
         cutree_cols = 2,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=6,
         fontsize_row = 6,)
rm(data_armadillo, armadillo_genes,armadillo_genes_vec,armadillo_genes,armadillo_genes_vec,add,remove_rownames)

##next microgamete candidate was found with the search term "female gametete".  
#the results showed a male-gamete-fusion protein refer Walker et al., 2015
gam_fusion_genes <- read.csv("Raw_Data/Primer_Design/Genes_male_fusion.csv")

##create subset in data' that list gam_fusion genes  
gam_fusion_genes_vec <- gam_fusion_genes$Gene.ID
data_gam_fusion <-data[data$Gene.ID %in% gam_fusion_genes_vec,]

##add COX1 and 60s genes to compare expressions
add2 <- data[data$Gene.ID %in% c("EfaB_PLUS_5788.g526", "EfaB_MINUS_39892.g2639","EfaB_PLUS_25458.g2106","EfaB_PLUS_15576.g1276"), ] 
data_gam_fusion <-bind_rows(data_gam_fusion,anti_join(add2,data_gam_fusion,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in gam_fusion_genes must be converted back to rows names
data_gam_fusion <- data_gam_fusion %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
#rename COX1 genes 
row.names(data_gam_fusion)[8] <- 'cytochrome C oxidase subunit_2'
row.names(data_gam_fusion)[9] <- 'ribosomal protein RPL13'
row.names(data_gam_fusion)[10] <- 'cytochrome c oxidase subunit'
row.names(data_gam_fusion)[11] <- 'ribosomal RNA large subunit '

##make columns numeric 
data_gam_fusion <- data_gam_fusion %>% mutate_if(is.character,as.numeric)

##plot pheatmap
pheatmap(data_gam_fusion [1:11,],
         display_numbers = TRUE,
         cutree_rows = 5,
         cutree_cols = 5)
rm(data_gam_fusion,gam_fusion_genes,gam_fusion_genes_vec,add2)
#In our case, oocyst shedding begins at dpi7 for 8 mice (dpi0 excluded)and peak is at dpi8 for all mice. 
##armadillo gene is not successful as gene expression variance across dpi - hence avoid and focus on flagellar genes
###DECIDE which gene to use based on expression level on map: refer Excel sheet 2


#####Thirdly, design 2 primers for macrogametes (serach term: egg cell (Failed:macrogamete)) 
##Load  E. falciformis macrogamete genes from ToxoDB
macro_genes <- read.csv("Raw_Data/Primer_Design/Genes_macrogametes.csv")

##create subset in data' that list  macrogamete genes  
macro_genes_vec <- macro_genes$Gene.ID
data_macro <-data[data$Gene.ID %in% macro_genes_vec,]

##add actin, C-oxidase and 18s genes to compare expressions
add <- data[data$Gene.ID %in% c("EfaB_PLUS_4114.g414", "EfaB_PLUS_31498.g2326","EfaB_MINUS_39892.g2639"), ] 
data_macro <-bind_rows(data_macro,anti_join(add,data_macro,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in OW_genes must be converted back to rows names
data_macro <- data_macro %>% remove_rownames %>% column_to_rownames(var="Gene.ID")

#rename reference genes 
row.names(data_macro)[315] <- 'actin'
row.names(data_macro)[313] <- '18S_ribosomal_RNA'
row.names(data_macro)[314] <- 'cytochrome_C_oxi'

##make columns numeric 
data_macro<- data_macro %>% mutate_if(is.character,as.numeric)

##oocyst shed peak after dpi7. Gametes must form after dpi3 -"sexual gametocytes can be observed from 4 to 6 dpi" (Jarquín-Díaz et al., 2022)
##hide dpi3 and oocyst and sporoziote expression to calculate sum of expression 
data_macro2 <- data_macro %>% select(-c('NMRI_1stInf_3dpi_rep1','NMRI_1stInf_3dpi_rep2', 'NMRI_oocysts_rep1', 'NMRI_sporozoites_rep1',
                                        'NMRI_2ndInf_3dpi_rep2','NMRI_oocysts_rep2','NMRI_1stInf_3dpi_rep1', 'NMRI_sporozoites_rep2',))
data_macro2$sum <- rowSums(data_macro2)

#select values with sum>2500 (mean value of ref. genes actin, 18s and C-Oxi) to better visualize data 
data_macro2 <- filter(data_macro2,(sum>2500))

#list selected genes and to subsut in data_macro
list_macro2 <- rownames(data_macro2)
data_macro <-data_macro[rownames(data_macro) %in% list_macro2,]

##plot pheatmap
pheatmap(data_macro[1:87,],
         display_numbers = TRUE,
         cutree_rows = 4,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=7,
         fontsize_row = 5,
         )
rm(data_macro,macro_genes,macro_genes_vec,add,data_macro2,list_macro2)
##NO FRUITFUL genes found Search term 'egg cell' failed.
#Gene 'EfaB_PLUS_22284.g1821' highly expressed speifically at dpi7 but orthology doesnot specift it to be a 
 #macrogamate stage specific gene (gcc3 in Toxoplasma and Plectin in other species )

##search term in toxoDB: "amiloride-sensitive amine" One interesting transcript upregulated in 
#E. tenella macrogametes is amiloride-sensitive amine oxidase (Martorelli Di Genova et al., 2020). 
keep_rows <- c("EfaB_PLUS_3469.g355", "EfaB_PLUS_4114.g414", "EfaB_MINUS_39892.g2639") 
enzy <- data[rownames(data) %in% keep_rows, ]
row.names(enzy)[row.names(enzy) == "EfaB_PLUS_4114.g414"] <- "actin"
row.names(enzy)[row.names(enzy) == "EfaB_MINUS_39892.g2639"] <- "cytochro_C_Oxi"
pheatmap(enzy[1:3,],display_numbers = TRUE,)
rm(enzy,keep_rows)
#EfaB_PLUS_3469.g355 a potential candidate







##next sporozoite candidate search term "sporozoite". 
##Load  E. falciformis 'sporozoite' genes from ToxoDB
sporo_genes <- read.csv("Raw_Data/Primer_Design/Genes_Sporozoite.csv")

##create subset in 'data_sporo' that list  genes  
sporo_genes_vec <- sporo_genes$Gene.ID
data_sporo <-data[data$Gene.ID %in% sporo_genes_vec,]

##add COX1 and 60s genes to compare expressions
add <- data[data$Gene.ID %in% c("EfaB_PLUS_15576.g1276", "EfaB_MINUS_39892.g2639","EfaB_PLUS_25458.g2106"), ] 
data_sporo <-bind_rows(data_sporo,anti_join(add,data_sporo,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in OW_genes must be converted back to rows names
data_sporo <- data_sporo %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
#rename COX1 genes 
row.names(data_sporo)[8] <- 'Cytochrome_c_oxidase_subunit_2'
row.names(data_sporo)[10] <- 'Cytochrome_c_oxidase_subunit'
row.names(data_sporo)[9] <- '60S ribosomal protein L16'

##make columns numeric 
data_sporo <- data_sporo %>% mutate_if(is.character,as.numeric)

##plot pheatmap
pheatmap(data_sporo[1:10,],
         display_numbers = TRUE,
         #cutree_rows = 3,
         )
rm(data_sporo,sporo_genes,sporo_genes_vec,add)
##UNSUCESSFUL: Gene EfaB_PLUS_21840.g1792 is a rhoptry related protein but according to heatmap expressed 
#at oocyst stages and not in sporozoite stages.  


#next try involes a 'rhoptry protein' <- search term based on Olajide et al., 2022 
##Load  E. falciformis 'sporozoite' genes from ToxoDB
ROP_genes <- read.csv("Raw_Data/Primer_Design/Genes_Rhoptry_Sprozoites.csv")

##create subset in 'data_sporo' that list  genes  
ROP_genes_vec <- ROP_genes$Gene.ID
data_ROP <-data[data$Gene.ID %in% ROP_genes_vec,]

##add actin, C-oxidase and 18s genes to compare expressions
add <- data[data$Gene.ID %in% c("EfaB_PLUS_4114.g414","EfaB_MINUS_39892.g2639","EfaB_PLUS_28105.g2216"), ] 
data_ROP <-bind_rows(data_ROP,anti_join(add,data_ROP,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in ROP_genes must be converted back to rows names
data_ROP <- data_ROP %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
#change reference genes names
row.names(data_ROP)[row.names(data_ROP) == "EfaB_PLUS_4114.g414"] <- "actin"
row.names(data_ROP)[row.names(data_ROP) == "EfaB_MINUS_39892.g2639"] <- "cytochro_C_Oxi"
row.names(data_ROP)[row.names(data_ROP) == "EfaB_PLUS_28105.g2216"] <- "60S_ribosome"

##make columns numeric 
data_ROP <- data_ROP %>% mutate_if(is.character,as.numeric)

##select sporozites to identify genes with high ROP expression 
data_ROPsub <- data_ROP %>% select(c('NMRI_sporozoites_rep1', 'NMRI_sporozoites_rep2',))
data_ROPsub$sum <- rowSums(data_ROPsub)
#sum for actin=2.667624e+02 and cytochrome oxidase=5053 and60s ribosome=3049 
#select values with sum>400 
data_ROPsub <- filter(data_ROPsub,(sum>1100))

#list selected genes and to subsut in data_macro
list_ROPsub <- rownames(data_ROPsub)
data_ROP <-data_ROP[rownames(data_ROP) %in% list_ROPsub,]

#exclude oocyst genes that are over expressed/ i.e. list all genes with oocyst expression under sum<900
data_ROPoocy <- data_ROP %>% select(c('NMRI_oocysts_rep1', 'NMRI_oocysts_rep2',))
data_ROPoocy$sum <- rowSums(data_ROPoocy)
data_ROPoocy <- filter(data_ROPoocy,(sum<600))
list_ROPoocy <- rownames(data_ROPoocy)
data_ROP <- data_ROP[rownames(data_ROP) %in% list_ROPoocy,]

#filter data_ROP for genes overexpressed in other stages
#keep genes of dpi 7 that are expressed under 1900 and so on......
data_ROP <- filter(data_ROP,(data_ROP$NMRI_1stInf_7dpi_rep1 <1900))
data_ROP <- filter(data_ROP,(data_ROP$NMRI_1stInf_5dpi_rep2 <1900))
data_ROP <- filter(data_ROP,(data_ROP$NMRI_1stInf_7dpi_rep2 <1900))
data_ROP <- filter(data_ROP,(data_ROP$NMRI_2ndInf_5dpi_rep1 <1900))
data_ROP <- filter(data_ROP,(data_ROP$Rag_2ndInf_5dpi_rep1 <1900))
data_ROP <- filter(data_ROP,(data_ROP$Rag_1stInf_5dpi_rep2 <1900))
data_ROP <- filter(data_ROP,(data_ROP$C57BL6_1stInf_5dpi_rep1 <1900))

##plot pheatmap
pheatmap(data_ROP[1:201,],
         display_numbers = TRUE,
         fontsize_number = 4,
         cutree_rows = 8,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=7,
         fontsize_row = 5,
         )
rm(add,data_ROP, data_ROPoocy,data_ROPsub, list_ROPoocy,list_ROPsub,ROP_genes_vec,ROP_genes)

##Just check which genes are highly expressed in sprozoites 
data <- data %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
data <- data %>% mutate_if(is.character,as.numeric)
data_Spo_high <- filter(data,(data$NMRI_sporozoites_rep2 >900 | data$NMRI_sporozoites_rep1>900 ))

pheatmap(data_Spo_high[1:189,],
         display_numbers = TRUE,
         fontsize_number = 4,
         cutree_rows = 8,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=7,
         fontsize_row = 5,
         )
rm(data_Spo_high)
#nothing additional found 





###Look for 'SAG' and 'Surface antigens' as search terms on ToxoDB => meroziotes
##Load  E. falciformis genes from ToxoDB
SAG_genes <- read.csv("Raw_Data/Primer_Design/Genes_SAG.csv")
SAG2_genes<- read.csv("Raw_Data/Primer_Design/Genes_Surface_Antigen.csv")

##create subset in 'data' that list SAG genes  
SAG_genes_vec <- SAG_genes$Gene.ID
SAG2_genes_vec <- SAG2_genes$Gene.ID
data_SAG <-data[data$Gene.ID %in% c(SAG_genes_vec, SAG2_genes_vec),]

##add actin, C-oxidase and 18s genes to compare expressions
add <- data[data$Gene.ID %in% c("EfaB_PLUS_4114.g414","EfaB_MINUS_39892.g2639","EfaB_PLUS_28105.g2216"), ] 
data_SAG <-bind_rows(data_SAG,anti_join(add,data_SAG,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in SAG_genes must be converted back to rows names
data_SAG <- data_SAG %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
#change reference genes names
row.names(data_SAG)[row.names(data_SAG) == "EfaB_PLUS_4114.g414"] <- "actin"
row.names(data_SAG)[row.names(data_SAG) == "EfaB_MINUS_39892.g2639"] <- "cytochro_C_Oxi"
row.names(data_SAG)[row.names(data_SAG) == "EfaB_PLUS_28105.g2216"] <- "60S_ribosome"

##make columns numeric 
data_SAG <- data_SAG %>% mutate_if(is.character,as.numeric)

#filter data_SAG for genes overexpressed in other stages and over over expressed genes for dpi3
data_SAG <- filter(data_SAG,(data_SAG$NMRI_1stInf_3dpi_rep2 >100))
data_SAG <-filter(data_SAG,(data_SAG$NMRI_sporozoites_rep2<2000))
data_SAG <-filter(data_SAG,(data_SAG$NMRI_sporozoites_rep1<2000))

##plot pheatmap
pheatmap(data_SAG[1:240,],
         display_numbers = TRUE,
         fontsize_number = 4,
         cutree_rows = 6,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=7,
         fontsize_row = 5,
         )
rm(data_SAG,SAG_genes,SAG_genes_vec,add,SAG2_genes,SAG2_genes_vec)
#identify a SAG gene not expressed in sporozoites but within dpi 3 and/or 5 (Fig.4 C Marugan-Hernandez et al., 2020)
#such a SAG protein most likely belongs to merozoites 
#check SAG_genes of sprozoites:
data_SAG <- filter(data_SAG,(data_SAG$NMRI_sporozoites_rep1 >400))
pheatmap(data_SAG[1:48,])
##no others found

###Look for 'MIC' and 'microneme' as search terms on ToxoDB => meroziotes
##Load  E. falciformis genes from ToxoDB
MIC_genes <- read.csv("Raw_Data/Primer_Design/Genes_Microneme.csv")

##create subset in data' that list oocyst wall genes (OW-genes) 
MIC_genes_vec <- MIC_genes$Gene.ID
data_MIC <-data[data$Gene.ID %in% MIC_genes_vec,]

##add actin, C-oxidase and 18s genes to compare expressions
add <- data[data$Gene.ID %in% c("EfaB_PLUS_4114.g414","EfaB_MINUS_39892.g2639","EfaB_PLUS_28105.g2216"), ] 
data_MIC <-bind_rows(data_MIC,anti_join(add,data_MIC,by="Gene.ID"))

##Now, to plot a pheatmap, the Gene.ID column in MIC_genes must be converted back to rows names
data_MIC <- data_MIC %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
#change reference genes names
row.names(data_MIC)[row.names(data_MIC) == "EfaB_PLUS_4114.g414"] <- "actin"
row.names(data_MIC)[row.names(data_MIC) == "EfaB_MINUS_39892.g2639"] <- "cytochro_C_Oxi"
row.names(data_MIC)[row.names(data_MIC) == "EfaB_PLUS_28105.g2216"] <- "60S_ribosome"
##make columns numeric 
data_MIC <- data_MIC %>% mutate_if(is.character,as.numeric)

##plot pheatmap
pheatmap(data_MIC[1:23,],
         display_numbers = TRUE,
         cutree_rows = 3,
         cutree_cols = 3)
rm(data_MIC,MIC_genes,MIC_genes_vec,add)




###Look for highly expressed genes for oocyst 
data <- data %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
data <- data %>% mutate_if(is.character,as.numeric)
data_Oocy_high <- filter(data,(data$NMRI_oocysts_rep1 >1500 | data$NMRI_oocysts_rep2 >1500 ))

pheatmap(data_Oocy_high[1:168,],
         display_numbers = TRUE,
         fontsize_number = 4,
         cutree_rows = 8,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=7,
         fontsize_row = 5,)
rm(data_Oocy_high)

#Look for highly expressed genes at dpi 3
data_dpi3_high <- filter(data,(data$NMRI_1stInf_3dpi_rep1 >1500 | data$NMRI_1stInf_3dpi_rep2 >1500 |
                               data$NMRI_2ndInf_3dpi_rep2>1000))
pheatmap(data_dpi3_high[1:74,],
         display_numbers = TRUE,
         fontsize_number = 4,
         cutree_rows = 6,
         cutree_cols = 5,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=7,
         fontsize_row = 5,)
rm(data_dpi3_high)

#Look for highly expressed genes at dpi 5
data_dpi5_high <- filter(data,(data$NMRI_1stInf_5dpi_rep1>1400 | data$NMRI_1stInf_5dpi_rep2 >1400 |
                                 data$NMRI_1stInf_5dpi_rep3>1400 | data$NMRI_2ndInf_5dpi_rep1>1400))

pheatmap(data_dpi5_high[1:75,],
         display_numbers = TRUE,
         fontsize_number = 4,
         cutree_rows = 6,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=7,
         fontsize_row = 5,)
rm(data_dpi3_high)

#Look for highly expressed genes at dpi 7
data_dpi7_high <- filter(data,(data$NMRI_1stInf_7dpi_rep1>1000 | data$NMRI_1stInf_7dpi_rep2 >1000 |
                                 data$NMRI_2ndInf_7dpi_rep1>1000 | data$NMRI_2ndInf_7dpi_rep2>1000))
pheatmap(data_dpi7_high[1:172,],
         display_numbers = TRUE,
         fontsize_number = 4,
         cutree_rows = 6,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=7,
         fontsize_row = 5,)
rm(data_dpi7_high)
