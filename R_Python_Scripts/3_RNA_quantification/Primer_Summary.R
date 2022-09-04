#Run from root of REPO Eimeria_Quant2
##Script to summarize stage specific genes and a house keeping gene used RNA primer design
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


##Catogorise genes 
 # oocyst genes (EXPRESSED ONLY IN OOCYST)
 # early oocyst genes (EXPRESSED only/mostly IN OOCYST & dpi7)
 # merozoites (HIGH EXPRESSION AT dpi5 and not expressed in sporozoites, oocyst and dpi7)
 #sporoites 
 #macrogamete (HIGH EXPRESSION at dpi7 and stage specific)
 #microgametes (flagella genes are specific to male gametes & has HIGH EXPRESSION at dpi7)
 #reference/house keeping genes actin ect.. 

#load selected genes
genes <- read.csv("Data/Primary_Clustering_Summary.csv")
#list genes
genes_vec <- genes$Gene.ID
Summary <-data[data$Gene.ID %in% genes_vec,]

#rename genes based using genes column 2 
Summary <- merge(Summary,genes,by="Gene.ID")
Summary = subset(Summary, select = -c(Gene.ID))

#convert protein column to rows
rownames(Summary) <- Summary$Protein

#remove row protein and make dataframe numeric 
Summary = subset(Summary, select = -c(Protein,cluster,stage))
Summary <- Summary %>% mutate_if(is.character,as.numeric)

#plot pheatmap
pheatmap(Summary[1:33,],
         display_numbers = TRUE,
         fontsize_number = 5.5,
         cutree_rows = 8,
         cutree_cols = 8,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=9,
         fontsize_row = 9,
         )




##RECHECK AND RESELECT USINF FIGURE 4 of EHRET et a. 2017 gene cluster 
##if selected genes are not in gene cluster -> avoid 

cluster <- read.csv("https://raw.githubusercontent.com/derele/Ef_RNAseq/master/output_data/Ef_hclustered_cycle.csv", header=FALSE)

#make first row header of data.frame cluster
names(cluster) <- cluster[1,]
cluster <- cluster[-1,]

#load Primer summary data
Primer_summary <- read.csv("Data/Primer_Summary.csv")

#select 30 genes from cluster using lsit to subset
list_cluster <- Primer_summary$Gene.ID
cluster <-cluster[cluster$Cluster %in% list_cluster,]

#merge
cluster <- rename(cluster, Gene.ID=Cluster)
Primer_summary <- merge(Primer_summary, cluster, by="Gene.ID", all.x=TRUE)


#Ehret et al., 2017: Cluster 1: oocyst 
#                  : Cluster 5: oocyst and late gamtetocyte 
#                  : Cluster 4: sporozoites 
#                  : Cluster 2: gametes and early oocyst 
#                  : Cluster 7: gametes and early oocyst
#                  : Cluster 3: merozoites
#                  : Cluster 6: merozoites
Primer_summary <- rename(Primer_summary, "cluster"="NA")
Primer_summary$stage <- ifelse(Primer_summary$cluster==3, "merozoites",
                               ifelse(Primer_summary$cluster==6, "merozoites",
                                      ifelse(Primer_summary$cluster==2, "gametes&EarlyOocyst",
                                             ifelse(Primer_summary$cluster==7, "gametes&EarlyOocyst",
                                                    ifelse(Primer_summary$cluster==4, "sporoziotes",
                                                           ifelse(Primer_summary$cluster==1, "oocyst",
                                                                  ifelse(Primer_summary$cluster==5, "oocyst&LateGametocytes",
                                                                         "None")))))))



##filter genes from Ehret et al., cluster 1 genes 
cluster_1 <- filter(cluster, cluster==1)

vec <- cluster_1$Gene.ID
data2 <-data[data$Gene.ID %in% vec,]

rownames(data2) <- data2$Gene.ID
data2 = subset(data2, select = -c(Gene.ID))
data2 <- data2 %>% mutate_if(is.character,as.numeric)

remove_rownames <- c("EfaB_rRNA_SSU","EfaB_rRNA_LSU")
data2 <- data2[!(row.names(data2) %in% remove_rownames), ]

pheatmap(data2[1:564,],
         #display_numbers = TRUE,
         fontsize_number = 5.5,
         cutree_rows = 8,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=3,
         fontsize_row = 3,)
##use genes EfaB_PLUS_9917 (Oocy_hypothetical protein) , EfaB_PLUS_34010.g2454 (transmembrane protein)




##filter genes from Ehret et al., cluster 4 genes 
cluster_1 <- filter(cluster, cluster==4)

vec <- cluster_1$Gene.ID
data2 <-data[data$Gene.ID %in% vec,]

rownames(data2) <- data2$Gene.ID
data2 = subset(data2, select = -c(Gene.ID))
data2 <- data2 %>% mutate_if(is.character,as.numeric)

remove_rownames <- c("EfaB_rRNA_SSU","EfaB_rRNA_LSU")
data2 <- data2[!(row.names(data2) %in% remove_rownames), ]

pheatmap(data2[1:731,],
         #display_numbers = TRUE,
         fontsize_number = 5.5,
         cutree_rows = 8,
         cutree_cols = 4,
         treeheight_col=0,
         treeheight_row = 0,
         fontsize_col=3,
         fontsize_row = 3,)


