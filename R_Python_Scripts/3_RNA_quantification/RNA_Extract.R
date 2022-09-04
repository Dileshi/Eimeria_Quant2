library(dplyr)
library(tidyverse)

RNA = read.csv("Output_Data/All_parameters_merged.csv")

RNA_Avaible = RNA[RNA$Faeces_for_RNA.Extraction == 'Available',] 

RNA_Avaible <- RNA %>% 
  dplyr::filter(Faeces_for_RNA.Extraction == "Available")


RNA_Ava <- RNA_Avaible[RNA_Avaible$EH_ID %in% c("LM0180", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229","LM0238", "LM0191",
                            "LM0240","LM0246","LM0247","LM0194","LM0190","LM0199","LM0197","LM0198","LM0188","LM0254","LM0192"),]

RNA_Ava <- RNA_Ava[!duplicated(RNA_Ava$labels), ]

CQW = RNA_Avaible %>% filter(labels == 'E57aCQW')


