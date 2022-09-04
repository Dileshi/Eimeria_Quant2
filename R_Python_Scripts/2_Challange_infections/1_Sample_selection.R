library(dplyr)
library(tidyverse)

data = read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")

challenge <- data[data$EH_ID %in% c("LM0180", "LM0182","LM0184","LM0185","LM0186","LM0195","LM0228","LM0229","LM0238", "LM0191",
                            "LM0240","LM0246","LM0247","LM0194","LM0190","LM0199","LM0197","LM0198","LM0188","LM0254","LM0192"),]
challenge <- challenge %>% 
  dplyr::filter(infection == "challenge")

write.csv(challenge,"Quant2_challenge.csv")






