#rename the entries under the column 'Sample_ID' in file Nanodrop_all to 
#E57aXXX to allow merge with file E88_Primary

Nanodrop_all$labels = paste('E57a',Nanodrop_all$Sample_ID, sep = '')

#merge files E88_Primary and Nanodrop_all
library(dplyr)
Merged = left_join(E88_Primary,Nanodrop_all)

#remove column 'Sample_ID'
DNA_Extractions = subset(Merged, select = -c(Sample_ID) )

#export DNA_Extractions file
write.csv(DNA_Extractions, 'DNA_Extractions.csv')