library(dplyr)

data_Fig2 <- read.csv("Untitled/Sheet 2-Table 1.csv")

data_Fig2$No..Of.chromosome = as.numeric(data_Fig2$No..Of.chromosome)

data_Fig2%>%
  ggplot(aes(x=Days, y=No..Of.chromosome))+
  geom_point(color="tomato", size=3)












