#Written by Alexa Cohn on January 16, 2020
#Last edited by Alexa Cohn on January 17, 2020

library(ggplot2)
library(tidyr)
library(dplyr)

GOterms <- read.csv("~/enriched_GOterms.csv", stringsAsFactors = FALSE, header = TRUE)
colnames(GOterms)
names(GOterms)[1] <- "GO.Term"

# subset based on type of GO term
BP <- subset(GOterms, Type == "BP", select = c(GO.Term, Cerro, Javiana, Typhimurium))
MF <- subset(GOterms, Type == "MF", select = c(GO.Term, Cerro, Javiana, Typhimurium))
CC <- subset(GOterms, Type == "CC", select = c(GO.Term, Cerro, Javiana, Typhimurium))

# rearrange data for ggplot
BP_df <- gather(BP, key = Serovar, value = counts, Cerro:Typhimurium)

MF_df <- gather(MF, key = Serovar, value = Counts, Cerro:Typhimurium)

CC_df <- gather(CC, key = Serovar, value = Counts, Cerro:Typhimurium)

# make plots of each category of GO terms

ggplot(data = BP_df, aes(x = stringr::str_wrap(GO.Term,50), y = counts, fill = Serovar))+
  geom_bar(position = "dodge", stat = "identity")+
  coord_flip()+
  labs(x = "Biological Process", y = "Counts of GO Terms")+
  scale_fill_manual(values = c("turquoise", "mediumorchid2", "lightpink1"))
ggsave("Figure2a.pdf", plot = last_plot(), width = 2.2, height = 2.9, units = "in" )

ggplot(data = MF_df, aes(x = stringr::str_wrap(GO.Term,30), y = Counts, fill = Serovar))+
  geom_bar(position = "dodge", stat = "identity")+
  coord_flip()+
  labs(x = "Molecular Function", y = "Counts of GO Terms")+
  scale_fill_manual(values = c("turquoise", "mediumorchid2", "lightpink1"))

ggplot(data = CC_df, aes(x = stringr::str_wrap(GO.Term,30), y = Counts, fill = Serovar))+
  geom_bar(position = "dodge", stat = "identity")+
  coord_flip()+
  labs(x = "Cellular Component", y = "Counts of GO Terms")+
  scale_fill_manual(values = c("turquoise", "mediumorchid2", "lightpink1"))
