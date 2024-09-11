# clear space
rm(list=ls()) 
# clear graphics
graphics.off() 

library(dplyr)
library(tidyr)
library(readxl) 
library(tidyverse) 
library(here)
library(AMR)
library(reshape2)
library(ggplot2)


df1 <- read_excel("heatmap_data3.xlsx")
df1$Class <- as.factor(df1$Class)
str(df1)

#reorder antibiotics
#df1$Antibiotic <- factor(df1$Antibiotic, levels=c("Amikacin", "Gentamicin",
#                                          "Imipenem", "Meropenem",
#                                          "Cefepime","Ceftazidime",
#                                          "Ceftriaxone","Ciprofloxacin",
#                                          "Levofloxacin","Ampicillin",
#                                          "Amoxicillin/clavulanic acid",
#                                          "Piperacillin/tazobactam",
#                                          "Colistin", "Trimethoprim/sulfamethoxazole"))


ggp <- ggplot(df1, aes(x=Source, y=Class, fill = tau)) +                           
  geom_tile(color = "white")+
  scale_fill_gradient2(low = 'cornsilk1',  mid = 'gold', high = 'forestgreen', #midpoint = 0.5, 
                       name = 'Kendall correlation\n coefficient tau') + 
  theme_classic() +
  theme(legend.position="right") +
  xlab("Sampling sources") +
  ylab("Antibiotic class") +
  facet_wrap(~Pathogen) +
  geom_text(aes(label = round((Mean_resistance), digits=2)), 
            color = "white", fontface='bold', size = 6) +
  theme(strip.text.x = element_text(size = 12, face = "italic"),
        axis.text.x = element_text( size = 12, face = "bold"),
        axis.text.y = element_text( size = 12, face = "bold"),
        axis.title = element_text( size = 16, face = "bold"),
        legend.text=element_text(size=10))

ggp


ggsave("heatmap3_tau.jpeg", width = 40, height = 20, units = "cm")

