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


df2 <- read_excel("heatmap_data_atb_class.xlsx")
df2$Class <- as.factor(df2$Class)
str(df2)

df1 <- df2 %>% 
  filter(Source != "Blood") # filter out blood source from mapping tau

ggp <- ggplot(df1, aes(x=Source, y=Class, fill = tau)) +                           
  geom_tile(color = "white")+
  scale_fill_gradient2(low = 'cornsilk1',  mid = 'gold', high = 'forestgreen', #midpoint = 0.5, 
                       name = 'Kendall correlation\n coefficient tau') + 
  theme_classic() +
  theme(legend.position="right") +
  xlab("Sampling sources") +
  ylab("Antibiotic class") +
  facet_wrap(~Pathogen) +
  geom_text(aes(label = round((tau), digits=2)), 
            color = "white", fontface='bold', size = 6) +
  theme(strip.text.x = element_text(size = 14, face = "italic"),
        axis.text.x = element_text( size = 12, face = "bold"),
        axis.text.y = element_text( size = 12, face = "bold"),
        axis.title = element_text( size = 16, face = "bold"),
        legend.text=element_text(size=10))

ggp

#ggsave("heatmap_tau.jpeg", width = 32, height = 20, units = "cm")

