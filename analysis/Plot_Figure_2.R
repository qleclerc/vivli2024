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
library(ggpubr)
library(cowplot)


df <- read.csv("final_AMR_dataset_edit_final.csv")
str(df)
unique(df$Antibiotic)

# 9 antibiotic classes of interest
# Aminopenicillin, Antipseudo_penicillin, Aminoglycoside, Carbapenem, 
# Cephalosporin_3rd, Cephalosporin_4th, 
# Fluoroquinolone, Sulfonamide

# rename Sulpha to Sulfonamide and Colistin to Polypeptide and select antibiotic classes of interest:
ATB_classes <- c("Aminopenicillin", "Antipseudo_penicillin", 
                 "Aminoglycoside", "Carbapenem",
                 "Cephalosporin_3rd", "Cephalosporin_4th",
                 "Fluoroquinolone", "Polypeptide",
                 "Sulfonamide")

df2 <- df %>%
  mutate(Antibiotic = replace(Antibiotic, Antibiotic =="Sulpha", "Sulfonamide")) %>%
  mutate(Antibiotic = replace(Antibiotic, Antibiotic =="Colistin", "Polypeptide")) %>%
  filter(Antibiotic %in% ATB_classes)   # filter out ATB classes from antibiotics

#proportion of resistance
df2$p <- df2$Resistant/df2$Total

####################################################################################
x = df2 %>%
  filter(Source == "Blood") %>% # source Blood for x axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup

y = df2 %>%
  #filter(Data != "GLASS") %>%
  filter(Source == "Gastro") %>% # source for y axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup 

x_y <- merge(x, y, by = c('Country', 'Region', 'Year', 'Pathogen', 'Antibiotic')) %>%
  filter(!is.nan(p.x)) %>% # remove NaN
  filter(!is.nan(p.y)) %>% # remove NaN
  filter(!is.na(Antibiotic)) %>%
  mutate(Total = Total.y+Total.x,
         Diff = p.y-p.x) %>%
  mutate(Antibiotic = factor(Antibiotic, # list 9 antibiotic classes of interest
                             levels = c("Aminopenicillin", "Antipseudo_penicillin", 
                                        "Aminoglycoside", "Carbapenem",
                                        "Cephalosporin_3rd", "Cephalosporin_4th",
                                        "Fluoroquinolone", "Polypeptide",
                                        "Sulfonamide")))

x_y_text = data.frame()

for(pathogen in unique(x_y$Pathogen)){
  x_y_p = x_y %>%
    filter(Pathogen == pathogen)
  
  x_y_text = rbind(x_y_text,
                   data.frame(Pathogen = pathogen,
                              Comparisons = nrow(x_y_p),
                              Data_points = sum(x_y_p$Total.y)
                   ))
}

####################################################################################
####################################################################################
x2 = df2 %>%
  filter(Source == "Blood") %>% # source Blood for x axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup

y2 = df2 %>%
  #filter(Data != "GLASS") %>%
  filter(Source == "Respiratory") %>% # source for y axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup 

x_y2 <- merge(x2, y2, by = c('Country', 'Region', 'Year', 'Pathogen', 'Antibiotic')) %>%
  filter(!is.nan(p.x)) %>% # remove NaN
  filter(!is.nan(p.y)) %>% # remove NaN
  mutate(Total = Total.y+Total.x,
         Diff = p.y-p.x) %>%
  mutate(Antibiotic = factor(Antibiotic, # list 9 antibiotic classes of interest
                             levels = c("Aminopenicillin", "Antipseudo_penicillin", 
                                        "Aminoglycoside", "Carbapenem",
                                        "Cephalosporin_3rd", "Cephalosporin_4th",
                                        "Fluoroquinolone", "Polypeptide",
                                        "Sulfonamide")))

x_y_text2 = data.frame()

for(pathogen in unique(x_y2$Pathogen)){
  x_y_p2 = x_y2 %>%
    filter(Pathogen == pathogen)
  
  x_y_text2 = rbind(x_y_text2,
                   data.frame(Pathogen = pathogen,
                              Comparisons = nrow(x_y_p2),
                              Data_points = sum(x_y_p2$Total.y)
                   ))
}

####################################################################################
x3 = df2 %>%
  filter(Source == "Blood") %>% # source Blood for x axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup

y3 = df2 %>%
  #filter(Data != "GLASS") %>%
  filter(Source == "Urine") %>% # source for y axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup 

x_y3 <- merge(x3, y3, by = c('Country', 'Region', 'Year', 'Pathogen', 'Antibiotic')) %>%
  filter(!is.nan(p.x)) %>% # remove NaN
  filter(!is.nan(p.y)) %>% # remove NaN
  mutate(Total = Total.y+Total.x,
         Diff = p.y-p.x) %>%
  mutate(Antibiotic = factor(Antibiotic, # list 9 antibiotic classes of interest
                             levels = c("Aminopenicillin", "Antipseudo_penicillin", 
                                        "Aminoglycoside", "Carbapenem",
                                        "Cephalosporin_3rd", "Cephalosporin_4th",
                                        "Fluoroquinolone", "Polypeptide",
                                        "Sulfonamide")))

x_y_text3 = data.frame()

for(pathogen in unique(x_y3$Pathogen)){
  x_y_p3 = x_y3 %>%
    filter(Pathogen == pathogen)
  
  x_y_text3 = rbind(x_y_text3,
                    data.frame(Pathogen = pathogen,
                               Comparisons = nrow(x_y_p3),
                               Data_points = sum(x_y_p3$Total.y)
                    ))
}

####################################################################################
x4 = df2 %>%
  filter(Source == "Blood") %>% # source Blood for x axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup

y4 = df2 %>%
  #filter(Data != "GLASS") %>%
  filter(Source == "Wound") %>% # source for y axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup 

x_y4 <- merge(x4, y4, by = c('Country', 'Region', 'Year', 'Pathogen', 'Antibiotic')) %>%
  filter(!is.nan(p.x)) %>% # remove NaN
  filter(!is.nan(p.y)) %>% # remove NaN
  mutate(Total = Total.y+Total.x,
         Diff = p.y-p.x) %>%
  mutate(Antibiotic = factor(Antibiotic, # list 9 antibiotic classes of interest
                             levels = c("Aminopenicillin", "Antipseudo_penicillin", 
                                        "Aminoglycoside", "Carbapenem",
                                        "Cephalosporin_3rd", "Cephalosporin_4th",
                                        "Fluoroquinolone", "Polypeptide",
                                        "Sulfonamide")))

x_y_text4 = data.frame()

for(pathogen in unique(x_y4$Pathogen)){
  x_y_p4 = x_y4 %>%
    filter(Pathogen == pathogen)
  
  x_y_text4 = rbind(x_y_text4,
                    data.frame(Pathogen = pathogen,
                               Comparisons = nrow(x_y_p4),
                               Data_points = sum(x_y_p4$Total.y)
                    ))
}

####################################################################
#################################### Correlation Plot, Figure 2


colors = c("African Region" = "red",
           "European Region" = "chartreuse3",
           "Region of the Americas" = "cyan",
           "South-East Asian Region" = "purple",
           "Western Pacific Region" = "gold",
           "Eastern Mediterranean Region" = "orange")


# PLOT CORRELATION:
p1 <- ggscatter(x_y, x = "p.x", y = "p.y", 
          add = "reg.line", conf.int = TRUE, cor.coef.size = 6,
          cor.coef = TRUE, cor.method = "kendall", # kendall
          xlab = "Resistance Proportion in Blood/Sterile Samples", 
          ylab = "Resistance Proportion in Gastro Samples") + 
  facet_wrap(~Pathogen) + 
  theme(legend.position = "none", # "none"
        strip.text.x = element_text(size = 12, face = "italic"),
        legend.text=element_text(size=12),
        #legend.title=element_blank(),
        legend.background=element_blank()) + 
  geom_point(aes(x = p.x, y = p.y, fill = factor(Region)), pch=21, colour="black", size = 2) + 
  scale_fill_manual(name = "WHO Region:", values = colors) 

p2 <- ggscatter(x_y2, x = "p.x", y = "p.y", 
                add = "reg.line", conf.int = TRUE, cor.coef.size = 6,
                cor.coef = TRUE, cor.method = "kendall", # kendall
                xlab = "Resistance Proportion in Blood/Sterile Samples", 
                ylab = "Resistance Proportion in Respiratory Samples") + 
  facet_wrap(~Pathogen) + 
  theme(legend.position = "none", # "none"
        strip.text.x = element_text(size = 12, face = "italic"),
        legend.text=element_text(size=12),
        #legend.title=element_blank(),
        legend.background=element_blank()) + 
  geom_point(aes(x = p.x, y = p.y, fill = factor(Region)), pch=21, colour="black", size = 2) + 
  scale_fill_manual(name = "WHO Region:", values = colors) 


p3 <- ggscatter(x_y3, x = "p.x", y = "p.y", 
                add = "reg.line", conf.int = TRUE, cor.coef.size = 6,
                cor.coef = TRUE, cor.method = "kendall", # kendall
                xlab = "Resistance Proportion in Blood/Sterile Samples", 
                ylab = "Resistance Proportion in Urine Samples") + 
  facet_wrap(~Pathogen) + 
  theme(legend.position = "none", # "none"
        strip.text.x = element_text(size = 12, face = "italic"),
        legend.text=element_text(size=12),
        #legend.title=element_blank(),
        legend.background=element_blank()) + 
  geom_point(aes(x = p.x, y = p.y, fill = factor(Region)), pch=21, colour="black", size = 2) + 
  scale_fill_manual(name = "WHO Region:", values = colors)


p4 <- ggscatter(x_y4, x = "p.x", y = "p.y", 
                add = "reg.line", conf.int = TRUE, cor.coef.size = 6,
                cor.coef = TRUE, cor.method = "kendall", # kendall
                xlab = "Resistance Proportion in Blood/Sterile Samples", 
                ylab = "Resistance Proportion in Wound Samples") + 
  facet_wrap(~Pathogen) + 
  theme(legend.position = "none", # "none"
        strip.text.x = element_text(size = 12, face = "italic"),
        legend.text=element_text(size=12),
        #legend.title=element_blank(),
        legend.background=element_blank()) + 
  geom_point(aes(x = p.x, y = p.y, fill = factor(Region)), pch=21, colour="black", size = 2) + 
  scale_fill_manual(name = "WHO Region:", values = colors) 


Fig2 <- plot_grid(p1, p2, p3, p4, 
                 labels = c("A", "B", "C", "D"), label_size = 24,
                 ncol = 2, nrow = 2)

#ggsave("Fig2.jpeg", width = 35, height = 28, units = "cm")
