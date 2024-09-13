
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

df = read.csv(here::here("data", "final_AMR_dataset.csv"))

df2 <- df %>%
  filter(Data != "GLASS") %>%
  filter(!(Source %in% c("HEENT", "Reproduction", "Other"))) %>%
  filter(Antibiotic %in% c("Aminoglycoside", "Aminopenicillin", "Carbapenem", "Cephalosporin_3rd", "Cephalosporin_4th", "Fluoroquinolone", "Sulfonamide_c", "Polypeptide")) %>%
  mutate(Antibiotic = replace(Antibiotic, Antibiotic=="Sulfonamide_c", "Sulfonamide"))

#proportion of resistance
df2$p <- df2$Resistant/df2$Total

####################################################################################
x = df2 %>%
  filter(Source == "Blood/sterile") %>% # source Blood for x axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup %>%
  filter(Total > 10)
  
y = df2 %>%
  #filter(Data != "GLASS") %>%
  filter(Source == "Gastro") %>% # source for y axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup  %>%
  filter(Total > 10)

x_y <- merge(x, y, by = c('Country', 'Region', 'Year', 'Pathogen', 'Antibiotic')) %>%
  filter(!is.nan(p.x)) %>% # remove NaN
  filter(!is.nan(p.y)) %>% # remove NaN
  filter(!is.na(Antibiotic)) %>%
  mutate(Total = Total.y+Total.x,
         Diff = p.y-p.x) %>%
  mutate(Antibiotic = factor(Antibiotic, # list 9 antibiotic classes of interest
                             levels = c("Aminopenicillin", "Aminoglycoside", "Carbapenem",
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
  filter(Source == "Blood/sterile") %>% # source Blood for x axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup %>%
  filter(Total > 10)

y2 = df2 %>%
  #filter(Data != "GLASS") %>%
  filter(Source == "Respiratory") %>% # source for y axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup  %>%
  filter(Total > 10)

x_y2 <- merge(x2, y2, by = c('Country', 'Region', 'Year', 'Pathogen', 'Antibiotic')) %>%
  filter(!is.nan(p.x)) %>% # remove NaN
  filter(!is.nan(p.y)) %>% # remove NaN
  mutate(Total = Total.y+Total.x,
         Diff = p.y-p.x) %>%
  mutate(Antibiotic = factor(Antibiotic, # list 9 antibiotic classes of interest
                             levels = c("Aminopenicillin", "Aminoglycoside", "Carbapenem",
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
  filter(Source == "Blood/sterile") %>% # source Blood for x axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup %>%
  filter(Total > 10)

y3 = df2 %>%
  #filter(Data != "GLASS") %>%
  filter(Source == "Urine") %>% # source for y axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup  %>%
  filter(Total > 10)

x_y3 <- merge(x3, y3, by = c('Country', 'Region', 'Year', 'Pathogen', 'Antibiotic')) %>%
  filter(!is.nan(p.x)) %>% # remove NaN
  filter(!is.nan(p.y)) %>% # remove NaN
  mutate(Total = Total.y+Total.x,
         Diff = p.y-p.x) %>%
  mutate(Antibiotic = factor(Antibiotic, # list 9 antibiotic classes of interest
                             levels = c("Aminopenicillin", "Aminoglycoside", "Carbapenem",
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
  filter(Source == "Blood/sterile") %>% # source Blood for x axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup %>%
  filter(Total > 10)

y4 = df2 %>%
  #filter(Data != "GLASS") %>%
  filter(Source == "Wound") %>% # source for y axis
  group_by(Country, Region, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup  %>%
  filter(Total > 10)

x_y4 <- merge(x4, y4, by = c('Country', 'Region', 'Year', 'Pathogen', 'Antibiotic')) %>%
  filter(!is.nan(p.x)) %>% # remove NaN
  filter(!is.nan(p.y)) %>% # remove NaN
  mutate(Total = Total.y+Total.x,
         Diff = p.y-p.x) %>%
  mutate(Antibiotic = factor(Antibiotic, # list 9 antibiotic classes of interest
                             levels = c("Aminopenicillin", "Aminoglycoside", "Carbapenem",
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

colors = brewer.pal(6, "Set1")

# PLOT CORRELATION:
x_y_label_e = paste0("r = ",
                   round(cor(x_y$p.x[x_y$Pathogen=="Escherichia coli"],
                             x_y$p.y[x_y$Pathogen=="Escherichia coli"],
                             method="kendall"), 2))
x_y_label_k = paste0("r = ",
                   round(cor(x_y$p.x[x_y$Pathogen=="Klebsiella pneumoniae"],
                             x_y$p.y[x_y$Pathogen=="Klebsiella pneumoniae"],
                             method="kendall"), 2))

p1 = ggplot(x_y) +
  facet_wrap(~Pathogen) + 
  geom_point(aes(x = p.x, y = p.y, fill = factor(Region)), 
             pch=21, size = 1) + 
  geom_smooth(aes(x=p.x, y=p.y), method="lm", colour="black") +
  geom_text(data=data.frame(Pathogen="Escherichia coli"), aes(x = 0.2, y=1, label = x_y_label_e)) +
  geom_text(data=data.frame(Pathogen="Klebsiella pneumoniae"), aes(x = 0.2, y=1, label = x_y_label_k)) +
  theme_bw() +
  theme(legend.position = "none", # "none"
        strip.text.x = element_text(size = 11, face = "italic"),
        axis.text.x=element_text(size=11, angle=45, hjust=1),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=11)) + 
  scale_fill_manual(name = "WHO Region:", values = colors) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x = "Res. prop. in Blood/sterile samples", 
       y = "Res. prop. in Gastro samples")


x_y2_label_e = paste0("r = ",
                     round(cor(x_y2$p.x[x_y2$Pathogen=="Escherichia coli"],
                               x_y2$p.y[x_y2$Pathogen=="Escherichia coli"],
                               method="kendall"), 2))
x_y2_label_k = paste0("r = ",
                     round(cor(x_y2$p.x[x_y2$Pathogen=="Klebsiella pneumoniae"],
                               x_y2$p.y[x_y2$Pathogen=="Klebsiella pneumoniae"],
                               method="kendall"), 2))

p2 = ggplot(x_y2) +
  facet_wrap(~Pathogen) + 
  geom_point(aes(x = p.x, y = p.y, fill = factor(Region)), 
             pch=21, size = 1) + 
  geom_smooth(aes(x=p.x, y=p.y), method="lm", colour="black") +
  geom_text(data=data.frame(Pathogen="Escherichia coli"), aes(x = 0.2, y=1, label = x_y2_label_e)) +
  geom_text(data=data.frame(Pathogen="Klebsiella pneumoniae"), aes(x = 0.2, y=1, label = x_y2_label_k)) +
  theme_bw() +
  theme(legend.position = "none", # "none"
        strip.text.x = element_text(size = 11, face = "italic"),
        axis.text.x=element_text(size=11, angle=45, hjust=1),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=12)) + 
  scale_fill_manual(name = "WHO Region:", values = colors) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x = "Res. prop. in Blood/sterile samples", 
       y = "Res. prop. in Respiratory samples")

x_y3_label_e = paste0("r = ",
                      round(cor(x_y3$p.x[x_y3$Pathogen=="Escherichia coli"],
                                x_y3$p.y[x_y3$Pathogen=="Escherichia coli"],
                                method="kendall"), 2))
x_y3_label_k = paste0("r = ",
                      round(cor(x_y3$p.x[x_y3$Pathogen=="Klebsiella pneumoniae"],
                                x_y3$p.y[x_y3$Pathogen=="Klebsiella pneumoniae"],
                                method="kendall"), 2))

p3 = ggplot(x_y3) +
  facet_wrap(~Pathogen) + 
  geom_point(aes(x = p.x, y = p.y, fill = factor(Region)), 
             pch=21, size = 1) + 
  geom_smooth(aes(x=p.x, y=p.y), method="lm", colour="black") +
  geom_text(data=data.frame(Pathogen="Escherichia coli"), aes(x = 0.2, y=1, label = x_y3_label_e)) +
  geom_text(data=data.frame(Pathogen="Klebsiella pneumoniae"), aes(x = 0.2, y=1, label = x_y3_label_k)) +
  theme_bw() +
  theme(legend.position = "none", # "none"
        strip.text.x = element_text(size = 11, face = "italic"),
        axis.text.x=element_text(size=11, angle=45, hjust=1),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=11)) + 
  scale_fill_manual(name = "WHO Region:", values = colors) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x = "Res. prop. in Blood/sterile samples", 
       y = "Res. prop. in Urine samples")

x_y4_label_e = paste0("r = ",
                      round(cor(x_y4$p.x[x_y4$Pathogen=="Escherichia coli"],
                                x_y4$p.y[x_y4$Pathogen=="Escherichia coli"],
                                method="kendall"), 2))
x_y4_label_k = paste0("r = ",
                      round(cor(x_y4$p.x[x_y4$Pathogen=="Klebsiella pneumoniae"],
                                x_y4$p.y[x_y4$Pathogen=="Klebsiella pneumoniae"],
                                method="kendall"), 2))

p4 = ggplot(x_y4) +
  facet_wrap(~Pathogen) + 
  geom_point(aes(x = p.x, y = p.y, fill = factor(Region)), 
             pch=21, size = 1) + 
  geom_smooth(aes(x=p.x, y=p.y), method="lm", colour="black") +
  geom_text(data=data.frame(Pathogen="Escherichia coli"), aes(x = 0.2, y=1, label = x_y4_label_e)) +
  geom_text(data=data.frame(Pathogen="Klebsiella pneumoniae"), aes(x = 0.2, y=1, label = x_y4_label_k)) +
  theme_bw() +
  theme(legend.position = "bottom", # "none"
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        strip.text.x = element_text(size = 11, face = "italic"),
        axis.text.x=element_text(size=11, angle=45, hjust=1),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=11)) + 
  scale_fill_manual(name = "WHO Region:", values = colors) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x = "Res. prop. in Blood/sterile samples", 
       y = "Res. prop. in Wound samples")


df2 = df2 %>%
  group_by(Country, Region, Source, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total),
            Resistant = sum(Resistant)) %>%
  mutate(#Data = "All",
    p = Resistant/Total) %>%
  ungroup %>%
  filter(Total > 10)


corr_results_all = data.frame()

for(bac in unique(df2$Pathogen)){
  
  df_bac = df2 %>% filter(Pathogen==bac)
  
  for(abx in unique(df_bac$Antibiotic)){
    
    df_abx = df_bac %>% filter(Antibiotic==abx)
    
    blood_data = df_abx %>% filter(Source == "Blood/sterile")
    
    for(type in unique(df_abx$Source)){
      
      df_source = df_abx %>% filter(Source==type)
      
      df_reg = inner_join(blood_data, df_source,
                          by=c("Country", "Region", "Year", "Pathogen", "Antibiotic"))
      
      cor_val = cor(df_reg$p.x, df_reg$p.y, method = "kendall")
      corr_results_all = rbind(corr_results_all,
                               data.frame(Pathogen=bac, Antibiotic=abx,
                                          Source=type, Coeff=cor_val))
      
    }
  }
}

p5 = ggplot(corr_results_all %>% filter(Source != "Blood/sterile"), aes(x=Source, y=Antibiotic, fill = Coeff)) +                           
  geom_tile(color = "white")+
  scale_fill_gradient2(low = 'cornsilk1',  mid = 'gold', high = 'forestgreen', #midpoint = 0.5, 
                       name = 'Correlation coefficient') + 
  theme_bw() +
  theme(legend.position="right") +
  xlab("Sampling sources") +
  ylab("Antibiotic class") +
  facet_wrap(~Pathogen, ncol=1) +
  geom_text(aes(label = round(Coeff, digits=2)), 
            color = "white", fontface='bold', size = 3) +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 11, angle=45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12))

plot_grid(plot_grid(plot_grid(p1, p2, p3, p4+theme(legend.position = "none"), 
                    labels = c("a)", "b)", "c)", "d)"),
                    hjust=0, ncol = 2, nrow = 2),
          get_plot_component(p4+guides(fill=guide_legend(override.aes = list(size=3))),
                             "guide-box-bottom", return_all = T),
          nrow=2, rel_heights = c(1,0.1)),
          p5+theme(legend.position="none"), labels = c("", "e)"),
          rel_widths = c(1,0.55))


ggsave(here("figures", "fig2.png"), height=8, width = 12, bg = "white")
