
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(here)

df_AMR = read.csv(here::here("data", "final_AMR_dataset.csv")) %>%
  filter(Data != "GLASS") %>%
  filter(!(Source %in% c("HEENT", "Reproduction", "Other"))) %>%
  filter(!is.na(Region)) %>%
  group_by(Region, Year, Source, Pathogen, Antibiotic) %>%
  summarise(Resistant = sum(Resistant),
            Total = sum(Total)) %>%
  filter(Total >= 10) %>%
  mutate(prop = Resistant/Total) %>%
  filter(!(is.nan(prop))) %>%
  filter(prop > 0) %>%
  filter(Antibiotic %in% c("Aminoglycoside", "Aminopenicillin", "Carbapenem", "Cephalosporin_3rd", "Cephalosporin_4th", "Fluoroquinolone", "Sulfonamide_c", "Polypeptide")) %>%
  mutate(Antibiotic = replace(Antibiotic, Antibiotic=="Sulfonamide_c", "Sulfonamide"))

df_sterile = df_AMR %>%
  filter(Source == "Blood/sterile")

df_for_reg = df_AMR %>%
  filter(Source!="Blood/sterile") %>%
  inner_join(df_sterile, by=c("Region", "Year", "Pathogen", "Antibiotic"))

ggplot(df_for_reg) +
  geom_point(aes(Resistant.x, Resistant.y)) +
  facet_grid(cols=vars(Pathogen), rows=vars(Source.x), scales = "free") +
  theme_bw()

bug = df_for_reg$Pathogen[1]
drug = df_for_reg$Antibiotic[1]
type = df_for_reg$Source.x[1]
region = df_for_reg$Region[1]

results_all = data.frame()

for(bug in unique(df_for_reg$Pathogen)){
  
  df_bug = df_for_reg %>% filter(Pathogen == bug)
  
  for(drug in unique(df_bug$Antibiotic)){
    
    df_drug = df_bug %>% filter(Antibiotic == drug)
    
    for(region in unique(df_drug$Region)){
      
      df_region = df_drug %>% filter(Region == region)
      
      for(type in unique(df_region$Source.x)){
        
        df_type = df_region %>% filter(Source.x == type)
        
        tt = glm(data = df_type, Resistant.y~Resistant.x, family = "poisson")
        tt = coef(summary(tt))
        
        results_all = rbind(results_all,
                            data.frame(Pathogen = bug, Antibiotic = drug,
                                       Region = region, Source = type,
                                       Value=row.names(tt),
                                       Estimate = tt[,"Estimate"],
                                       std = tt[,"Std. Error"],
                                       pval = tt[,"Pr(>|z|)"],
                                       row.names = NULL))
      }
    }
  }
}

results_all$Estimate = exp(results_all$Estimate)

results_all$cat_Estimate = cut(results_all$Estimate, breaks = c(0,50,100,200,9999),
                               labels = c("0-50","50-100","100-200",">200"))

pa=ggplot(results_all %>% filter(Value=="(Intercept)" & pval < 0.05), aes(x=Source, y=Antibiotic, fill = cat_Estimate)) +                           
  geom_tile(color = "white")+
  scale_fill_discrete(type = c('grey50', '#87CEEB', "#4682B4", '#003366'), 
                      name = 'Intercept') + 
  theme_bw() +
  xlab("Sampling sources") +
  ylab("Antibiotic class") +
  facet_grid(cols=vars(Region), rows=vars(Pathogen)) +
  geom_text(aes(label = round(Estimate, digits=0)), 
            color = "white", fontface='bold', size = 3) +
  theme(legend.position="bottom",
        strip.text.y = element_text(size = 11, face = "italic"),
        strip.text.x = element_text(size=11),
        axis.text.x = element_text(size = 11, angle=45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12))

results_all$cat_Estimate = cut(results_all$Estimate, breaks = c(0.5,0.9,1.1,1.2,1.4,9999),
                               labels = c("0.5-0.9","0.9-1.1","1.1-1.2","1.2-1.4",">1.4"))

pb=ggplot(results_all %>% filter(Value=="Resistant.x" & pval < 0.05), aes(x=Source, y=Antibiotic, fill = cat_Estimate)) +                           
  geom_tile(color = "white")+
  scale_fill_discrete(type = c("orange2", 'grey50', '#87CEEB', "#4682B4", '#003366'), 
                      name = 'Regression coefficient') + 
  theme_bw() +
  xlab("Sampling sources") +
  ylab("Antibiotic class") +
  facet_grid(cols=vars(Region), rows=vars(Pathogen)) +
  geom_text(aes(label = round(Estimate, digits=2)), 
            color = "white", fontface='bold', size = 3) +
  theme(legend.position="bottom",
        strip.text.y = element_text(size = 11, face = "italic"),
        strip.text.x = element_text(size=11),
        axis.text.x = element_text(size = 11, angle=45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12))

plot_grid(pa,pb, labels=c("a)", "b)"), hjust=0, nrow=2)

ggsave(here("figures","fig3.png"), height=12, width=15)

