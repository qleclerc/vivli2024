
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)

df_AMR = read.csv(here::here("data", "final_AMR_dataset.csv")) %>%
  filter(Data != "GLASS") %>%
  filter(!(Source %in% c("HEENT", "Reproduction", "Other"))) %>%
  filter(!is.na(Region)) %>%
  group_by(Region, Year, Source, Pathogen, Antibiotic) %>%
  summarise(Resistant = sum(Resistant),
            Total = sum(Total)) %>%
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


ggplot(results_all %>% filter(Value=="Resistant.x" & pval < 0.05)) +
  # geom_point(aes(Antibiotic, exp(Estimate), colour=Region),
  #                 position = position_dodge(0.5), size=1) +
  geom_pointrange(aes(Antibiotic, exp(Estimate), ymin=exp(Estimate-std), ymax=exp(Estimate+std), colour=Region),
             position = position_dodge(0.5), size=0.15) +
  geom_hline(yintercept = 1, linetype="dashed", alpha=0.2) +
  facet_grid(cols=vars(Pathogen), rows=vars(Source)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "bottom")

ggsave("reg_by_reg_exp.png", height=6, width=7)

ggplot(results_all %>% filter(Value=="(Intercept)" & pval < 0.05)) +
  # geom_point(aes(Antibiotic, exp(Estimate), colour=Region),
  #                 position = position_dodge(0.5), size=1) +
  geom_pointrange(aes(Antibiotic, exp(Estimate), ymin=exp(Estimate-std), ymax=exp(Estimate+std), colour=Region),
                  position = position_dodge(0.5), size=0.15) +
  geom_hline(yintercept = 1, linetype="dashed", alpha=0.2) +
  facet_grid(cols=vars(Pathogen), rows=vars(Source)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "bottom")

ggsave("reg_by_reg_exp_beta0.png", height=6, width=7)


