library(dplyr)
library(tidyr)
library(readxl)
library(tidyverse)
library(here)
library(AMR)
library(reshape2)

df_AMR = read.csv(here::here("data", "final_full_AMR_dataset.csv"))

df_AMR2 = df_AMR %>%
  filter(Data != "GLASS") %>%
  filter(Total != 0) %>%
  filter(!(Antibiotic %in% c("Aminoglycoside", "Aminopenicillin", "Antipseudo_penicillin", "Carbapenem", "Cephalosporin_3rd", "Cephalosporin_4th", "Fluoroquinolone", "Sulpha"))) %>%
  group_by(Country, Source, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total), Resistant = sum(Resistant)) %>%
  ungroup

df_AMR3 = df_AMR2 %>%
  group_by(Source, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total), Resistant = sum(Resistant)) %>%
  ungroup %>%
  filter(Total >= 1000) %>%
  group_by(Pathogen, Antibiotic) %>%
  summarise(Source = paste0(unique(Source), collapse = "-"),
            Total = sum(Total), Resistant = sum(Resistant)) %>%
  ungroup %>%
  filter(grepl("Blood", Source)) %>%
  filter(grepl("-", Source))

write.csv(df_AMR3, here("data", "other_combinations_possible.csv"), row.names = F)

unique(df_AMR3$Pathogen)
unique(df_AMR3$Antibiotic)
