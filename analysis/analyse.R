
library(dplyr)
library(tidyr)
library(readxl)
library(tidyverse)
library(here)
library(AMR)
library(reshape2)

df_AMR = read.csv(here::here("data", "final_AMR_dataset.csv"))

df_AMR2 = df_AMR %>%
  filter(Total != 0) %>%
  group_by(Country, Source, Year, Pathogen, Antibiotic) %>%
  summarise(Total = sum(Total), Resistant = sum(Resistant)) %>%
  ungroup %>%
  filter(Resistant >= 10) %>%
  select(Source, Pathogen, Antibiotic) %>%
  distinct() %>%
  group_by(Pathogen, Antibiotic) %>%
  summarise(Source = paste0(unique(Source), collapse = "-")) %>%
  ungroup %>%
  filter(grepl("Blood", Source)) %>%
  filter(grepl("-", Source)) 

unique(df_AMR2$Pathogen)
unique(df_AMR2$Antibiotic)

# 5 bacteria and 37 antibiotics for which there seems to be some data (97 combinations)
