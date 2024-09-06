
################################################################################

################################################################################
######################### AMR Vivli data challenge #############################
################################################################################

################################################################################

#libraries required
library(dplyr)
library(tidyr)
library(readxl)
library(tidyverse)
library(here)
library(AMR)
library(reshape2)

# helper function to define class resistance
class_resistance = function(dataset){
  
  # Fluoroquinolones
  
  cnames = c("CIP", "LVX", "OFX")
  cnames = cnames[!cnames%in%names(dataset)]
  if(length(cnames)!=0) dataset[cnames] = NA
  
  dataset = dataset %>%
    mutate(Fluoroquinolone = CIP) %>%
    mutate(Fluoroquinolone = ifelse(is.na(Fluoroquinolone), LVX, Fluoroquinolone)) %>%
    mutate(Fluoroquinolone = ifelse(is.na(Fluoroquinolone), OFX, Fluoroquinolone))
  
  if(all(is.na(dataset$Fluoroquinolone))) warning("Fluoroquinolone not added - no susceptibility data available")
  
  # Aminoglycosides
  
  cnames = c("GEN", "TOB", "AMK", "KAN", "NET", "STR1")
  cnames = cnames[!cnames%in%names(dataset)]
  if(length(cnames)!=0) dataset[cnames] = NA
  
  dataset = dataset %>%
    mutate(Aminoglycoside = GEN) %>%
    mutate(Aminoglycoside = ifelse(is.na(Aminoglycoside), TOB, Aminoglycoside)) %>%
    mutate(Aminoglycoside = ifelse(is.na(Aminoglycoside), AMK, Aminoglycoside)) %>%
    mutate(Aminoglycoside = ifelse(is.na(Aminoglycoside), KAN, Aminoglycoside)) %>%
    mutate(Aminoglycoside = ifelse(is.na(Aminoglycoside), NET, Aminoglycoside)) %>%
    mutate(Aminoglycoside = ifelse(is.na(Aminoglycoside), STR1, Aminoglycoside))
  
  if(all(is.na(dataset$Aminoglycoside))) warning("Aminoglycoside not added - no susceptibility data available")
  
  # Cephalosporins_3rd
  
  cnames = c("CTX", "CRO", "CZO", "CDR", "CFM")
  cnames = cnames[!cnames%in%names(dataset)]
  if(length(cnames)!=0) dataset[cnames] = NA
  
  dataset = dataset %>%
    mutate(Cephalosporin_3rd = CTX) %>%
    mutate(Cephalosporin_3rd = ifelse(is.na(Cephalosporin_3rd), CRO, Cephalosporin_3rd)) %>%
    mutate(Cephalosporin_3rd = ifelse(is.na(Cephalosporin_3rd), CZO, Cephalosporin_3rd)) %>%
    mutate(Cephalosporin_3rd = ifelse(is.na(Cephalosporin_3rd), CDR, Cephalosporin_3rd)) %>%
    mutate(Cephalosporin_3rd = ifelse(is.na(Cephalosporin_3rd), CFM, Cephalosporin_3rd))
  
  if(all(is.na(dataset$Cephalosporin_3rd))) warning("Cephalosporin_3rd not added - no susceptibility data available")
  
  # Carbapenems
  
  cnames = c("IPM", "MEM", "ETP")
  cnames = cnames[!cnames%in%names(dataset)]
  if(length(cnames)!=0) dataset[cnames] = NA
  
  dataset = dataset %>%
    mutate(Carbapenem = IPM) %>%
    mutate(Carbapenem = ifelse(is.na(Carbapenem), MEM, Carbapenem)) %>%
    mutate(Carbapenem = ifelse(is.na(Carbapenem), ETP, Carbapenem))
  
  if(all(is.na(dataset$Carbapenem))) warning("Carbapenem not added - no susceptibility data available")
  
  # Cephalosporins_4th
  
  cnames = c("FEP")
  cnames = cnames[!cnames%in%names(dataset)]
  if(length(cnames)!=0) dataset[cnames] = NA
  
  dataset = dataset %>%
    mutate(Cephalosporin_4th = FEP) 
  
  if(all(is.na(dataset$Cephalosporin_4th))) warning("Cephalosporin_4th not added - no susceptibility data available")
  
  # Antipseudo_penicillins
  
  cnames = c("AMC", "SAM")
  cnames = cnames[!cnames%in%names(dataset)]
  if(length(cnames)!=0) dataset[cnames] = NA
  
  dataset = dataset %>%
    mutate(Antipseudo_penicillin = AMC) %>%
    mutate(Antipseudo_penicillin = ifelse(is.na(Antipseudo_penicillin), SAM, Antipseudo_penicillin))
  
  if(all(is.na(dataset$Antipseudo_penicillin))) warning("Antipseudo_penicillin not added - no susceptibility data available")
  
  # Sulpha
  
  cnames = c("SXT", "SSS")
  cnames = cnames[!cnames%in%names(dataset)]
  if(length(cnames)!=0) dataset[cnames] = NA
  
  dataset = dataset %>%
    mutate(Sulpha = SXT) %>%
    mutate(Sulpha = ifelse(is.na(Sulpha), SSS, Sulpha))
  
  if(all(is.na(dataset$Sulpha))) warning("Sulpha not added - no susceptibility data available")
  
  # Aminopenicillins
  
  cnames = c("AMP", "AMX", "PIP", "MEC")
  cnames = cnames[!cnames%in%names(dataset)]
  if(length(cnames)!=0) dataset[cnames] = NA
  
  dataset = dataset %>%
    mutate(Aminopenicillin = AMP) %>%
    mutate(Aminopenicillin = ifelse(is.na(Aminopenicillin), AMX, Aminopenicillin)) %>%
    mutate(Aminopenicillin = ifelse(is.na(Aminopenicillin), PIP, Aminopenicillin)) %>%
    mutate(Aminopenicillin = ifelse(is.na(Aminopenicillin), MEC, Aminopenicillin))
  
  if(all(is.na(dataset$Aminopenicillin))) warning("Aminopenicillin not added - no susceptibility data available")
  
  # remove added columns with only NA
  dataset = dataset[, !apply(dataset, 2, function(x) all(is.na(x)))]
  
  return(dataset)
}


################################################################################
########################### Load AMR datasets ##################################
################################################################################

if(!dir.exists(here("data", "raw"))) stop("Please place the industry surveillance datasets in a new subfolder of the \"data\" folder, named \"raw\"")

#industry
df_ATLAS <- read.csv(here("data", "raw", "ATLAS.csv"))
df_SIDERO <- read_excel(here("data", "raw", "SIDERO.xlsx"))
df_GEARS <- read_excel(here("data", "raw", "GEARS.xlsx"))
df_KEYSTONE <- read_excel(here("data", "raw", "KEYSTONE.xlsx"))

#public
#I'm doing some GLASS reformatting here already, easier than doing it later
df_GLASS = read.csv(here::here("data", "glass_combined.csv"))
df_GLASS = df_GLASS %>%
  mutate(PathogenName = replace(PathogenName, grepl("Salmonella", PathogenName), "Salmonella spp"),
         PathogenName = replace(PathogenName, grepl("Acinetobacter", PathogenName), "Acinetobacter baumannii"),
         PathogenName = replace(PathogenName, grepl("Shigella", PathogenName), "Shigella sonnei"))
df_GLASS = df_GLASS %>%
  mutate(AbTargets = replace(AbTargets, AbTargets %in% c("Meticillin", "Cefoxitin"), "Oxacillin")) %>%
  group_by(Iso3, CountryTerritoryArea, WHORegionName, Year, Specimen, PathogenName, AbTargets) %>%
  mutate(TotalSpecimenIsolates = sum(TotalSpecimenIsolates),
         InterpretableAST = sum(InterpretableAST),
         Resistant = sum(Resistant)) %>%
  mutate(PercentResistant = Resistant/InterpretableAST*100) %>%
  ungroup


source_index = read_xlsx(here::here("data", "source_index.xlsx"), skip=3)

################################################################################
############### Choose variable of interests across datasets ###################
################################################################################

#Countries

#Years
# years_of_interest = c(2018:2019)
years_of_interest = c() #c(2017:2020)

#Bacterial species 
# bacteria_of_interest = c()
bacteria_of_interest = as.mo(c("E. coli", "K pneumoniae"))

#Antibiotics tested
#ESBL (3GC) and carbapenems (CBP)
antibiotics_of_interest = as.ab(c("Ampicillin", "Cefepime", "Cefotaxime", "Ceftazidime",
                                  "Ceftriaxone", "Ciprofloxacin", "Co-trimoxazole",
                                  "Colistin", "Ertapenem", "Imipenem", "Levofloxacin",
                                  "Meropenem", "Nitrofurantoin", "Amoxicillin/clavulanic acid",
                                  "Piperacillin/tazobactam", "Fosfomycin", "Gentamicin", "Amikacin"))
# antibiotics_of_interest = c()


################################################################################
############### Get same variables for all AMR datasets ########################
################################################################################

#### ATLAS ####
df_ATLAS_2 = df_ATLAS[,-c(112:135, 108)] #col 108 is ceftibuten avibactam, removing otherwise picked up as dup of ceftibuten, pre clinical anyways
df_ATLAS_2 <- df_ATLAS_2[,c(grep(pattern = "_I", colnames(df_ATLAS_2), invert=T))]
# MRSA harmonisation - MRSA status indicated by "MRSA" label in Phenotype, not necessarily oxa MIC
# so, set all MRSA phenotypes to have oxa MIC of 8
df_ATLAS_2$Oxacillin[df_ATLAS_2$Phenotype == "MSSA"] = 0.5
df_ATLAS_2$Oxacillin[df_ATLAS_2$Phenotype == "MRSA"] = 8
# ESBL harmonisation - ESBL status indicated by "ESBL" label in Phenotype, not necessarily MIC
# so, set all ESBL phenotypes to have ampicillin MIC of 32
df_ATLAS_2$Ampicillin[df_ATLAS_2$Phenotype == "ESBL"] = 32
# Acinetobacter harmonisation - assume all spp are baumannii
df_ATLAS_2$Species[grepl("Acinetobacter spp", df_ATLAS_2$Species)] = "Acinetobacter baumannii"
# Salmonella harmonisation - assume all are spp
df_ATLAS_2$Species[grepl("Salmonella", df_ATLAS_2$Species)] = "Salmonella spp"

df_ATLAS_2 <- df_ATLAS_2[,c(5,7,8,10,12,3,14:ncol(df_ATLAS_2))]
colnames(df_ATLAS_2)[-c(1:6)] = as.ab(colnames(df_ATLAS_2)[-c(1:6)])
df_ATLAS_2$Species = as.mo(df_ATLAS_2$Species)

#get antibiotics
if(!(all(is.na(antibiotics_of_interest)))) df_ATLAS_2 <- df_ATLAS_2[,c(1:6, which(colnames(df_ATLAS_2) %in% antibiotics_of_interest))]
#get bacteria
if(!(all(is.na(bacteria_of_interest)))) df_ATLAS_2 <- df_ATLAS_2 %>% filter(Species %in% bacteria_of_interest)
#get years
if(length(years_of_interest) != 0) df_ATLAS_2 <- df_ATLAS_2 %>% filter(Year %in% years_of_interest)

if(nrow(df_ATLAS_2) > 0){
  #format resistances
  df_ATLAS_2 = df_ATLAS_2 %>%
    mutate_at(vars(-Country, -Gender, -Age.Group, -Source, -Year, -Species), as.mic)
  df_ATLAS_2 = df_ATLAS_2[,which(apply(df_ATLAS_2, 2, function(x) !(all(is.na(x)))))]
  
  df_ATLAS_2b = df_ATLAS_2 %>%
    mutate_if(is.mic, as.sir, mo = .$Species, guideline = "CLSI")
  
  # NEW: I'm now applying CLSI breakpoints to all datasets in a first instance
  # however, sometimes there are no CLSI breakpoints, so for those I try again with
  # EUCAST breakpoints. This process is the same for all datasets.
  if(any(apply(df_ATLAS_2b, 2, function(x) all(is.na(x))))){
    cat("\nRetrying some columns with EUCAST guidelines...\n")
    to_retry = which(apply(df_ATLAS_2b, 2, function(x) all(is.na(x))))
    df_ATLAS_2b[,to_retry] = df_ATLAS_2 %>%
      select(all_of(c(6,to_retry))) %>%
      mutate_if(is.mic, as.sir, mo = .$Species, guideline = "EUCAST") %>%
      select(-1)
  }
  df_ATLAS_2 = df_ATLAS_2b
  rm(df_ATLAS_2b)
}

df_ATLAS_2$Source = str_to_lower(df_ATLAS_2$Source)
for(i in 1:nrow(source_index)){
  df_ATLAS_2$Source[str_which(df_ATLAS_2$Source, source_index$contains[i])] = source_index$match[i]
}

df_ATLAS_2$Age.Group[df_ATLAS_2$Age.Group %in% c("0 to 2 Years", "13 to 18 Years", "3 to 12 Years")] = "0 to 18 Years"
df_ATLAS_2$Age.Group[df_ATLAS_2$Age.Group %in% c("65 to 84 Years", "85 and Over")] = "65 and Over"

df_ATLAS_2$Gender = str_sub(df_ATLAS_2$Gender, 1, 1)
df_ATLAS_2$Gender[df_ATLAS_2$Gender == ""] = "Unknown"

df_ATLAS_2 = class_resistance(df_ATLAS_2)

#####

#### GEARS ####
df_GEARS_2 <- df_GEARS[,c(6,8,7,9,2,3,11:ncol(df_GEARS))]
colnames(df_GEARS_2)[colnames(df_GEARS_2)=="C_MIC"] = "CZT_MIC"
colnames(df_GEARS_2)[-c(1:6)] = as.ab(colnames(df_GEARS_2)[-c(1:6)])
df_GEARS_2$Organism = as.mo(df_GEARS_2$Organism)
# no MRSA here so no harmonisation needed

#get antibiotics
if(!(all(is.na(antibiotics_of_interest)))) df_GEARS_2 <- df_GEARS_2[,c(1:6, which(colnames(df_GEARS_2) %in% antibiotics_of_interest))]
#get bacteria
if(!(all(is.na(bacteria_of_interest)))) df_GEARS_2 <- df_GEARS_2 %>% filter(Organism %in% bacteria_of_interest)
#get years
if(length(years_of_interest) != 0) df_GEARS_2 <- df_GEARS_2 %>% filter(Year %in% years_of_interest)


if(nrow(df_GEARS_2) > 0){
  #format resistances
  df_GEARS_2 = df_GEARS_2 %>%
    mutate_at(vars(-Country, -Gender, -Age, -BodySite, -Year, -Organism), as.numeric) %>%
    mutate_at(vars(-Country, -Gender, -Age, -BodySite, -Year, -Organism), round, digits = 3) %>%
    mutate_at(vars(-Country, -Gender, -Age, -BodySite, -Year, -Organism), as.mic)
  df_GEARS_2 = df_GEARS_2[,which(apply(df_GEARS_2, 2, function(x) !(all(is.na(x)))))]
  
  df_GEARS_2b = df_GEARS_2 %>%
    mutate_if(is.mic, as.sir, mo = .$Organism, guideline = "CLSI")
  
  if(any(apply(df_GEARS_2b, 2, function(x) all(is.na(x))))){
    cat("\nRetrying some columns with EUCAST guidelines...\n")
    to_retry = which(apply(df_GEARS_2b, 2, function(x) all(is.na(x))))
    df_GEARS_2b[,to_retry] = df_GEARS_2 %>%
      select(all_of(c(6,to_retry))) %>%
      mutate_if(is.mic, as.sir, mo = .$Organism, guideline = "EUCAST") %>%
      select(-1)
  }
  df_GEARS_2 = df_GEARS_2b
  rm(df_GEARS_2b)
  
}

df_GEARS_2$BodySite = str_to_lower(df_GEARS_2$BodySite)
for(i in 1:nrow(source_index)){
  df_GEARS_2$BodySite[str_which(df_GEARS_2$BodySite, source_index$contains[i])] = source_index$match[i]
}

df_GEARS_2$Age[df_GEARS_2$Age %in% c(0:18)] = 0
df_GEARS_2$Age[df_GEARS_2$Age %in% c(19:64)] = 19
df_GEARS_2$Age[df_GEARS_2$Age %in% c(65:150)] = 65

df_GEARS_2$Age[df_GEARS_2$Age %in% c(-1)] = "Unknown"
df_GEARS_2$Age[df_GEARS_2$Age %in% c(0)] = "0 to 18 Years"
df_GEARS_2$Age[df_GEARS_2$Age %in% c(19)] = "19 to 64 Years"
df_GEARS_2$Age[df_GEARS_2$Age %in% c(65)] = "65 and Over"
df_GEARS_2$Age[df_GEARS_2$Age %in% "NG"] = "Unknown"

df_GEARS_2$Gender[df_GEARS_2$Gender == "N"] = "Unknown"

df_GEARS_2 = class_resistance(df_GEARS_2)

#####

#### KEYSTONE ####
df_KEYSTONE_2 <- df_KEYSTONE[,c(34, 38, 37, 41, 2, 3, c(4:32))]
colnames(df_KEYSTONE_2)[-c(1:6)] = as.ab(colnames(df_KEYSTONE_2)[-c(1:6)])
colnames(df_KEYSTONE_2)[colnames(df_KEYSTONE_2) == "Age"] = "Age_group"
# MRSA harmonisation - here, only oxacillin is used so okay
# Acinetobacter harmonisation - assume baumannii-calcoaceticus species complex is baumannii
df_KEYSTONE_2$Organism[grepl("Acinetobacter baumannii", df_KEYSTONE_2$Organism)] = "Acinetobacter baumannii"
# Salmonella harmonisation - assume all are spp
df_KEYSTONE_2$Organism[grepl("Salmonella", df_KEYSTONE_2$Organism)] = "Salmonella spp"
df_KEYSTONE_2$Organism = as.mo(df_KEYSTONE_2$Organism)


#get antibiotics
if(!(all(is.na(antibiotics_of_interest)))) df_KEYSTONE_2 <- df_KEYSTONE_2[,c(1:6, which(colnames(df_KEYSTONE_2) %in% antibiotics_of_interest))]
#get bacteria
if(!(all(is.na(bacteria_of_interest)))) df_KEYSTONE_2 <- df_KEYSTONE_2 %>% filter(Organism %in% bacteria_of_interest)
#get years
if(length(years_of_interest) != 0) df_KEYSTONE_2 <- df_KEYSTONE_2 %>% filter(`Study Year` %in% years_of_interest)

if(nrow(df_KEYSTONE_2) > 0){
  #format resistances
  # WARNING: several columns use TRUE/FALSE for resistance. Assuming TRUE means
  # resistant, but this leads to some inconsistent results (eg very high S aureus
  # vancomycin resistance). Overall, there are even indications that the interpretation
  # might differ by column, since in any case many vancomycin-res S aureus are MSSA,
  # which seems highly unlikely)
  # I tried to compare to published KEYSTONE estimates, but I should have 100% sus
  # for S aureus vanc in 2019, and that is not the case
  df_KEYSTONE_2 = df_KEYSTONE_2 %>%
    mutate_if(is.logical, as.numeric)
  
  df_KEYSTONE_2 = df_KEYSTONE_2 %>%
    mutate_at(colnames(df_KEYSTONE_2)[
      (sapply(colnames(df_KEYSTONE_2), nchar) == 3) &
        (sapply(sapply(df_KEYSTONE_2, class),"[[",1) == "character")], as.mic)
  
  df_KEYSTONE_2 = df_KEYSTONE_2 %>%
    mutate_at(colnames(df_KEYSTONE_2)[
      (sapply(colnames(df_KEYSTONE_2), nchar) == 3) &
        (sapply(sapply(df_KEYSTONE_2, class),"[[",1) != "mic")], as.character)
  
  df_KEYSTONE_2 = df_KEYSTONE_2 %>%
    mutate_at(colnames(df_KEYSTONE_2)[
      (sapply(colnames(df_KEYSTONE_2), nchar) == 3) &
        (sapply(sapply(df_KEYSTONE_2, class),"[[",1) == "character")], ~replace(., . == "0", "S"))
  
  df_KEYSTONE_2 = df_KEYSTONE_2 %>%
    mutate_at(colnames(df_KEYSTONE_2)[
      (sapply(colnames(df_KEYSTONE_2), nchar) == 3) &
        (sapply(sapply(df_KEYSTONE_2, class),"[[",1) == "character")], ~replace(., . == "1", "R"))
  
  df_KEYSTONE_2 = df_KEYSTONE_2 %>%
    mutate_at(colnames(df_KEYSTONE_2)[
      (sapply(colnames(df_KEYSTONE_2), nchar) == 3) &
        (sapply(sapply(df_KEYSTONE_2, class),"[[",1) != "mic")], as.sir)
  
  df_KEYSTONE_2 = df_KEYSTONE_2[,which(apply(df_KEYSTONE_2, 2, function(x) !(all(is.na(x)))))]
  
  df_KEYSTONE_2b = df_KEYSTONE_2 %>%
    mutate_if(is.mic, as.sir, mo = .$Organism, guideline = "CLSI")
  
  if(any(apply(df_KEYSTONE_2b, 2, function(x) all(is.na(x))))){
    cat("\nRetrying some columns with EUCAST guidelines...\n")
    to_retry = which(apply(df_KEYSTONE_2b, 2, function(x) all(is.na(x))))
    df_KEYSTONE_2b[,to_retry] = df_KEYSTONE_2 %>%
      select(all_of(c(6,to_retry))) %>%
      mutate_if(is.mic, as.sir, mo = .$Organism, guideline = "EUCAST") %>%
      select(-1)
  }
  df_KEYSTONE_2 = df_KEYSTONE_2b
  rm(df_KEYSTONE_2b)
  
}

df_KEYSTONE_2$`Infection Type` = str_to_lower(df_KEYSTONE_2$`Infection Type`)
for(i in 1:nrow(source_index)){
  df_KEYSTONE_2$`Infection Type`[str_which(df_KEYSTONE_2$`Infection Type`, source_index$contains[i])] = source_index$match[i]
}

df_KEYSTONE_2$Age_group[is.na(df_KEYSTONE_2$Age_group)] = -1
df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group > 1900] = df_KEYSTONE_2$`Study Year`[df_KEYSTONE_2$Age_group > 1900]-df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group > 1900]
df_KEYSTONE_2 = df_KEYSTONE_2[!(df_KEYSTONE_2$Age_group > 120),]

df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group %in% c(0:18)] = 0
df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group %in% c(19:64)] = 19
df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group %in% c(65:150)] = 65

df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group %in% c(-1)] = "Unknown"
df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group %in% c(0)] = "0 to 18 Years"
df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group %in% c(19)] = "19 to 64 Years"
df_KEYSTONE_2$Age_group[df_KEYSTONE_2$Age_group %in% c(65)] = "65 and Over"

df_KEYSTONE_2$Gender[is.na(df_KEYSTONE_2$Gender)] = "Unknown"

df_KEYSTONE_2 = class_resistance(df_KEYSTONE_2)

#####

#### SIDERO ####
df_SIDERO_2 <- df_SIDERO[,c(3,6,5,1,7:20)]
# gender and age not available... filling with NA
df_SIDERO_2$Gender = "Unknown"
df_SIDERO_2$Age = "Unknown"
df_SIDERO_2 = df_SIDERO_2[,c(1,19,20,2:18)]
df_SIDERO_2 = df_SIDERO_2 %>% select(-"Meropenem/ Vaborbactam at 8")
colnames(df_SIDERO_2) = sapply(colnames(df_SIDERO_2), function(x) unlist(strsplit(x, "/"))[[1]])
colnames(df_SIDERO_2)[-c(1:6)] = as.ab(colnames(df_SIDERO_2)[-c(1:6)])
# no MRSA here so no harmonisation needed
# Acinetobacter harmonisation - assume spp is baumannii
df_SIDERO_2$`Organism Name`[grepl("Acinetobacter baumannii", df_SIDERO_2$`Organism Name`)] = "Acinetobacter baumannii"
df_SIDERO_2$`Organism Name`[grepl("Acinetobacter sp.", df_SIDERO_2$`Organism Name`)] = "Acinetobacter baumannii"
df_SIDERO_2$`Organism Name` = as.mo(df_SIDERO_2$`Organism Name`)

#get antibiotics
if(!(all(is.na(antibiotics_of_interest)))) df_SIDERO_2 <- df_SIDERO_2[,c(1:6, which(colnames(df_SIDERO_2) %in% antibiotics_of_interest))]
#get bacteria
if(!(all(is.na(bacteria_of_interest)))) df_SIDERO_2 <- df_SIDERO_2 %>% filter(`Organism Name` %in% bacteria_of_interest)
#get years
if(length(years_of_interest) != 0) df_SIDERO_2 <- df_SIDERO_2 %>% filter(`Year Collected` %in% years_of_interest)

if(nrow(df_SIDERO_2) > 0){
  #format resistances
  df_SIDERO_2 = df_SIDERO_2 %>%
    mutate_at(vars(-Country, -Gender, -Age, -`Body Location`, -`Year Collected`, -`Organism Name`), as.numeric) %>%
    mutate_at(vars(-Country, -Gender, -Age, -`Body Location`, -`Year Collected`, -`Organism Name`), round, digits = 3) %>%
    mutate_at(vars(-Country, -Gender, -Age, -`Body Location`, -`Year Collected`, -`Organism Name`), as.mic)
  df_SIDERO_2 = df_SIDERO_2[,which(apply(df_SIDERO_2, 2, function(x) !(all(is.na(x)))))]
  
  df_SIDERO_2b = df_SIDERO_2 %>%
    mutate_if(is.mic, as.sir, mo = .$`Organism Name`, guideline = "CLSI")
  
  if(any(apply(df_SIDERO_2b, 2, function(x) all(is.na(x))))){
    cat("\nRetrying some columns with EUCAST guidelines...\n")
    to_retry = which(apply(df_SIDERO_2b, 2, function(x) all(is.na(x))))
    df_SIDERO_2b[,to_retry] = df_SIDERO_2 %>%
      select(all_of(c(6,to_retry))) %>%
      mutate_if(is.mic, as.sir, mo = .$`Organism Name`, guideline = "EUCAST") %>%
      select(-1)
  }
  df_SIDERO_2 = df_SIDERO_2b
  rm(df_SIDERO_2b)
  
}

df_SIDERO_2$`Body Location` = str_to_lower(df_SIDERO_2$`Body Location`)
for(i in 1:nrow(source_index)){
  df_SIDERO_2$`Body Location`[str_which(df_SIDERO_2$`Body Location`, source_index$contains[i])] = source_index$match[i]
}

df_SIDERO_2 = class_resistance(df_SIDERO_2)

#####

#### GLASS ####
df_GLASS_2 <- df_GLASS[,c("CountryTerritoryArea", "Year", "Specimen", "PathogenName", "AbTargets", "InterpretableAST", "Resistant")]

#get antibiotics
df_GLASS_2$AbTargets = as.ab(df_GLASS_2$AbTargets)
if(!(all(is.na(antibiotics_of_interest)))) df_GLASS_2 <- df_GLASS_2 %>% filter(AbTargets %in% antibiotics_of_interest)
#get bacteria
df_GLASS_2$PathogenName = as.mo(df_GLASS_2$PathogenName)
if(!(all(is.na(bacteria_of_interest)))) df_GLASS_2 <- df_GLASS_2 %>% filter(PathogenName %in% bacteria_of_interest)
#get years
if(length(years_of_interest) != 0) df_GLASS_2 <- df_GLASS_2 %>% filter(Year %in% years_of_interest)

df_GLASS_2 = df_GLASS_2 %>%
  mutate(Age = "Unknown",
         Gender = "Unknown") %>%
  group_by(CountryTerritoryArea, Gender, Age, Specimen, Year, PathogenName, AbTargets) %>%
  summarise(InterpretableAST = sum(InterpretableAST),
            Resistant = sum(Resistant)) %>%
  ungroup() 

df_GLASS_2$Specimen = str_to_lower(df_GLASS_2$Specimen)
for(i in 1:nrow(source_index)){
  df_GLASS_2$Specimen[str_which(df_GLASS_2$Specimen, source_index$contains[i])] = source_index$match[i]
}

################################################################################
####### Compute total and R isolates by country/year/bacteria/antibiotic #######
################################################################################

## Format like in GLASS
final_column_names = c("Country", "Gender", "Age", "Source", "Year", "Pathogen", "Antibiotic", "Total", "Resistant")

#### ATLAS ####

if(nrow(df_ATLAS_2) > 0 & !(all(is.na(df_ATLAS_2)[,-c(1:6)]))){
  df_ATLAS_3 = df_ATLAS_2 %>%
    melt(id.vars = c("Country", "Year", "Species", "Gender", "Age.Group", "Source")) %>%
    mutate(value = as.sir(value)) %>%
    filter(!is.na(value)) %>%
    group_by(Country, Year, Species, Gender, Age.Group, Source, variable, .drop = FALSE) %>%
    summarise(Total = n(), Resistant = count_resistant(value)) %>%
    ungroup %>%
    rename(Pathogen = Species,
           Antibiotic = variable,
           Age = Age.Group) %>%
    mutate(Antibiotic = as.character(Antibiotic)) %>%
    select(all_of(final_column_names))
  
  ## Add dataset name
  df_ATLAS_3$Data <- 'ATLAS'
} else df_ATLAS_3 = data.frame()

#####

##### GEARS ####

if(nrow(df_GEARS_2) > 0 & !(all(is.na(df_GEARS_2)[,-c(1:6)]))){
  df_GEARS_3 = df_GEARS_2 %>%
    melt(id.vars = c("Country", "Year", "Organism", "Gender", "Age", "BodySite")) %>%
    mutate(value = as.sir(value)) %>%
    filter(!is.na(value)) %>%
    group_by(Country, Year, Organism, Gender, Age, BodySite, variable, .drop = FALSE) %>%
    summarise(Total = n(), Resistant = count_resistant(value)) %>%
    ungroup %>%
    rename(Pathogen = Organism,
           Antibiotic = variable,
           Source = BodySite) %>%
    mutate(Antibiotic = as.character(Antibiotic)) %>%
    select(all_of(final_column_names))
  
  ## Add dataset name
  df_GEARS_3$Data <- 'GEARS'
} else df_GEARS_3 = data.frame()

#####

#### KEYSTONE ####

if(nrow(df_KEYSTONE_2) > 0 & !(all(is.na(df_KEYSTONE_2)[,-c(1:6)]))){ 
  df_KEYSTONE_3 = df_KEYSTONE_2 %>%
    melt(id.vars = c("Country", "Study Year", "Organism", "Gender", "Age_group", "Infection Type")) %>%
    mutate(value = as.sir(value)) %>%
    filter(!is.na(value)) %>%
    group_by(Country, `Study Year`, Organism, Gender, Age_group, `Infection Type`, variable, .drop = FALSE) %>%
    summarise(Total = n(), Resistant = count_resistant(value)) %>%
    ungroup %>%
    rename(Pathogen = Organism,
           Antibiotic = variable,
           Year = `Study Year`,
           Age = Age_group,
           Source = `Infection Type`) %>%
    mutate(Antibiotic = as.character(Antibiotic)) %>%
    select(all_of(final_column_names))
  
  ## Add dataset name
  df_KEYSTONE_3$Data <- 'KEYSTONE'
} else df_KEYSTONE_3 = data.frame()

#####

#### SIDERO ####

if(nrow(df_SIDERO_2) > 0 & !(all(is.na(df_SIDERO_2)[,-c(1:6)]))){
  df_SIDERO_3 = df_SIDERO_2 %>%
    melt(id.vars = c("Country", "Year Collected", "Organism Name", "Gender", "Age", "Body Location")) %>%
    mutate(value = as.sir(value)) %>%
    filter(!is.na(value)) %>%
    group_by(Country, `Year Collected`, `Organism Name`, Gender, Age, `Body Location`, variable, .drop = FALSE) %>%
    summarise(Total = n(), Resistant = count_resistant(value)) %>%
    ungroup %>%
    rename(Pathogen = `Organism Name`,
           Antibiotic = variable,
           Year = `Year Collected`,
           Source = `Body Location`) %>%
    mutate(Antibiotic = as.character(Antibiotic)) %>%
    select(all_of(final_column_names))
  
  ## Add dataset name
  df_SIDERO_3$Data <- 'SIDERO'
} else df_SIDERO_3 = data.frame()

#####  

#### GLASS ####

if(nrow(df_GLASS_2) > 0){
  df_GLASS_3 = df_GLASS_2
  colnames(df_GLASS_3) <- final_column_names
  df_GLASS_3$Antibiotic = as.character(df_GLASS_3$Antibiotic)
    
  ## Add dataset name
  df_GLASS_3$Data <- 'GLASS'
} else df_GLASS_3 = data.frame()

#####

################################################################################
############################# Join AMR datasets ################################
################################################################################

## Join all datasets

df_AMR <- rbind(df_GLASS_3,
                df_ATLAS_3,
                df_GEARS_3,
                df_KEYSTONE_3,
                df_SIDERO_3)

if(nrow(df_AMR) == 0){
  cat("The requested year-drug-bug combination does not exist in any dataset!\n")
} else{
  
  for(dataset in c("df_GLASS_3",
                   "df_ATLAS_3",
                   "df_GEARS_3",
                   "df_KEYSTONE_3",
                   "df_SIDERO_3")){
    if(nrow(get(dataset)) == 0) cat("The requested year-drug-bug combination does not exist in",
                                    strsplit(dataset, "_")[[1]][2], "\n")
  }
  
  ## Change countries names so they correspond between datasets
  df_AMR$Country <- str_replace_all(df_AMR$Country,
                                    c("United States" = "USA",
                                      "US$" = "USA",
                                      "United Kingdom of Great Britain and Northern Ireland" = "UK",
                                      "United Kingdom" = "UK",
                                      "Korea, South" = "South Korea",
                                      "Republic of Korea" = "South Korea",
                                      "Russian Federation" = "Russia",
                                      "Lao People's Democratic Republic" = "Laos",
                                      "Slovak Republic" = "Slovakia",
                                      "INDIA" = "India",
                                      "PHILIPPINES" = "Philippines",
                                      "Iran (Islamic Republic of)" = "Iran",
                                      "Czechia" = "Czech Republic"))
  
  ## Add WHO Region ##
  df_GLASS$CountryTerritoryArea =  str_replace_all(df_GLASS$CountryTerritoryArea,
                                                   c("United States" = "USA",
                                                     "US$" = "USA",
                                                     "United Kingdom of Great Britain and Northern Ireland" = "UK",
                                                     "United Kingdom" = "UK",
                                                     "Korea, South" = "South Korea",
                                                     "Republic of Korea" = "South Korea",
                                                     "Russian Federation" = "Russia",
                                                     "Lao People's Democratic Republic" = "Laos",
                                                     "Slovak Republic" = "Slovakia",
                                                     "INDIA" = "India",
                                                     "PHILIPPINES" = "Philippines",
                                                     "Iran (Islamic Republic of)" = "Iran",
                                                     "Czechia" = "Czech Republic"))
  
  df_AMR = df_AMR %>%
    left_join(df_GLASS %>%
                select(CountryTerritoryArea, WHORegionName) %>%
                distinct() %>%
                setNames(c("Country", "Region")))
  
  
  ## Restore nicer bacteria and antibiotic names
  df_AMR$Antibiotic = as.character(df_AMR$Antibiotic)
  df_AMR$Antibiotic[nchar(df_AMR$Antibiotic) == 3] = ab_name(as.ab(df_AMR$Antibiotic[nchar(df_AMR$Antibiotic) == 3]))
  df_AMR$Pathogen = mo_name(as.mo(df_AMR$Pathogen))
  
  ## Summary
  cat("Combined datasets cover", length(unique(df_AMR$Country)), "countries\n",
      "with a total of", sum(df_AMR$Total), "data points\n",
      "(one data point is one test for one year-country-bacteria-antibiotic-dataset combination)\n")
  cat("Countries per datasets:")
  print(df_AMR %>% group_by(Data) %>% summarise(n_countries = length(unique(Country))))
  cat("Data points per datasets:")
  print(df_AMR %>% group_by(Data) %>% summarise(n_datapoints = sum(Total)))
  
  # Save it into csv
  write.csv(x = df_AMR,
            file = here("data", "final_AMR_dataset.csv"),
            row.names = F)
  
}

