library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(data.table)

col_pal = brewer.pal(6, "Dark2")
  
df_AMR = read.csv(here::here("data", "final_AMR_dataset.csv")) %>%
  filter(Data != "GLASS") %>%
  filter(!(Source %in% c("HEENT", "Reproduction", "Other"))) %>%
  group_by(Country, Region, Year, Source, Pathogen, Antibiotic) %>%
  summarise(Resistant = sum(Resistant),
            Total = sum(Total)) %>%
  mutate(prop = Resistant/Total) %>%
  filter(Total>=10) %>%
  filter(!(is.nan(prop))) %>%
  filter(Antibiotic %in% c("Aminoglycoside", "Aminopenicillin", "Carbapenem", "Cephalosporin_3rd", "Cephalosporin_4th", "Fluoroquinolone", "Sulfonamide_c", "Polypeptide")) %>%
  mutate(Antibiotic = replace(Antibiotic, Antibiotic=="Sulfonamide_c", "Sulfonamide"))

df_AMR_agg = df_AMR %>%
  group_by(Region, Year, Source, Pathogen, Antibiotic) %>%
  summarise(Resistant = sum(Resistant),
            Total = sum(Total)) %>%
  mutate(prop = Resistant/Total) %>%
  ungroup


df_AMR_all = df_AMR %>%
  group_by(Source, Pathogen, Antibiotic) %>%
  summarise(prop=mean(prop))

pa = ggplot(df_AMR_all, aes(x=Source, y=Antibiotic, fill = prop)) +                           
  geom_tile(color = "white")+
  scale_fill_gradient2(low = 'cornsilk1',  mid = 'gold', high = 'darkred', #midpoint = 0.5, 
                       name = 'Average resistance proportion') + 
  theme_bw() +
  theme(legend.position="right") +
  xlab("Sampling sources") +
  ylab("Antibiotic class") +
  facet_grid(rows=vars(Pathogen)) +
  geom_text(aes(label = round(prop, digits=2)), 
            color = "white", fontface='bold', size = 3) +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 11, angle=45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12))

abx = "Fluoroquinolone"
bac = "Escherichia coli"

pb = ggplot(df_AMR_agg %>% filter(Pathogen==bac & Antibiotic==abx & Region %in% c("European Region", "African Region"))) +
  geom_line(aes(Year, prop, colour = Region), linewidth=1) +
  geom_point(data = df_AMR %>% 
               filter(Pathogen==bac &
                        Antibiotic==abx &
                        Region == "European Region"),
             aes(Year, prop, colour = Region, group = Country), alpha = 0.05, size=2) +
  geom_line(data = df_AMR %>%
              filter(Pathogen==bac &
                       Antibiotic==abx &
                       Region == "European Region"),
            aes(Year, prop, colour = Region, group = Country), alpha = 0.05, linewidth=1) +
  geom_point(data = df_AMR %>% 
              filter(Pathogen==bac &
                       Antibiotic==abx &
                       Region == "African Region"),
            aes(Year, prop, colour = Region, group = Country), alpha = 0.2, size=2) +
  geom_line(data = df_AMR %>% 
               filter(Pathogen==bac &
                        Antibiotic==abx &
                        Region == "African Region"),
             aes(Year, prop, colour = Region, group = Country), alpha = 0.2, linewidth=1) +
  facet_grid(cols = vars(Source)) +
  scale_color_discrete(type=col_pal) +
  theme_bw() +
  theme(strip.text.y=element_text(size=11),
        axis.text.x=element_text(angle=45,hjust=1,size=11),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x=element_text(size=11),
        legend.text=element_text(size=11),
        legend.title = element_text(size=12),
        legend.position = "bottom") +
  labs(x = "Year", y = "Resistance proportion", colour="WHO Region:")


plot_grid(pa+theme(legend.position="none"), pb, rel_widths = c(0.6,1),
          labels=c("a)", "b)"), hjust=0)

ggsave(here("figures", "fig1.png"), height=5, width = 11, bg = "white")
