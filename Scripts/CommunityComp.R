## Benthic community data ##
## Created by Nyssa Silbiger ###
## Created on 2024-04-03 ##

### load libraries #####
library(here)
library(tidyverse)
library(lubridate)
library(calecopal)
library(patchwork)

## read in data #####
sessile <- read_csv(here("Data","Community_Composition","sessile_data.csv"))
# the functional data
sessile_function <-read_csv(here("Data","Community_Composition","SessileFunctional.csv"))

mobile <- read_csv(here("Data","Community_Composition","mobile_data.csv"))
mobile_function <-read_csv(here("Data","Community_Composition","MobileFunctional.csv"))



## sum across the canopy and understory ###

sessile_total <-sessile %>%
  group_by(plot, season)%>%
  summarise_at(.vars = vars(Silvetia_compressa:Chondracanthus_canaliculatus), .funs = sum)%>% 
  mutate(Total_hairs = rowSums(across(Silvetia_compressa:Chondracanthus_canaliculatus))) %>% # count total crosshairs
  ungroup()%>% # calculate percent cover
  mutate(across(c(Silvetia_compressa:Chondracanthus_canaliculatus),  ~ . /Total_hairs*100))%>%
  rename(BenthicID = plot, Season = season)%>% # keep consistent with other data
  mutate(Season = factor(Season, levels = c("summer2023","fall2023","winter2024")))%>%
  mutate(Benthos = case_when( # add benthos name
    grepl("B", BenthicID) ~ "Barnacles",
    grepl("M", BenthicID) ~"Mussels",
    grepl("P", BenthicID) ~"Surfgrass",
    grepl("R", BenthicID) ~"Rockweed"), 
    .after = BenthicID)


## Make a stacked barplot of all the species across seasons
sessile_total %>%
  group_by(Benthos, Season)%>% # calculate the average for each season
  summarise(across(c(Silvetia_compressa:Chondracanthus_canaliculatus),  ~mean(.,na.rm = TRUE))) %>%
  pivot_longer(cols = c(Silvetia_compressa:Chondracanthus_canaliculatus), names_to = "Species", values_to = "Percent_Cover")%>%
  left_join(sessile_function)%>%
  ggplot(aes(x = Season,y = Percent_Cover, fill = FunctionalGroup ))+
    geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = cal_palette("tidepool", n = 13, type = "continuous"))+
  labs(x = "",
       y = "% Cover",
       color = "")+
  theme_bw()+
  facet_wrap(~Benthos)

ggsave(here("Output","SessileComp.png"), width = 10, height = 6)

############# Mobiles #######

## calculate denisity per m2 (used 0.25m2 quad)
mobile_density <-
  mobile %>%
  rename(BenthicID = plot, Season = season)%>% # keep consistent with other data
  mutate(Season = factor(Season, levels = c("summer2023","fall2023","winter2024")))%>%
  mutate(Benthos = case_when( # add benthos name
    grepl("B", BenthicID) ~ "Barnacles",
    grepl("M", BenthicID) ~"Mussels",
    grepl("P", BenthicID) ~"Surfgrass",
    grepl("R", BenthicID) ~"Rockweed"), 
    .after = BenthicID)%>%
  select(-c(Littorina_sp_raw, Lottia_sp_raw))%>%
  mutate(across(c(Lottia_sp_scaled:Strongylocentrotus_purpuratus),  ~./0.25))
  

mobile_density %>%
  group_by(Benthos, Season)%>% # calculate the average for each season
  summarise(across(c(Lottia_sp_scaled:Strongylocentrotus_purpuratus),  ~mean(.,na.rm = TRUE))) %>%
  pivot_longer(cols = c(Lottia_sp_scaled:Strongylocentrotus_purpuratus), names_to = "Species", values_to = "Density_m2")%>%
  left_join(mobile_function)%>% #join with functional data
  ggplot(aes(x = Season,y = Density_m2, fill = FunctionalGroup ))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = cal_palette("tidepool", n = 9, type = "continuous"))+
  labs(x = "",
       y = "Mean Density (m2)",
       color = "")+
  theme_bw()+
  facet_wrap(~Benthos, scale = "free")

ggsave(here("Output","MobileComp.png"), width = 10, height = 6)
