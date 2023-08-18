## process pH for Orion
# Created by Nyssa Silbiger
# Edited on 07/27/2021

library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)
library(calecopal)
library(ggridges)

## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","Biogeochemistry","TrisCalibrationLog.csv")) %>%
  mutate(TrisCalDate = ymd(TrisCalDate))

pHData<-read_csv(here("Data","Biogeochemistry","pHProbe_Data.csv"))%>%
  mutate(TrisCalDate = ymd(TrisCalDate),
         Sampling_Date = mdy(Sampling_Date))

# Needed for phosphate data
# NutData<-read_csv(here("Data","March2022","Nutrients","Nutrients_watersampling_Mar22.csv")) %>%
#   select(CowTagID, SeepCode, Day_Night, Tide, Date, Phosphate_umolL, NN_umolL=Nitrite_umolL, Ammonia_umolL, Silicate_umolL) %>%
#   mutate(Date = mdy(Date))
# 
# pHData<-left_join(pHData,NutData) 

NoTemp<-which(is.na(pHData$TempInSitu))

pHData$TempInSitu[NoTemp]<-19

## take the mV calibration files by each date and use them to calculate pH
pHSlope<-pHcalib %>%
  filter(HOBO_Orion =="Orion") %>% # extract only the orion data
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = data))) %>% # linear regression of mV and temp of the tris
  summarise(broom::tidy(fitpH)) %>% # make the output tidy
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  right_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = TempInLab*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  drop_na(TempInSitu)%>%
  drop_na(mV) %>%
   mutate(pH = pH(Ex=mV,Etris=mVTris,S=Salinity,T=TempInLab))  # calculate pH of the samples using the pH seacarb function


# The TA data is missing from some bottles because of lack of water, but it really doesnt affect the pH calc.
# I am replacing the missing values with 2300 for the pH calculation then converting it back to NA

NoTA<-which(is.na(pHSlope$TA))

pHSlope$TA[NoTA]<-2300

NoPO<-which(is.na(pHSlope$Phosphate_umolL))
pHSlope$Phosphate_umolL[NoPO]<-0

# If the salinity is missing from the TA measurements take the salinity from the field
pHSlope$Salinity<-ifelse(is.na(pHSlope$Salinity), pHSlope$Salinity_In_Lab, pHSlope$Salinity)


#Now calculate pH
pHSlope <-pHSlope%>%
  mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
                            S = Salinity,Pt = 0.1, k1k2 = "m10", kf = "dg")) %>%
  # mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
  #                           S = Salinity,Pt = Phosphate_umolL, k1k2 = "m10", kf = "dg")) %>%
  select(!pH) %>% # I only need the in situ pH calculation so remove this
  rename(pH = pH_insitu) %>% # rename it 
  ungroup() %>%
  select(Sampling_Date, Sampling_Time, Day_Night,Benthos, Quad_ID, UniqueID, Day_Night, Salinity, pH, TempInSitu, TA,Processing_DateTime, Notes) # keep what I want

pHSlope$TA[NoTA]<-NA # make TA na again for the missing values

  #select(Date, CowTagID,Tide, Day_Night, SamplingTime,Salinity,pH, pH_insitu, TempInSitu) ## need to calculate pH insi then it is done

## write the data
write_csv(x = pHSlope, file = here("Data","Biogeochemistry","pHProbe_Data_calculated.csv"))



summ<-pHSlope %>% 
  group_by(Benthos,Quad_ID) %>%
  summarise(deltapH = abs(pH[Day_Night == "Day"] - pH[Day_Night == "Night"])) 


summ %>%  
ggplot(aes(x = Benthos, y = deltapH, color = Benthos))+
  geom_boxplot()+
  geom_jitter(width = 0.1)

anova(lm(deltapH~Benthos, data = summ))


pHSlope %>%
  ggplot(aes(x = Benthos, y = pH, color = Day_Night, shape = factor(Sampling_Date)))+
  geom_boxplot(width = .1)+
  facet_grid(~factor(Sampling_Date))

pHSlope %>%
  ggplot(aes(x = Sampling_Time, y = pH))+
  geom_point()+
  facet_wrap(~Benthos)
#  geom_smooth()


pHSlope %>%
  ggplot(aes(x = pH, fill = Benthos))+
  geom_density(alpha = 0.5)+
  facet_wrap(~Day_Night)
#  geom_smooth()

pHSlope<-pHSlope %>% 
  mutate(Sampling_DateTime = ymd_hms(paste(Sampling_Date, Sampling_Time))) %>%
  group_by(Sampling_DateTime) %>%
  mutate(deltapH = pH - pH[Benthos == "Open Ocean"]) %>%
  ungroup() 

pHSlope%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = Benthos, y = deltapH, fill = Day_Night))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_boxplot()+
  geom_jitter(position = position_dodge(width = .8))+
  
  labs(y = "pH (Difference from Open Ocean Sample)",
       x = "",
       fill = "Day/Night")+
  scale_fill_manual(values = cal_palette("chaparral1"))+
  #scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

ggsave(here("Output","pHdifference.png"), width = 8, height = 6)

## sample plot but means and errorbars instead
pHSlope%>%
  filter(Benthos != "Open Ocean")%>%
  group_by(Benthos, Day_Night)%>%
  summarise(mean_pHdiff = mean(deltapH, na.rm = TRUE),
            se_pHdiff = sd(deltapH, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = Benthos, y = mean_pHdiff, color = Day_Night))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 4)+
  geom_jitter(data = pHSlope%>%
                filter(Benthos != "Open Ocean"), aes(x = Benthos, y = deltapH),alpha = 0.3, width = 0.1)+
  geom_errorbar(aes(x = Benthos, ymin = mean_pHdiff-se_pHdiff, ymax =mean_pHdiff+se_pHdiff), width = 0.01, size = 1.2)+
  labs(y = "pH (Difference from Open Ocean Sample)",
       x = "",
       color = "Day/Night")+
  scale_color_manual(values = cal_palette("chaparral1"))+
  #scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(here("Output","pHdifference_means.png"), width = 8, height = 6)


pHSlope%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = deltapH, y = Benthos, fill = Benthos))+
  geom_vline(xintercept = 0)+
  geom_density_ridges(alpha = 0.5)+
  scale_fill_manual(values = cal_palette("tidepool"))+
 # scale_fill_manual(values = cal_palette("chaparral1"))+
  labs(x = "pH (Difference from Open Ocean Sample)",
       y = "")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

  
