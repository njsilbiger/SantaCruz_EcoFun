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
library(jtools)
library(interactions)
library(sandwich)
library(patchwork)

## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","Biogeochemistry","TrisCalibrationLog.csv")) %>%
  mutate(TrisCalDate = ymd(TrisCalDate))

pHData<-read_csv(here("Data","Biogeochemistry","pHProbe_Data.csv"))%>%
  mutate(TrisCalDate = ymd(TrisCalDate),
         Sampling_Date = mdy(Sampling_Date))

# Needed for phosphate data


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
#  drop_na(TempInSitu)%>%
  drop_na(mV) %>%
   mutate(pH = pH(Ex=mV,Etris=mVTris,S=Salinity,T=TempInLab))  # calculate pH of the samples using the pH seacarb function


# The TA data is missing from some bottles because of lack of water, but it really doesnt affect the pH calc.
# I am replacing the missing values with 2300 for the pH calculation then converting it back to NA

NoTA<-which(is.na(pHSlope$TA))
pHSlope$TA[NoTA]<-2300

NoPO<-which(is.na(pHSlope$PO4))
pHSlope$PO4[NoPO]<-0

NoTemp<-which(is.na(pHSlope$TempInSitu))
pHSlope$TempInSitu[NoTemp]<-15

#Now calculate pH
pHSlope <-pHSlope%>%
  mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
                            S = Salinity,Pt = PO4, k1k2 = "m10", kf = "dg")) %>%
  # mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
  #                           S = Salinity,Pt = Phosphate_umolL, k1k2 = "m10", kf = "dg")) %>%
  select(!pH) %>% # I only need the in situ pH calculation so remove this
  rename(pH = pH_insitu) %>% # rename it 
  ungroup() %>%
  select(Sampling_Date, Sampling_Time, Day_Night,Benthos, Quad_ID, UniqueID, Day_Night, Salinity, pH, TempInSitu, TA,Processing_DateTime,NO3, PO4, SIL, NO2, NH4, Notes) # keep what I want

pHSlope$TA[NoTA]<-NA # make TA na again for the missing values
pHSlope$PO4[NoPO]<-NA # make PO na again for the missing values
pHSlope$TempInSitu[NoTemp]<-NA # make TA na again for the missing values

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
  mutate(Sampling_Time = as.character(Sampling_Time)) %>%
  arrange(row_number() < match('06:00:00', Sampling_Time)) %>%
  mutate(Sampling_Time = factor(Sampling_Time, unique(Sampling_Time))) %>%
  ggplot(aes(x = Sampling_Time, y = pH))+
  geom_point()+
 # geom_smooth()+
  facet_wrap(~Benthos)+
  theme(axis.text.x = element_text(hjust = 1, angle=90,
                                   face="bold", size=12), 
        axis.text.y = element_text( face="bold", size=12),
        strip.text = element_text(size=12, face="bold"))

p1_pH<-pHSlope %>% 
  filter(Benthos != "Open Ocean")%>%
  mutate(Sampling_Time_Fake = ymd_hms(paste("2023-10-10",as.character(Sampling_Time)))) %>%
  mutate(Sampling_Time_Fake  = if_else(hour(Sampling_Time_Fake) > 6, Sampling_Time_Fake - days(1), Sampling_Time_Fake)) %>% # so that we can reorder the plot to start with 6am
  ggplot(aes(x = Sampling_Time_Fake, y = pH, color = Benthos))+
  geom_point()+
  geom_smooth(method = "lm", formula = "y~poly(x,2)")+
  scale_x_datetime(date_labels = "%H:%M",
                   limits = c(ymd_hms("2023-10-09 05:00:00", ymd_hms("2023-10-10 04:00:00"))) )+
  scale_color_manual(values = cal_palette("chaparral1"))+  facet_wrap(~Benthos, nrow=1)+
  labs(x = "",
       y = "pH")+
  theme_bw() +
  theme(legend.position = "none",
                   axis.title = element_text(size = 16),
                   axis.text = element_text(size = 14),
                   strip.background = element_blank(),
                   strip.text = element_text(size = 16),
        axis.text.x = element_blank()
  )


ggsave(here("Output","pHcontinuous.png"), width = 8, height = 4)
# 

pHSlope %>%
  ggplot(aes(x = pH, fill = Benthos))+
  geom_density(alpha = 0.5)
#  geom_smooth()

## pull out open ocean data and average it by date
ocean<-pHSlope%>% 
  filter(Benthos == "Open Ocean") %>%
  mutate(Sampling_DateTime = ymd_hms(paste(Sampling_Date, Sampling_Time))) %>%
  group_by(Sampling_DateTime) %>%
  summarise(OceanpH = mean(pH))

pHSlope<-pHSlope %>% 
  mutate(Sampling_DateTime = ymd_hms(paste(Sampling_Date, Sampling_Time))) %>%
  left_join(ocean)%>%
 # group_by(Sampling_DateTime) %>%
  mutate(deltapH = pH - OceanpH) # subtract from average ocean sample

### add a new column for time groupings
pHSlope<-pHSlope %>%
  mutate( Samplinghour = as.numeric(hour(Sampling_Time)),
    group_time = 
           case_when(Samplinghour >=6 & Samplinghour < 12 ~ "Morning",
                     Samplinghour >=12 &Samplinghour< 18 ~ "Afternoon",
                     Samplinghour >=18 & Samplinghour < 24 ~ "Evening",
                     Samplinghour >= 0 & Samplinghour< 6 ~ "Night"
                     )) %>%
  mutate(group_time = factor(group_time, levels = c("Morning","Afternoon","Evening","Night")))
#   mutate(deltapH = pH - pH[Benthos == "Open Ocean"]) %>%
#   ungroup() 

pHSlope%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = Benthos, y = deltapH, fill = group_time))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_boxplot()+
  geom_jitter(position = position_dodge(width = .8))+
  
  labs(y = "pH (Difference from Ocean Sample)",
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

### do it continunous
p2_pH<-pHSlope %>% 
  mutate(Sampling_Time_Fake = ymd_hms(paste("2023-10-10",as.character(Sampling_Time)))) %>%
  mutate(Sampling_Time_Fake  = if_else(hour(Sampling_Time_Fake) >= 5, Sampling_Time_Fake - days(1), Sampling_Time_Fake)) %>% # so that we can reorder the plot to start with 6am
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = Sampling_Time_Fake, y = deltapH, color = Benthos))+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_smooth(method = "lm", formula = "y~poly(x,2)")+
  scale_x_datetime(date_labels = "%H:%M",
                   limits = c(ymd_hms("2023-10-09 05:00:00", ymd_hms("2023-10-10 04:00:00"))) )+
  scale_color_manual(values = cal_palette("chaparral1"))+
  facet_wrap(~Benthos, nrow=1)+
  labs(x = "",
       y = expression(paste(Delta,"pH from Ocean")))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_blank()
        #strip.text = element_text(size = 16)
        )

p1_pH/p2_pH

ggsave(here("Output","pHcomposite.png"), width = 10, height = 6)


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

pHSlope<-pHSlope %>%
  drop_na(TempInSitu)


### Calculate DIC from seacarb ####
AllCO2<-carb(8, pHSlope$pH, pHSlope$TA/1000000, S=pHSlope$Salinity, T=pHSlope$TempInSitu, Patm=1, P=0, Pt=0, Sit=0,
             k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential", 
             warn="y", eos="eos80")

AllCO2 <- AllCO2 %>%
  mutate(ALK = ALK*1000000,
         CO2 = CO2*1000000,
         CO3 = CO3*1000000,
         DIC = DIC*1000000,
         HCO3 = HCO3*1000000) %>% # convert everything back to umol %>%
select(DIC, pCO2 = pCO2insitu, CO2, CO3, HCO3, OmegaCalcite, OmegaAragonite) 

AllCO2 <- pHSlope %>%
  bind_cols(AllCO2) %>% # bring together with the original data
  filter(deltapH < 0.14) %>% # There is a huge mussel outlier in the pH data
  mutate(TA_norm = TA*Salinity/33,
         DIC_norm = DIC*Salinity/33,# salinity normalize
         TA_DIC = TA/DIC,
         Month = month(Sampling_Date)) %>% # TA divided by DIC
  mutate(Season = case_when(Month %in% c(7,8)~"Summer",
                            Month %in% c(10,11)~ "Fall"))

# make the TA vs DIC plots
p_slopes<-AllCO2 %>%
 # filter(Benthos !="Open Ocean") %>%
  filter(TA <2300) %>%
ggplot(aes(x = DIC, y = TA, color = Season))+
  geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  labs(x = expression(paste("DIC (",mu,"mol kg"^-1,")")),
       y = expression(paste("TA (",mu,"mol kg"^-1,")")))+
  scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(#legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))+
  facet_wrap(~Benthos)

ggsave(here("Output","TADIC.png"), width = 6, height = 4)
  #facet_wrap(~Benthos)

AllCO2 %>%
  # filter(Benthos !="Open Ocean") %>%
  ggplot(aes(x = DIC, y = TA, color = Season))+
  geom_point()+
  geom_smooth(method = "lm",se = FALSE)

# run an ancova to see if the slopes are different
TADICmod<-lm(TA~DIC*Benthos, data = AllCO2 %>%filter(Benthos !="Open Ocean"))
anova(TADICmod)
summary(TADICmod)

# calculate marginal effects
ss <- sim_slopes(TADICmod, pred = DIC, modx = Benthos, johnson_neyman = FALSE)
plot(ss)

# extract the individual slopes
slopes<-tibble(ss$slopes)
colnames(slopes)[1]<-"Benthos"

# make the plot
P_estimate<-slopes %>%
  ggplot(aes(y = Benthos, x = Est., color = Benthos))+
  geom_point(size = 3)+
  geom_errorbarh(aes(xmin = Est. - S.E.,xmax = Est. + S.E.  ), height = 0.01)+
  scale_color_manual(values = cal_palette("chaparral1"))+
  labs(x = "TA/DIC Slopes",
       y = "")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
# Type 2 linear regression

# bring the TA DIC plots together
p_slopes+P_estimate

ggsave(here("Output","TADICComposite.png"), width = 8, height = 4)

# plot 1 by 1 and extract the confidence intervals
type2<- lmodel2(DIC ~ TA, data=AllCO2 %>% filter(Benthos == "Rockweed"), nperm = 999)
plot(type2, "MA")
type2


##### Plot average aragonit and calc sat state #####

AllCO2 %>%
  ggplot(aes(x = Benthos, y = TA/DIC))+
  geom_point()

TADICmod<-lm(TA_DIC~Benthos*Season, data = AllCO2)
anova(TADICmod)
summary(TADICmod)

AllCO2<-AllCO2  %>%
  mutate(TA_DIC = TA/DIC) %>%
  group_by(Sampling_DateTime) %>%
  mutate(deltaTADIC = TA_DIC - TA_DIC[Benthos == "Open Ocean"],
         deltaNO2 = NO2 - NO2[Benthos == "Open Ocean"],
         deltaNO3 = NO3 - NO3[Benthos == "Open Ocean"],
         deltaPO4 = PO4 - PO4[Benthos == "Open Ocean"],
         deltaNH4 = NH4 - NH4[Benthos == "Open Ocean"]
        # deltaSIL = SIL - SIL[Benthos == "Open Ocean"]
         ) %>% # difference from open ocean
  ungroup() 


AllCO2 %>% 
  filter(Benthos != "Open Ocean") %>%
  ggplot(aes(x = Benthos, y = deltaTADIC))+
  geom_boxplot()

modTADICdelta<-lm(deltaTA_DIC~Benthos, data = AllCO2 %>% 
                    filter(Benthos != "Open Ocean") )
anova(modTADICdelta)
summary(modTADICdelta)

AllCO2 %>% 
  filter(Benthos != "Open Ocean") %>%
  ggplot(aes(x = Benthos, y = deltaPO4))+
  geom_boxplot()

AllCO2 %>% 
 # filter(Benthos != "Open Ocean") %>%
  ggplot(aes(x = Benthos, y = PO4))+
  geom_boxplot()

AllCO2 %>%
  ggplot(aes(x = TempInSitu, y = deltapH, color = Benthos))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Benthos)

## delta pH strongly correlated with temp
mod_pH_Temp<-lm(deltapH~TempInSitu*Benthos, data = AllCO2 %>% filter(Benthos != "Open Ocean"))
anova(mod_pH_Temp)
summary(mod_pH_Temp)

### raw pH NOT correlated with temp, meaning that the temp is causing differences in the CHANGE in pH which is due to metabolism
AllCO2 %>%
  ggplot(aes(x = TempInSitu, y = pH, color = Benthos))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Benthos)

