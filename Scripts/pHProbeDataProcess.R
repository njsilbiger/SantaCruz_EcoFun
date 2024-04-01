## process pH for Orion
# Created by Nyssa Silbiger
# Edited on 03/31/2024

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
library(ggtext)

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
  reframe(broom::tidy(fitpH)) %>% # make the output tidy
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
  reframe(deltapH = abs(pH[Day_Night == "Day"] - pH[Day_Night == "Night"])) 


summ %>%  
ggplot(aes(x = Benthos, y = deltapH, color = Benthos))+
  geom_boxplot()+
  geom_jitter(width = 0.1)

anova(lm(deltapH~Benthos, data = summ))


pHSlope %>%
  ggplot(aes(x = Benthos, y = pH, color = Day_Night #shape = factor(Sampling_Date)
             ))+
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
# ocean<-pHSlope%>% 
#   filter(Benthos == "Open Ocean") %>%
#   mutate(Sampling_DateTime = ymd_hms(paste(Sampling_Date, Sampling_Time))) %>%
#   group_by(Sampling_DateTime, Day_Night) %>%
#   summarise(OceanpH = mean(pH))

ocean<-pHSlope%>% 
  filter(Benthos == "Open Ocean") %>%
 # mutate(Sampling_DateTime = ymd_hms(paste(Sampling_Date, Sampling_Time))) %>%
  group_by(Sampling_Date, Day_Night) %>%
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
  mutate(group_time = factor(group_time, levels = c("Morning","Afternoon","Evening","Night"))) # TA divided by DIC

#   mutate(deltapH = pH - pH[Benthos == "Open Ocean"]) %>%
#   ungroup() 

pHSlope%>%
  mutate( Month = month(Sampling_Date)) %>% # TA divided by DIC
  mutate(Season = case_when(Month %in% c(7,8)~"Summer",
                            Month %in% c(10,11)~ "Fall")) %>%
  filter(Benthos != "Open Ocean", Season == "Fall")%>%
  ggplot(aes(x = Benthos, y = deltapH, fill = Day_Night))+
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


pHSlope%>%
  mutate( Month = month(Sampling_Date),
          Benthos = factor(Benthos, levels = c("Open Ocean","Barnacle","Mussel","Rockweed","Surfgrass"))) %>% # TA divided by DIC
  mutate(Season = case_when(Month %in% c(7,8)~"Summer",
                            Month %in% c(10,11)~ "Fall")) %>%
  filter(Season == "Fall",)%>%
  ggplot(aes(x = Benthos, y = pH, fill = Day_Night))+
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


meanocean<-pHSlope %>%
  mutate(Benthos = factor(Benthos, levels = c("Open Ocean","Barnacle","Mussel","Rockweed","Surfgrass")),
    Month = month(Sampling_Date))%>%
  mutate(Season = case_when(Month %in% c(7,8)~"Summer",
                            Month %in% c(10,11)~ "Fall")) %>%
  filter(
         Season == "Fall",
         PO4<1,
         NH4 <4)%>%
  group_by(Benthos, Day_Night)%>%
  summarise(pH_mean = mean(pH, na.rm = TRUE),
            pH_se = sd(pH, na.rm = TRUE)/sqrt(n()),
            NO3_mean = mean(NO3, na.rm = TRUE),
            NO3_se = sd(NO3, na.rm = TRUE)/sqrt(n()),
            PO4_mean = mean(PO4, na.rm = TRUE),
            PO4_se = sd(PO4, na.rm = TRUE)/sqrt(n()),
            NH4_mean = mean(NH4, na.rm = TRUE),
            NH4_se = sd(NH4, na.rm = TRUE)/sqrt(n()),
            SIL_mean = mean(SIL, na.rm = TRUE),
            SIL_se = sd(SIL, na.rm = TRUE)/sqrt(n()),
            TA_mean = mean(TA, na.rm = TRUE),
            TA_se = sd(TA, na.rm = TRUE)/sqrt(n())
            
  )


## sample plot but means and errorbars instead - shaded reagoin is mean and error for the ocean sample
pH_plotmean<-meanocean%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = as.numeric(Benthos), y = pH_mean, color = Day_Night))+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 7.90-0.00829 ,
                ymax = 7.90+0.00829) ,
                fill = "#DCC27A", alpha = 0.1, show.legend = FALSE)+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 8.09-0.0109 ,
                ymax = 8.09+0.0109) ,
            fill = "#B0B9BE", alpha = 0.1, show.legend = FALSE)+
  geom_hline(yintercept = 7.9, lty = 2)+
  
  geom_hline(yintercept = 8.09, lty = 2)+
  geom_point(size = 4)+
  geom_errorbar(aes(x = as.numeric(Benthos), ymin = pH_mean-pH_se, ymax =pH_mean+pH_se), width = 0.01, size = 1.2)+
  labs(y = "pH ",
       x = "",
       color = "")+
  scale_color_manual(values = cal_palette("chaparral1"))+
  scale_x_continuous(breaks=c(2,3,4,5),
                   labels=c("Barnacle", "Mussel", "Rockweed","Surfgrass"))+
  #scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        )


#ggsave(filename = here("Output","pH_plotmean.png"),pH_plotmean, width = 6, height = 4)


NO3_plotmean<-meanocean%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = as.numeric(Benthos), y = NO3_mean, color = Day_Night))+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 1.44-0.320  ,
                ymax = 1.44+0.320 ) ,
            fill = "#DCC27A", alpha = 0.1, show.legend = FALSE)+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 2.18 -0.225  ,
                ymax = 2.18 +0.225 ) ,
            fill = "#B0B9BE", alpha = 0.1, show.legend = FALSE)+
  geom_hline(yintercept = 1.44, lty = 2)+
  
  geom_hline(yintercept = 2.18, lty = 2)+
  geom_point(size = 4)+
  geom_errorbar(aes(x = as.numeric(Benthos), ymin = NO3_mean-NO3_se, ymax =NO3_mean+NO3_se), width = 0.01, size = 1.2)+
  labs(y = "NO3 ",
       x = "",
       color = "")+
  scale_color_manual(values = cal_palette("chaparral1"))+
  scale_x_continuous(breaks=c(2,3,4,5),
                     labels=c("Barnacle", "Mussel", "Rockweed","Surfgrass"))+
  #scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
  )


PO4_plotmean<-meanocean%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = as.numeric(Benthos), y = PO4_mean, color = Day_Night))+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 0.375 -0.0184   ,
                ymax = 0.375 +0.0184  ) ,
            fill = "#DCC27A", alpha = 0.1, show.legend = FALSE)+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 0.438 -0.0261   ,
                ymax = 0.438 +0.0261  ) ,
            fill = "#B0B9BE", alpha = 0.1, show.legend = FALSE)+
  geom_hline(yintercept = 0.375, lty = 2)+
  
  geom_hline(yintercept = 0.438, lty = 2)+
  geom_point(size = 4)+
  geom_errorbar(aes(x = as.numeric(Benthos), ymin = PO4_mean-PO4_se, ymax =PO4_mean+PO4_se), width = 0.01, size = 1.2)+
  labs(y = "PO4 ",
       x = "",
       color = "")+
  scale_color_manual(values = cal_palette("chaparral1"))+
  scale_x_continuous(breaks=c(2,3,4,5),
                     labels=c("Barnacle", "Mussel", "Rockweed","Surfgrass"))+
  #scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
  )

NH4_plotmean<-meanocean%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = as.numeric(Benthos), y = NH4_mean, color = Day_Night))+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 1.21-  0.147   ,
                ymax = 1.21 + 0.147  ) ,
            fill = "#DCC27A", alpha = 0.1, show.legend = FALSE)+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 1.23 - 0.286   ,
                ymax = 1.23 + 0.286  ) ,
            fill = "#B0B9BE", alpha = 0.1, show.legend = FALSE)+
  geom_hline(yintercept = 1.21, lty = 2)+
  
  geom_hline(yintercept = 1.23, lty = 2)+
  geom_point(size = 4)+
  geom_errorbar(aes(x = as.numeric(Benthos), ymin = NH4_mean-NH4_se, ymax =NH4_mean+NH4_se), width = 0.01, size = 1.2)+
  labs(y = "NH4 ",
       x = "",
       color = "")+
  scale_color_manual(values = cal_palette("chaparral1"))+
  scale_x_continuous(breaks=c(2,3,4,5),
                     labels=c("Barnacle", "Mussel", "Rockweed","Surfgrass"))+
  #scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
  )

SIL_plotmean<-meanocean%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = as.numeric(Benthos), y = SIL_mean, color = Day_Night))+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 4.78 - 0.343   ,
                ymax = 4.78 + 0.343  ) ,
            fill = "#DCC27A", alpha = 0.1, show.legend = FALSE)+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 5.28 - 0.355   ,
                ymax = 5.28 + 0.355  ) ,
            fill = "#B0B9BE", alpha = 0.1, show.legend = FALSE)+
  geom_hline(yintercept = 4.78, lty = 2)+
  
  geom_hline(yintercept = 5.28, lty = 2)+
  geom_point(size = 4)+
  geom_errorbar(aes(x = as.numeric(Benthos), ymin = SIL_mean-SIL_se, ymax =SIL_mean+SIL_se), width = 0.01, size = 1.2)+
  labs(y = "Silicate ",
       x = "",
       color = "")+
  scale_color_manual(values = cal_palette("chaparral1"))+
  scale_x_continuous(breaks=c(2,3,4,5),
                     labels=c("Barnacle", "Mussel", "Rockweed","Surfgrass"))+
  #scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
  )

TA_plotmean<-meanocean%>%
  filter(Benthos != "Open Ocean")%>%
  ggplot(aes(x = as.numeric(Benthos), y = TA_mean, color = Day_Night))+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 2222. - 1.76   ,
                ymax = 2222.  +1.76  ) ,
            fill = "#DCC27A", alpha = 0.1, show.legend = FALSE)+
  geom_rect(aes(xmin = 1.9 , 
                xmax = 5.1 ,
                ymin = 2204- 4.61   ,
                ymax = 2204+  4.61 ) ,
            fill = "#B0B9BE", alpha = 0.1, show.legend = FALSE)+
  geom_hline(yintercept = 2222., lty = 2)+
  
  geom_hline(yintercept = 2204, lty = 2)+
  geom_point(size = 4)+
  geom_errorbar(aes(x = as.numeric(Benthos), ymin = TA_mean-TA_se, ymax =TA_mean+TA_se), width = 0.01, size = 1.2)+
  labs(y = "TA ",
       x = "",
       color = "")+
  scale_color_manual(values = cal_palette("chaparral1"))+
  scale_x_continuous(breaks=c(2,3,4,5),
                     labels=c("Barnacle", "Mussel", "Rockweed","Surfgrass"))+
  #scale_color_manual(values = cal_palette("chaparral1"))+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
  )


(pH_plotmean+TA_plotmean)/(NO3_plotmean+NH4_plotmean)/(SIL_plotmean+ PO4_plotmean) +plot_layout(guides = 'collect') +plot_annotation(tag_levels = "a")

ggsave(filename = here("Output","compositeChem.png"), width = 10, height = 10)

## sample plot but means and errorbars instead
pHSlope%>%
  filter(Benthos != "Open Ocean")%>%
  group_by(Benthos, group_time)%>%
  summarise(mean_pHdiff = mean(deltapH, na.rm = TRUE),
            se_pHdiff = sd(deltapH, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = Benthos, y = mean_pHdiff, color = group_time))+
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

# pHSlope<-pHSlope %>%
#   drop_na(TempInSitu)


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
TADICmod<-lm(TA~DIC*Benthos*Season, data = AllCO2 %>%filter(Benthos !="Open Ocean") )
anova(TADICmod)
summary(TADICmod)

# calculate marginal effects - add in season interaction here... 
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
#type2<- lmodel2(DIC ~ TA, data=AllCO2 %>% filter(Benthos == "Rockweed"), nperm = 999)
#plot(type2, "MA")
#type2


##### Plot average aragonit and calc sat state #####

AllCO2 %>%
  ggplot(aes(x = Benthos, y = TA/DIC))+
  geom_point()

TADICmod<-lm(TA_DIC~Benthos*Season, data = AllCO2)
anova(TADICmod)
summary(TADICmod)

AllCO2<-AllCO2  %>%
  mutate(TA_DIC = TA/DIC) %>%
#  group_by(Sampling_DateTime) %>%
  group_by(Sampling_Date, Day_Night) %>%
  mutate(deltaTADIC = TA_DIC - TA_DIC[Benthos == "Open Ocean"],
         deltaNO2 = NO2 - NO2[Benthos == "Open Ocean"],
         deltaNO3 = NO3 - NO3[Benthos == "Open Ocean"],
         deltaPO4 = PO4 - PO4[Benthos == "Open Ocean"],
         deltaNH4 = NH4 - NH4[Benthos == "Open Ocean"],
         deltaSIL = SIL - SIL[Benthos == "Open Ocean"]
         ) %>% # difference from open ocean
  ungroup() 


AllCO2 %>% 
  filter(Benthos != "Open Ocean") %>%
  ggplot(aes(x = Benthos, y = deltaTADIC))+
  geom_boxplot()

# modTADICdelta<-lm(deltaTA_DIC~Benthos, data = AllCO2 %>% 
#                     filter(Benthos != "Open Ocean") )
# anova(modTADICdelta)
# summary(modTADICdelta)

AllCO2 %>% 
  filter(Benthos != "Open Ocean",
         deltaPO4<4,
         deltaPO4 > -4) %>% # remove one huge outlier
  ggplot(aes(x = Benthos, y = deltaPO4, fill = Day_Night))+
  geom_boxplot()+
  #geom_jitter(posi)+
  facet_wrap(~Season, scales = "free")

AllCO2 %>% 
 # filter(Benthos != "Open Ocean") %>%
  ggplot(aes(x = Benthos, y = PO4, fill = Season))+
  geom_boxplot()

AllCO2 %>%
  ggplot(aes(x = TempInSitu, y = deltapH, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Benthos)

#ANCOVA temp to pH
Muss<-lm(deltapH ~ TempInSitu*Season, data = AllCO2 %>% filter(Benthos == "Mussel"))
anova(Muss)

Rocks<-lm(deltapH ~ TempInSitu*Season, data = AllCO2 %>% filter(Benthos == "Rockweed"))
anova(Rocks)

Barn<-lm(deltapH ~ TempInSitu*Season, data = AllCO2 %>% filter(Benthos == "Barnacle"))
anova(Barn)

Surf<-lm(deltapH ~ TempInSitu*Season, data = AllCO2 %>% filter(Benthos == "Surfgrass"))
anova(Surf)

Oce<-lm(deltapH ~ TempInSitu*Season, data = AllCO2 %>% filter(Benthos == "Open Ocean"))
anova(Oce)

AllCO2 %>%
  ggplot(aes(x = TempInSitu, y = pH, color = Season))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Benthos)

## delta pH strongly correlated with temp
mod_pH_Temp<-lm(deltapH~TempInSitu*Benthos*Season, data = AllCO2 %>% filter(Benthos != "Open Ocean"))
anova(mod_pH_Temp)
summary(mod_pH_Temp)

### 
AllCO2 %>%
  filter(NO3<6,
         PO4<1,
         NH4 <6)%>%
  ggplot(aes(x = TempInSitu, y = (NH4), color = Benthos))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Benthos, scale = "free")

AllCO2 %>%
  filter(NO3<6,
         PO4<1,
         NH4 <6)%>%
  ggplot(aes(x = TempInSitu, y = NO2))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Benthos, scale = "free")


AllCO2 %>%
  mutate(Ocean_not = ifelse(Benthos == "Open Ocean", "Ocean","Benthic")) %>%
  ggplot(aes(x = TempInSitu, y = deltaSIL, color = Benthos))+ ## eating diatoms at night??
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Benthos, scale = "free")


### run a pca of water chem

chem_only<-AllCO2 %>%
  drop_na(pH, TA, DIC, pCO2, CO2, CO3, HCO3) %>%
  #filter(NO3 < 3)%>%
 # drop_na( NO3, PO4, SIL, NO2, NH4, pH, TA)%>%
  select(pH, TA, DIC, pCO2, CO2, CO3, HCO3)

pca<-prcomp(chem_only,scale. = TRUE, center = TRUE )

# calculate percent explained by each PC
perc.explained<-round(100*pca$sdev/sum(pca$sdev),1)

# Extract the scores and loadings
PC_scores<-as_tibble(pca$x[,1:2])


PC_loadings<-as_tibble(pca$rotation) %>%
  bind_cols(labels = rownames(pca$rotation))


pca_all<-AllCO2 %>%
 # filter(NO3 < 3)%>%
  drop_na( pH, TA, DIC, pCO2, CO2, CO3, HCO3)%>%
  bind_cols(PC_scores)

scores_plot<-pca_all %>%
  mutate(ocean_not = ifelse(Benthos == "Open Ocean", "Ocean","Benthic"))%>%
  ggplot(aes(x = PC1, y = PC2, color =  ocean_not))+
  # coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  # scale_shape_manual(values = c(1, 22,15,16))+
  # scale_colour_manual(values = c("#D64550","#EA9E8D"))+
  # scale_fill_manual(values = c("#D64550","#EA9E8D"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(#fill = Tide, 
      label = paste(ocean_not), color =ocean_not), 
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap=0, tol = 0.05)+
  geom_point(size = 2) +
  labs(
    x = paste0("PC1 ","(",perc.explained[1],"%)"),
    #x = "",
    y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
       # axis.text.x = element_blank(),
      #  axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
        strip.text = element_blank())+
  facet_wrap(~Season, nrow=2)


## loadings plots
p2_loadings<-PC_loadings %>%
  ggplot(aes(x=PC1+0.1, y=PC2+0.1, label=labels))+
  geom_text(aes(x = PC1, y = PC2 ), show.legend = FALSE, size = 5) +
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1,yend=PC2),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
 # coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(color ="",
      # y = "",
    x = paste0("PC1 ","(",perc.explained[1],"%)"),
    y = paste0("PC2 ","(",perc.explained[2],"%)"))+
#  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.75, 0.75),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

scores_plot+p2_loadings

## PCA with nuts
chem_only<-AllCO2 %>%
  drop_na(NO3, NO2, NH4, PO4, SIL) %>%
  filter(NO3 < 3,
         PO4 < 3)%>%
  # drop_na( NO3, PO4, SIL, NO2, NH4, pH, TA)%>%
  select(NO3, NO2, NH4, PO4, SIL)

pca<-prcomp(chem_only,scale. = TRUE, center = TRUE )

# calculate percent explained by each PC
perc.explained<-round(100*pca$sdev/sum(pca$sdev),1)

# Extract the scores and loadings
PC_scores<-as_tibble(pca$x[,1:2])


PC_loadings<-as_tibble(pca$rotation) %>%
  bind_cols(labels = rownames(pca$rotation))


pca_all<-AllCO2 %>%
   filter(NO3 < 3,
          PO4 < 3)%>%
  drop_na( NO3, NO2, NH4, PO4, SIL)%>%
  bind_cols(PC_scores)

scores_plot<-pca_all %>%
  mutate(ocean_not = ifelse(Benthos == "Open Ocean", "Ocean","Benthic"))%>%
  ggplot(aes(x = PC1, y = PC2, color =  ocean_not))+
  # coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  # scale_shape_manual(values = c(1, 22,15,16))+
  # scale_colour_manual(values = c("#D64550","#EA9E8D"))+
  # scale_fill_manual(values = c("#D64550","#EA9E8D"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(#fill = Tide, 
      label = paste(ocean_not), color =ocean_not), 
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap=0, tol = 0.05)+
  geom_point(size = 2) +
  labs(
    x = paste0("PC1 ","(",perc.explained[1],"%)"),
    #x = "",
    y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
        # axis.text.x = element_blank(),
        #  axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
        strip.text = element_blank())+
  facet_wrap(~Season, nrow=2)


## loadings plots
p2_loadings<-PC_loadings %>%
  ggplot(aes(x=PC1+0.1, y=PC2+0.1, label=labels))+
  geom_text(aes(x = PC1, y = PC2 ), show.legend = FALSE, size = 5) +
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1,yend=PC2),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
  # coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(color ="",
       # y = "",
       x = paste0("PC1 ","(",perc.explained[1],"%)"),
       y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  #  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.75, 0.75),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

scores_plot+p2_loadings
