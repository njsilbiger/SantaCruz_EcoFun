#### Temperature processing Data#####
#### Created by Nyssa Silbiger ####
#### created on 8/10/2023 #####


### Load Libraries####
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(calecopal)
library(patchwork)


## read in the data ####

# pull out the file path
path1<-here("Data","Temperature","Summer_2023")
path2<-here("Data","Temperature","Fall_2023")
path3<-here("Data","Temperature","Winter_2024")

# extract the files
files1 <- dir(path = path1,pattern = ".csv", full.names = TRUE)
files2 <- dir(path = path2,pattern = ".csv", full.names = TRUE)
files3 <- dir(path = path3,pattern = ".csv", full.names = TRUE)

files<-c(files1, files3)
# create a function to read in and clean the files
read_fun2<-function(name){
  read_csv({{name}}) %>%
    clean_names()%>%
    rename(Temperature = ch_1_temperature_c,
           Lux = ch_2_light_lux,
           DateTime = date_time_pdt) %>%
    select(DateTime,Temperature,Lux) %>%
    mutate(DateTime = case_when(
      grepl("^0", DateTime) ~ mdy_hm(DateTime),
      .default = mdy_hms(DateTime))) %>% # some files are hm and others are hms
    slice(100:n()) # delete the first 100 points
}

read_fun<-function(name){
  read_csv({{name}}, skip = 100, col_names = c("X1","DateTime", "Temperature","Lux" )) %>%
    select(DateTime,Temperature,Lux) %>%
    mutate(DateTime = as.character(DateTime))%>%
    mutate(DateTime = case_when(
      grepl("^0", DateTime) ~ mdy_hms(DateTime),
      .default = mdy_hm(DateTime))) # some files are hm and others are hms
  } 

data2<-files2 %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_fun2,.id = "filename") 
  

# read in all the files into one dataframe
data<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_fun,.id = "filename") %>%
  bind_rows(data2)%>%
  mutate(filename = str_split_i(filename, "/",10),# extract the actual filename
         BenthicID = str_split_i(filename, "_",1),
         Season = str_split_i(filename, "_",2)) %>%# extract the plotID
  mutate(Benthos = case_when(grepl("B", BenthicID) ~ "Barnacles",
                           grepl("M", BenthicID) ~"Mussels",
                           grepl("P", BenthicID) ~"Surfgrass",
                           grepl("R", BenthicID) ~"Rockweed",
                           grepl("O", BenthicID) ~"Open Ocean"))
         

p1<-data %>%
  filter(Benthos != "Open Ocean") %>%
  ggplot(aes(x = DateTime, y = Temperature, group = BenthicID, color = Benthos))+
  geom_line()+
  scale_color_manual(values = cal_palette("tidepool"))+
  theme_bw()+
  facet_wrap(Season~Benthos, scales = "free_x", ncol = 4)+
  labs(x = "",
       y = expression("Temperature ("*degree*"C)"))+
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16))

p2<-data %>%
  filter(Benthos != "Open Ocean",
         BenthicID != "P4" 
         ) %>%
  group_by(Benthos, BenthicID, Season)%>%
  summarise(Temp_max = max(Temperature, na.rm = TRUE))%>%
  ungroup()%>%
  group_by(Benthos, Season) %>%
  summarise(Max_mean = mean(Temp_max, na.rm = TRUE),
            Max_se =  sd(Temp_max, na.rm = TRUE)/sqrt(n())) %>%
ggplot(aes(x = Benthos, y = Max_mean, color = Benthos, shape = Season))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = Max_mean - Max_se, ymax = Max_mean+Max_se, x = Benthos), width =.1)+
  scale_color_manual(values = cal_palette("tidepool"))+
  labs(y=expression("Average max temperature ("*degree*"C)"),
       x = "")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
  
p1/p2
ggsave(here("Output","Temperaturefig.png"), width = 12, height = 10)
