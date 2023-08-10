#### Temperature processing Data#####
#### Created by Nyssa Silbiger ####
#### created on 8/10/2023 #####


### Load Libraries####
library(tidyverse)
library(here)
library(janitor)



## read in the data ####

# pull out the file path
path<-here("Data","Temperature","Summer_2023")

# extract the files
files <- dir(path = path,pattern = ".csv", full.names = TRUE)

# create a function to read in and clean the files
read_fun<-function(name){
  read_csv({{name}}) %>%
    clean_names()%>%
    rename(Temperature = ch_1_temperature_c,
           Lux = ch_2_light_lux,
           DateTime = date_time_pdt) %>%
    select(DateTime,Temperature,Lux) %>%
    mutate(DateTime = mdy_hms(DateTime))
}

# read in all the files into one dataframe
data<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_fun,.id = "filename") %>%
  mutate(filename = str_split_i(filename, "/",10),# extract the actual filename
         BenthicID = str_split_i(filename, "_",1)) %>%# extract the plotID
  mutate(Benthos = case_when(grepl("B", BenthicID) ~ "Barnacles",
                           grepl("M", BenthicID) ~"Mussels",
                           grepl("P", BenthicID) ~"Surfgrass",
                           grepl("R", BenthicID) ~"Rockweed"))
         




data %>%
  drop_na(Benthos)%>%
  ggplot(aes(x = DateTime, y = Temperature, group = BenthicID, color = Benthos))+
  geom_line()

data %>%
  drop_na(Benthos)%>%
  ggplot(aes(x = Benthos, y = Temperature, color = Benthos))+
  geom_boxplot()

