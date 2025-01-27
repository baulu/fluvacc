#load libraries
library(dplyr)
library(tidyr)


#Import data (PROVISORISCH)
send_items <- read.csv2("/Users/lu/Desktop/Home/DS5752_VersandUK_Serum_2025-01.csv")
pidfile <- read.csv2("//Users/lu/Desktop/Home/FluVacc_PIDs.csv")

#Create PID-file with Serum samples for later HAI/Microneut.-Results
joint_file <- send_items %>% 
  mutate(PID = Subject) %>% 
  select(PID, Position, CryotubeID, SamplingDt) %>% 
  left_join(pidfile, by = join_by(PID)) 

joint_file2 <- pidfile %>% 
  mutate(Subject = PID) %>% 
  left_join(send_items, by = join_by(Subject)) %>% 
  mutate(PID = Subject) %>% 
  select(PID, Pat.ID, Position, CryotubeID, SamplingDt, Day.0, Day.7, Day.28., Drop.Out.) %>%
  # filter(Drop.Out. == "YES") #10 Dropouts, with only one providing a sample at baseline
  filter(!is.na(SamplingDt)) # Keep only IDs with an available sample

# Identify PIDs with only one sample (n=3)
joint_file2 %>% 
  group_by(PID) %>% 
  arrange(PID, SamplingDt) %>% 
  mutate(Sampling_number = row_number()) %>%
  select(PID, Pat.ID, CryotubeID, SamplingDt, Day.0, Day.7, Day.28., Drop.Out.,Sampling_number) %>%
  arrange(PID, Sampling_number) %>% 
  slice_tail() %>% 
  filter(Sampling_number == 1) %>%  # FluV_CO_041, FluV_MS_024, FluV_HIV_027(DropOUT)
  ungroup() 

  