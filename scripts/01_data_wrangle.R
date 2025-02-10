#load libraries
library(dplyr)
library(tidyr)
library(haven)
library(here)
library(lubridate)
library (openxlsx)

#Import data (PROVISORISCH)
send_items <- read.csv2("data/provisorisch/DS5752_VersandUK_Serum_2025-01.csv") %>% 
  mutate(SamplingDt = dmy(SamplingDt))
pidfile <- read.csv2("data/provisorisch/FluVacc_PIDs.csv")
basefile <- read.csv2("data/provisorisch/IMVASU2_DATA_2025-02-04_0848.csv")

#Create PID-file with Serum samples for later HAI/Microneut.-Results
joint_file <- send_items %>% #148 Subjects
  mutate(PID = Subject) %>% 
  select(PID, Studyname, Requestid, RackID, Position, CryotubeID, Material, Storage, SamplingDt) %>% 
  left_join(pidfile, by = join_by(PID)) %>% 
  group_by(PID) %>% 
  arrange(PID, SamplingDt) %>% 
  mutate(Sampling_number = row_number()) %>%
  ungroup()

joint_file2 <- pidfile %>% #156 inclusions (with 9 drop-outs and 1 drop-out with sample at baseline = 147 of send_items list)
  mutate(Subject = PID) %>% 
  left_join(send_items, by = join_by(Subject)) %>% 
  mutate(PID = Subject) %>% 
  select(PID, Pat.ID, Position, CryotubeID, SamplingDt, Day.0, Day.7, Day.28., Drop.Out.) %>%
  # filter(Drop.Out. == "YES") #10 Dropouts, with only one providing a sample at baseline
  filter(!is.na(SamplingDt)) # Keep only IDs with an available sample

# Identify PIDs with only one serum-sample (n=3)
joint_file2 %>% 
  group_by(PID) %>% 
  arrange(PID, SamplingDt) %>% 
  mutate(Sampling_number = row_number()) %>%
  select(PID, Pat.ID, CryotubeID, SamplingDt, Day.0, Day.7, Day.28., Drop.Out.,Sampling_number) %>%
  arrange(PID, Sampling_number) %>% 
  slice_tail() %>% 
  filter(Sampling_number == 1) %>%  # FluV_CO_041, FluV_MS_024, FluV_HIV_027(DropOUT)
  ungroup() 

##### Create list for Francis Crick with age group etc.
# patient characteristics file
basefile_select <- basefile %>% 
  select(study_id, redcap_event_name, gen_sex, gen_age, sv1_studygroup) %>% 
  group_by(study_id) %>% 
  arrange(redcap_event_name) %>% 
  slice_head() %>% 
  ungroup() %>% 
  mutate(Pat.ID = study_id) %>% 
  mutate(study_group = case_when(sv1_studygroup == 1 ~ "control",
                                 sv1_studygroup == 2 ~ "onco",
                                 sv1_studygroup == 3 ~ "ms",
                                 sv1_studygroup == 4 ~ "hiv",
                                 sv1_studygroup == 5 ~ "reheuma")) %>% 
  mutate(age_group = case_when(gen_age <= 64 ~ "<65",
                               gen_age >= 64 ~ "≥65"))
  
#append patient characteristics to sent-file for Crick institute
FluVac_HAI_Samples_ext_02_25 <-  joint_file %>% 
  left_join(basefile_select, join_by(Pat.ID)) %>% 
  group_by(PID) %>%
  arrange(PID, SamplingDt) %>% 
  mutate(date_diff = as.numeric(SamplingDt - lag(SamplingDt))) %>% #calculate days between sample 1 and 3
  mutate(PID = cur_group_id()) %>%
  ungroup() %>% 
  select(Studyname, Requestid, RackID, Position, CryotubeID, Material, Storage, Sampling_number, study_group, age_group, date_diff) %>% 
  print()


write.xlsx(FluVac_HAI_Samples_ext_02_25, file="processed/FluVac_HAI_Samples_ext_02_25.xlsx", overwrite = TRUE, asTable = TRUE)

################ PBMC analysis
send_items_pbmc <- read.csv2("/Users/lu/Desktop/Home/PBMCs_forSend.csv") 

joint_file_pbmc <- pidfile %>% #148 -> correct
  mutate(Subject = PID) %>% 
  left_join(send_items_pbmc, by = join_by(Subject)) %>% 
  mutate(PID = Subject) %>% 
  select(PID, Pat.ID, Timepoint.Check, LabScan, StorageDt, CryoID, SampleID, Status,Drop.Out.) %>% 
  # filter(Drop.Out. == "YES") #10 Dropouts
  filter(!is.na(LabScan)) # Keep only IDs with an available sample

# Identify PIDs with only one PBMC-Sample (n=?)
joint_file_pbmc %>% 
  group_by(PID) %>% 
  arrange(PID, Timepoint.Check) %>% 
  slice_tail() %>%
  filter(Timepoint.Check == 1) %>%  # One with only baseline data FluV_HIV_027
  ungroup() 
  