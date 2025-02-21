#load libraries
library(dplyr)
library(tidyr)
library(haven)
library(here)
library(lubridate)
library (openxlsx)
library(ggplot2)
library(stringr)
library(readr)

#Import data (PROVISORISCH)
send_items <- read_csv2("data/provisorisch/DS5752_VersandUK_Serum_2025-01.csv") %>% 
  mutate(SamplingDt = dmy(SamplingDt))
pidfile <- read_csv2("data/provisorisch/FluVacc_PIDs.csv")
basefile <- read_csv2("data/provisorisch/IMVASU2_DATA_2025-02-04_0848.csv")
biobank_report <- read.csv2("data/provisorisch/Report Study 5752 - FluVacc2025-02-03.csv")# %>% 
  filter(Cryo.Barcode != "FD31518743" & Cryo.Barcode != "FD31518744" & Cryo.Barcode !="FD31518757" & Cryo.Barcode !="FD31518758")# Exlude for empty tubes for FluV_CO_32
microneut_results_raw <- read_csv("data/final_results_microneut_2025-02-20_FluVacc.csv") %>% 
  mutate(CryotubeID = sample_barcode)



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
  mutate(syntheticID = cur_group_id()) %>% 
  ungroup() %>% 
  select(Studyname, Requestid, RackID, Position, CryotubeID, Material, Storage, syntheticID, Sampling_number, study_group, age_group, date_diff) %>% 
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
  
# Analysis of Samples - how many?

biobank_report_sum <- biobank_report %>% 
  select(Primary.PID, Primary.Proband.ID, Lab.Scan.Dt, Lab.Order.ID, Cryo.Sampletype, Cryo.Barcode, Cryo.Volume.ul, Cryo.Status, Cryo.Storagestatus, Cryo.Concentration, Cryo.Concentrationunits) %>% 
  mutate(Lab.Scan.Dt = parse_date_time(Lab.Scan.Dt,  "mdYHMS")) %>% 
  mutate(Lab.Scan.Dt = as_date(Lab.Scan.Dt)) %>% 
  filter(Primary.Proband.ID == "") %>% 
  filter(Cryo.Storagestatus == "In Circulation") %>% 
  group_by(Primary.PID) %>% 
  arrange(Primary.PID, Lab.Scan.Dt) %>% 
  mutate(VisitNR = dense_rank(Lab.Scan.Dt)) %>%  #Assign a visitNumber
  ungroup() %>% 
  group_by(Primary.PID, VisitNR, Cryo.Sampletype) %>% 
  arrange(Primary.PID, VisitNR, Cryo.Sampletype) %>% 
  mutate(aliquotsPerVis = row_number()) %>% 
  slice_tail() %>% 
  ungroup() %>% 
  select(Primary.PID, Lab.Scan.Dt, Cryo.Sampletype, VisitNR, aliquotsPerVis) %>% print() 


ggplot(biobank_report_sum, aes(x = VisitNR)) +
  geom_histogram() +
  labs(title = "Histogram of VisitNR", x = "Visit Number", y = "Count") +
  theme_minimal()

#PBMC Nr. of Aliquots per Visit
PBMC_aliquotNR <- biobank_report_sum %>% 
  filter(Cryo.Sampletype == "PBMC") %>% 
  filter(VisitNR <= 3) %>% 
  count(VisitNR, aliquotsPerVis, Cryo.Sampletype) %>% 
  group_by(VisitNR) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  ungroup() %>% 
  mutate(combined = str_c(n, " (", round(percent, 1), "%)")) %>% 
  select(VisitNR, aliquotsPerVis, combined) %>% 
  pivot_wider(names_from = aliquotsPerVis, values_from = combined) %>% 
  select(Visit_NR = VisitNR, One_Aliquote = "1", Two_Aliquotes = "2", Three_Aliquotes = "3", Four_Aliquotes = "4") %>% 
  mutate(Visit_NR = as.character(Visit_NR))

#Nr. of Patients with only 1 or two PBMC-aliquots at V1 and/or 3
PBMC_aliquotNR_min <- biobank_report_sum %>% 
  filter(Cryo.Sampletype == "PBMC") %>% 
  filter(VisitNR == 1 | VisitNR == 3) %>% 
  group_by(Primary.PID) %>% 
  arrange(Primary.PID, aliquotsPerVis) %>% 
  slice_head() %>% 
  ungroup() %>% 
  count(aliquotsPerVis) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  mutate(combined = str_c(n, " (", round(percent, 1), "%)")) %>% 
  select(aliquotsPerVis, combined) %>% 
  pivot_wider(names_from = aliquotsPerVis, values_from = combined) %>% 
  mutate(Visit_NR = "minimalAliquot (Only Visit 1&3)") %>% 
  select(Visit_NR, One_Aliquote = `1`, Two_Aliquotes = `2`, Three_Aliquotes = `3`, Four_Aliquotes = `4`) %>% 
  mutate(across(everything(), as.character))  # Convert all variables to character

PBMC_aliquotNR_comb <- PBMC_aliquotNR %>% 
  bind_rows(PBMC_aliquotNR_min)

write.xlsx(PBMC_aliquotNR_comb, file="processed/PBMC_aliquotNR.xlsx", overwrite = TRUE, asTable = TRUE)


#Serum Nr. of Aliquots per Visit
Serum_aliquotNR <- biobank_report_sum %>% 
  filter(Cryo.Sampletype == "Serum") %>% 
  filter(VisitNR <= 3) %>% 
  count(VisitNR, aliquotsPerVis, Cryo.Sampletype) %>% 
  group_by(VisitNR) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  ungroup() %>% 
  mutate(combined = str_c(n, " (", round(percent, 1), "%)")) %>% 
  select(VisitNR, aliquotsPerVis, combined) %>% 
  pivot_wider(names_from = aliquotsPerVis, values_from = combined) %>% 
  select(Visit_NR = VisitNR, One_Aliquote = "1", Two_Aliquotes = "2", Three_Aliquotes = "3")

write.xlsx(Serum_aliquotNR, file="processed/Serum_aliquotNR.xlsx", overwrite = TRUE, asTable = TRUE)

#### Looking at Microneutralisation data------------------------------------------------
## combining files
microneut_combined <- joint_file %>%
  left_join(microneut_results_raw, join_by(CryotubeID))

## raw first look at data
microneut_analysis_raw <- microneut_combined %>% 
  select(PID, CryotubeID, Pat_ID =`Pat-ID`, SamplingDt, Sampling_number, FluA_H1_ic50 = `FluA/H1_ic50`, FluA_H3_ic50 = `FluA/H3_ic50`,FluA_Vic_ic50 = `FluB/Vic_ic50`) %>% 
  arrange(Pat_ID, PID, Sampling_number) %>% 
  group_by(PID) %>% 
  mutate(FluV_H1_fold = FluA_H1_ic50 / lag(FluA_H1_ic50)) %>%
  mutate(FluV_H3_fold = FluA_H3_ic50 / lag(FluA_H3_ic50)) %>%
  mutate(FluV_Vic_fold = FluA_Vic_ic50 / lag(FluA_Vic_ic50)) %>%
  ungroup() %>% 
  mutate(pat_group = case_when(
    str_detect(Pat_ID, "CO")  ~ "Control",
    str_detect(Pat_ID, "HIV") ~ "HIV",
    str_detect(Pat_ID, "MS") ~ "MS", 
    str_detect(Pat_ID, "ONK") ~ "ONK",
    str_detect(Pat_ID, "RH") ~ "RH", # Example: add more conditions
    TRUE ~ "Pat_ID" # Default: Assigns NA if no match 
  )) %>% 
  mutate(H1_ov4 = FluV_H1_fold > 4) %>%  
  mutate(H3_ov4 = FluV_H3_fold > 4) %>%  
  mutate(Vic_ov4 = FluV_Vic_fold > 4) %>%  
  print()

micneut_vic <- microneut_analysis_raw %>% 
  group_by(pat_group) %>% 
  count(Vic_ov4) %>% filter(!is.na(Vic_ov4))%>% 
  select(pat_group, Over4 = Vic_ov4, Vic = n) %>% 
  print()

micneut_H1 <- microneut_analysis_raw %>% 
  group_by(pat_group) %>% 
  count(H1_ov4) %>% filter(!is.na(H1_ov4))%>% 
  select(pat_group, Over4 = H1_ov4, H1 = n) %>% 
  print()

micneut_H3 <- microneut_analysis_raw %>% 
  group_by(pat_group) %>% 
  count(H3_ov4) %>% filter(!is.na(H3_ov4))%>% 
  select(pat_group, Over4 = H3_ov4, H3 = n) %>% 
  print()

micneut_foldoverview <- micneut_vic %>% 
  left_join(micneut_H1) %>% 
  left_join(micneut_H3)

write.xlsx(microneut_analysis_raw, file="processed/microneut_analysis_raw.xlsx", overwrite = TRUE, asTable = TRUE)
write.xlsx(micneut_foldoverview, file="processed/micneut_foldoverview.xlsx", overwrite = TRUE, asTable = TRUE)

