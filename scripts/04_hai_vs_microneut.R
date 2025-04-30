#Index
# 1) Load libraries and data 
# 2) Correlation of HAI and microneutralisation assays
# 3) Comparing "response" 
    #3a) Categorise FC <2 and <4 and high/low baseline titers
    #3b) Preparation for Heatmap plots in quarto file
    #3c) 

# 1) -------------------------------------------------------------------------------------------------------------------
#load libraries
library(dplyr)
library(tidyr)
library(haven)
library(here)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(car)
library(broom)
library(forcats)
library(stringr)

#load data
microneut_analysis_raw  <- read.xlsx(here::here("processed", "microneut_analysis_raw.xlsx"))
hai_analysis_raw  <- read.xlsx(here::here("processed", "hai_analysis_raw.xlsx"))
influenza_antibody_results <- read.xlsx(here::here("processed", "influenza_antibody_results.xlsx"))
basefile_withPID <- read.xlsx(here::here("processed", "FluVac_basefile_withPID.xlsx"))

# 2) -------------------------------------------------------------------------------------------------------------------
#### Comparing HAI and Microneutralisation assay
#preparation
microneut_analysis_raw_p <- microneut_analysis_raw %>% #CAVE: USED IN OTHER SCRIPTS - e.g. 02_hai_micro...
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50, PID) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50")  

#preparation
hai_analysis_raw_p <- hai_analysis_raw %>% #CAVE: USED IN OTHER SCRIPTS - e.g. 02_hai_micro...
  select(pat_group, Sampling_number, H1N1, H3N2, BVic, BYam, PID) %>% 
  pivot_longer(cols = 3:6, names_to = "strain", values_to = "result")  


#harmonize dataset for analysis - same variables, row_append hai
influenza_antibody_harmonised <- microneut_analysis_raw_p %>% 
  mutate(strain = case_when( strain == "H1" ~ "H1N1",
                             strain == "H3" ~ "H3N2",
                             strain == "Vic" ~ "BVic")) %>% 
  mutate(type = "ic50_microneutr") %>% 
  select(pat_group, Sampling_number, strain, result = ic50, type, PID) %>%  
  rows_append(hai_analysis_raw_p %>% mutate(type = "hai")) %>% print()

tf <- influenza_antibody_harmonised %>% 
  filter(type == "ic50_microneutr") %>% 
  select(pat_group, Sampling_number, strain, ic50res = result, PID, type)

tf2 <- influenza_antibody_harmonised %>% 
  filter(type == "hai") %>% 
  select(pat_group, Sampling_number, strain, haires = result, PID, type)

tf3 <- tf2 %>%  
  left_join(tf, by = join_by(PID, pat_group, Sampling_number, strain))  %>% 
  filter(!is.na(ic50res)) %>% 
  select(pat_group, Sampling_number, strain, haires, PID, ic50res)

#using corr by strain
strain_corrs <- tf3 %>%
  group_by(strain) %>%
  summarise(
    cor_result = list(cor.test(haires, ic50res, method = "spearman")),
    .groups = "drop"
  )

#print correlation table
strain_corrs_tidy <- strain_corrs %>%
  mutate(tidy_result = lapply(cor_result, broom::tidy)) %>%
  unnest(tidy_result) %>%
  select(strain, estimate, statistic, p.value, method)

#rank values by strain and type(Hai vs. Microneut) for visual linear correlation - ranking only within the same strain
influenza_antibody_ranked <- influenza_antibody_harmonised %>% 
  group_by(type, strain) %>% 
  mutate(rank_result = rank(result)) %>%
  ungroup() %>% 
  print()

df <- influenza_antibody_ranked %>% 
  filter(type == "ic50_microneutr") %>% 
  select(PID, pat_group, Sampling_number, strain, ic50_rank = rank_result)

df2 <- influenza_antibody_ranked %>% 
  filter(type == "hai") %>% 
  select(PID, pat_group, Sampling_number, strain, hai_rank = rank_result)

df3 <- df2 %>%  
  left_join(df, by = join_by(PID, pat_group, Sampling_number, strain))

hai_vs_ic50_ranksumplot <- ggplot(data = df3, aes(x = hai_rank, y = ic50_rank)) +
  geom_jitter(aes(color = strain), size = 2, alpha = 0.5) +  # Correct jitter + color mapping
  geom_smooth(aes(color = strain), method = "lm", se = TRUE) +   # Optional: consistent regression line
  labs(title = "Rangsummenkorrelation: HAI vs IC50",
       x = "HAI Rank",
       y = "IC50 Rank",
       color = "Strain") +
  theme_minimal()

# 3) -------------------------------------------------------------------------------------------------------------------
###Comparing "responders" 
# 3) Categorise FC defined by FC <2, 2-4 and >4 and high/low baseline titers
influenza_folds_categ <-  influenza_antibody_results %>%
  mutate(across(
    c(hai_H1_fold, hai_H3_fold, hai_BVic_fold, hai_BYam_fold,
      FluV_H1_fold, FluV_H3_fold, FluV_Vic_fold),
    ~ case_when(
      .x < 2 ~ 1,
      .x >= 2 & .x < 4 ~ 3,
      .x >= 4 ~ 4,
      TRUE ~ 99
    )
  )) %>%
  filter(Sampling_number == 2) %>%
  select(Pat_ID, pat_group, Sampling_number, SamplingDt,
         hai_H1_fold, FluV_H1_fold, hai_H3_fold, FluV_H3_fold,
         hai_BVic_fold, FluV_Vic_fold, hai_BYam_fold) 

# 3b) 
#Reshape Data to Long Format for Plotting
influenza_folds_categ_long <- influenza_folds_categ %>%
  pivot_longer(
    cols = -c(Pat_ID, Sampling_number, SamplingDt, pat_group),
    names_to = "Antibody_Type",
    values_to = "Fold_Category"
  ) %>% 
  mutate(Pat_ID = as.character(Pat_ID)) %>% 
  group_by(Pat_ID) %>% 
  mutate(totalfold = sum(Fold_Category)) %>% 
  ungroup() %>% 
  mutate(
    Fold_Category = factor(Fold_Category),
    Pat_ID = fct_reorder(Pat_ID, totalfold, .desc = TRUE)  # Reorder by totalfold (high to low)
) %>% 
  mutate(
    Antibody_Type = factor(
      Antibody_Type,
      levels = c("hai_H1_fold", "FluV_H1_fold",
                 "hai_H3_fold", "FluV_H3_fold",
                 "hai_BVic_fold", "FluV_Vic_fold",
                 "hai_BYam_fold")
    )
  )


influenza_folds_custom <- influenza_folds_categ_long %>%
  mutate(x = case_when(
    Antibody_Type == "hai_H1_fold"   ~ 1,
    Antibody_Type == "FluV_H1_fold"  ~ 1.5,
    Antibody_Type == "hai_H3_fold"   ~ 3,
    Antibody_Type == "FluV_H3_fold"  ~ 3.5,
    Antibody_Type == "hai_BVic_fold" ~ 5,
    Antibody_Type == "FluV_Vic_fold" ~ 5.5,
    Antibody_Type == "hai_BYam_fold" ~ 7,
    TRUE ~ NA_real_
  ))

x_label_map <- c(
  "1"   = "H1 (HAI)",
  "1.5" = "H1 (Micr)",
  "3"   = "H3 (HAI)",
  "3.5" = "H3 (Micr)",
  "5"   = "B/Vic (HAI)",
  "5.5" = "B/Vic (Micr)",
  "7"   = "B/Yam (HAI)"
)

# 3c) 
influenza_folds_categ_long_str <- influenza_folds_categ_long %>%
  mutate(
    strain = case_when(
      Antibody_Type %in% c("hai_H1_fold", "FluV_H1_fold") ~ "H1",
      Antibody_Type %in% c("hai_H3_fold", "FluV_H3_fold") ~ "H3",
      Antibody_Type %in% c("hai_BVic_fold", "FluV_Vic_fold") ~ "B/Victoria",
      Antibody_Type == "hai_BYam_fold" ~ "B/Yamagata",
      TRUE ~ "Unknown"
    )
  ) %>% 
  mutate(assay = case_when(
    Antibody_Type %in% c("hai_H1_fold", "hai_H3_fold", "hai_BVic_fold", "hai_BYam_fold" ) ~ "hai",
    Antibody_Type %in% c("FluV_H1_fold", "FluV_H3_fold", "FluV_Vic_fold" ) ~ "microneut",
    TRUE ~ "Unknown"
  )) 

#calculate assaydif-variable by microneut - hai (giving 0s for no difference, 1/-1 or 2/-2 for one category difference, and 3/-3 for two category difference) 
influenza_folds_dif  <- influenza_folds_categ_long_str %>% 
  select(Pat_ID, Fold_Category, assay, strain, Fold_Category) %>% 
  mutate(Fold_Category= as.character(Fold_Category)) %>% 
  mutate(Fold_Category= as.numeric(Fold_Category)) %>% 
  pivot_wider(names_from = assay, values_from = Fold_Category) %>% 
  mutate(microneut = as.numeric(microneut)) %>% 
  mutate(hai = as.numeric(hai)) %>% 
  mutate(assaydif = microneut - hai) %>% 
  select(Pat_ID, strain, assaydif) %>% 
  filter(!is.na(assaydif)) %>% 
  print()


influenza_folds_customdif <-  influenza_folds_custom %>% 
  mutate(strain = case_when(Antibody_Type %in% c("hai_H1_fold", "FluV_H1_fold") ~ "H1",
      Antibody_Type %in% c("hai_H3_fold", "FluV_H3_fold") ~ "H3",
      Antibody_Type %in% c("hai_BVic_fold", "FluV_Vic_fold") ~ "B/Victoria",
      Antibody_Type == "hai_BYam_fold" ~ "B/Yamagata",
      TRUE ~ "Unknown"
    ))%>% 
  left_join(influenza_folds_dif)


influenza_overall_hai<- influenza_folds_categ_long_str %>% 
  select(Pat_ID, pat_group, Fold_Category, strain, assay) %>% 
  filter(assay == "hai") %>% 
  group_by(Pat_ID) %>% 
  count(Fold_Category) %>% 
  ungroup(Pat_ID) %>% 
  mutate(fold_category = case_when(Fold_Category == 4 ~ "high",
                                   Fold_Category == 3 ~ "moderate",
                                   Fold_Category == 1 ~ "low",
                                   TRUE ~ "Undefined")) %>% 
  select(Pat_ID, n, fold_category) %>% 
  left_join(influenza_folds_categ_long_str %>% select(Pat_ID, pat_group) %>% distinct(Pat_ID, pat_group)) %>% 
  pivot_wider(names_from = fold_category, values_from = n) %>% 
  mutate(overall_response_hai = case_when(high == 4 ~ "all high",
                                      moderate == 4 ~ "all moderate",
                                      low == 4 ~ "all low",
                                      !is.na(high) & !is.na(moderate) & is.na(low) ~ "mixed high/moderate",
                                      is.na(high) & !is.na(moderate) & !is.na(low) ~ "mixed moderate/low",
                                      !is.na(high) & is.na(moderate) & !is.na(low) ~ "mixed high/low",
                                      !is.na(high) & !is.na(moderate) & !is.na(low) ~ "mixed all",
                                      TRUE ~ "Undefined")) %>% print()

influenza_overall_micr<- influenza_folds_categ_long_str %>% 
  select(Pat_ID, pat_group, Fold_Category, strain, assay) %>% 
  filter(assay == "microneut") %>% 
  group_by(Pat_ID) %>% 
  count(Fold_Category) %>% 
  ungroup(Pat_ID) %>% 
  mutate(fold_category = case_when(Fold_Category == 4 ~ "high",
                                   Fold_Category == 3 ~ "moderate",
                                   Fold_Category == 1 ~ "low",
                                   TRUE ~ "Undefined")) %>% 
  select(Pat_ID, n, fold_category) %>% 
  left_join(influenza_folds_categ_long_str %>% select(Pat_ID, pat_group) %>% distinct(Pat_ID, pat_group)) %>% 
  pivot_wider(names_from = fold_category, values_from = n) %>% 
  mutate(overall_response_micr = case_when(high == 3 ~ "all high",
                                      moderate == 3 ~ "all moderate",
                                      low == 3 ~ "all low",
                                      !is.na(high) & !is.na(moderate) & is.na(low) ~ "mixed high/moderate",
                                      is.na(high) & !is.na(moderate) & !is.na(low) ~ "mixed moderate/low",
                                      !is.na(high) & is.na(moderate) & !is.na(low) ~ "mixed high/low",
                                      !is.na(high) & !is.na(moderate) & !is.na(low) ~ "mixed all",
                                      TRUE ~ "Undefined")) %>% print()

influenza_overall <-  influenza_overall_micr %>% 
  select(Pat_ID, pat_group, overall_response_micr) %>% 
  left_join(influenza_overall_hai) %>% 
  select(Pat_ID, pat_group,overall_response_hai, overall_response_micr)

influenza_overall %>%  count(overall_response_micr) %>% 
  print()

#response for different strains
influenza_folds_customdif_wide <- influenza_folds_customdif %>% 
  select(Pat_ID, pat_group, Antibody_Type, Fold_Category) %>% 
  mutate(Fold_Category = as.character(Fold_Category)) %>%  
  mutate(Fold_Category = as.numeric(Fold_Category)) %>%  
  pivot_wider(names_from = Antibody_Type, values_from = Fold_Category)  %>% print() 
 

count_h1_hai <- influenza_folds_customdif_wide %>% 
  count(hai_H1_fold) %>% 
  select(fold_category = hai_H1_fold, H1_HAI = n)

count_h1_mic <- influenza_folds_customdif_wide %>% 
  count(FluV_H1_fold) %>% 
  select(fold_category = FluV_H1_fold, H1_mic = n)

count_h3_hai <- influenza_folds_customdif_wide %>% 
  count(hai_H3_fold) %>% 
  select(fold_category = hai_H3_fold, H3_HAI = n)

count_h3_mic <- influenza_folds_customdif_wide %>% 
  count(FluV_H3_fold) %>% 
  select(fold_category = FluV_H3_fold, H3_mic = n)

count_bvic_hai <- influenza_folds_customdif_wide %>% 
  count(hai_BVic_fold) %>% 
  select(fold_category = hai_BVic_fold, BVic_HAI = n)

count_bvic_mic <- influenza_folds_customdif_wide %>% 
  count(FluV_Vic_fold) %>% 
  select(fold_category = FluV_Vic_fold, BVic_mic = n)

count_byam_hai <- influenza_folds_customdif_wide %>% 
  count(hai_BYam_fold) %>% 
  select(fold_category = hai_BYam_fold, BYam_HAI = n)

count_table_response <-  count_h1_hai %>% 
  left_join(count_h1_mic) %>% 
  left_join(count_h3_hai) %>% 
  left_join(count_h3_mic) %>% 
  left_join(count_bvic_hai) %>% 
  left_join(count_bvic_mic) %>% 
  left_join(count_byam_hai) %>% 
  mutate(fold_category = as.character(case_when(fold_category == 1 ~ "low response",
                                                fold_category == 3 ~ "moderate response",
                                                fold_category == 4 ~ "high response",
                                                TRUE ~ "undefinded")))


count_table_response2 <-count_table_response %>%
  pivot_longer(cols = 2:8, names_to = "strain", values_to = "count") %>%
  mutate(assay = case_when(
    str_detect(strain, "HAI") ~ "hai",
    str_detect(strain, "mic") ~ "mic",
    TRUE ~ "undefined"
  )) %>% 
  mutate(strain = case_when(
    str_detect(strain, "H1") ~ "H1",
    str_detect(strain, "H3") ~ "H3",
    str_detect(strain, "BVic") ~ "B/Victoria",
    str_detect(strain, "BYam") ~ "B/Yamagata",
    TRUE ~ "undefined"
  )) %>% 
  mutate(fold_category = factor(fold_category, levels = c("low response", "moderate response", "high response")))

count_table_response2_percent <- count_table_response2 %>%
  group_by(strain, assay) %>%
  mutate(percentage = count / sum(count) * 100)
