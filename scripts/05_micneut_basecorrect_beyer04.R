#Index
# 1) Load libraries and data 
# 2) Correction function of outcome titers to baselin titers (according to beyer 2004)
# 3) Create dataset for appliing correction function
# 4) Apply function to dataset

# 1) -------------------------------------------------------------------------------------------------------------------
#load libraries
library(dplyr)
library(tidyr)
library(haven)
library(here)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(stringr)
library(readr)
library(scales)
library(patchwork)

#load data
microneut_analysis_raw  <- read.xlsx(here::here("processed", "microneut_analysis_raw.xlsx"))
hai_analysis_raw  <- read.xlsx(here::here("processed", "hai_analysis_raw.xlsx"))
influenza_antibody_results <- read.xlsx(here::here("processed", "influenza_antibody_results.xlsx"))
basefile_withPID <- read.xlsx(here::here("processed", "FluVac_basefile_withPID.xlsx"))

# 2) -------------------------------------------------------------------------------------------------------------------
# Correction function of outcome titers to baselin titers (according to beyer 2004)
correct_T_post <- function(outcome_titers, baseline_titers, baseline_constant, slope) {
  log2(outcome_titers) - slope * (log2(baseline_titers) - baseline_constant)
}

#3) -------------------------------------------------------------------------------------------------------------------
# Create dataset for appliing correction function
basetiters <-  microneut_analysis_raw %>%
  filter(Sampling_number == 1) %>% 
  select(PID, H3_baseline = FluA_H3_ic50, H1_baseline = FluA_H1_ic50, B_baseline = FluA_Vic_ic50) %>%
  as.data.frame()

outcometiters <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, H3_outcome = FluA_H3_ic50, H1_outcome = FluA_H1_ic50, B_outcome = FluA_Vic_ic50) %>%
  as.data.frame()
  
titerset_combined <- basetiters %>% 
  left_join(outcometiters) %>% 
  na.omit() 

titerset_b <- titerset_combined %>% 
  select(PID, baseline_titers = H3_baseline, outcome_titers = H3_outcome)

titerset_h1n1 <- titerset_combined %>% 
  select(PID, baseline_titers = H1_baseline, outcome_titers = H1_outcome)

titerset_h3n2 <- titerset_combined %>% 
  select(PID, baseline_titers = H3_baseline, outcome_titers = H3_outcome)

#slopes: 
# b: 0.3063  
b_model_titers <- lm(log2(outcome_titers) ~ log2(baseline_titers), data = titerset_b)

# Set constants
baseline_constant <- log2(39)   # Correction level, log2(39) for seronegative as 39 is entered for <40
slope <- as.numeric(0.3063)             # slope of lin. regression

# 4)-------------------------------------------------------------------------------------------------------------------
# Apply the function
titerset_b_corrected <- titerset_b %>%
  mutate(
    outcome_titers_corrected = correct_T_post(outcome_titers, baseline_titers, baseline_constant, slope),
    outcome_titers_corrected_unlogged = 2^outcome_titers_corrected
  ) 

# look at result
b_titers_uncorrected<- ggplot(data = titerset_b_corrected, aes(y = log2(outcome_titers), x = log2(baseline_titers))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5, 13)) +
  scale_y_continuous(limits = c(4, 15)) +
  labs(x = "Baseline Titers (log2)",
       y = "Raw Outcome Titers (log2)",
       title = "B Victoria")  +
  theme_classic()

b_titers_corrected<- ggplot(data = titerset_b_corrected, aes(y = outcome_titers_corrected, x = log2(baseline_titers))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5, 13)) +
  scale_y_continuous(limits = c(4, 15)) +
  labs(x = "Baseline Titers (log2)",
       y = "Corrected Outcome Titers (log2)",
       title = "B Victoria") +
  theme_classic()


# multivar. model
b_titerset_b_corrected_meta <- titerset_b_corrected %>% 
  left_join(basefile_withPID)

b_model3 <- lm(outcome_titers_corrected ~ log2(baseline_titers), data = b_titerset_b_corrected_meta)
sum_b_model3 <- summary(b_model3) # +Age, sex (als factor), gruppe (als factor)
b_model4 <- lm(outcome_titers_corrected ~ log2(baseline_titers) + gen_age + as.factor(gen_sex) + as.factor(study_group), data = b_titerset_b_corrected_meta)
sum_b_model4 <- summary(b_model4) # +Age, sex (als factor), gruppe (als factor)



