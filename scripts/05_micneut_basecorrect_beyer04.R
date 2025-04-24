#Index
# 1) Load libraries and data 
# 2) Correction function of outcome titers to baselin titers (according to beyer 2004)
# 3) Create dataset for appliing correction function
# 4) Slopes and medians
# 5) B Victoria
# 6) H1N1
# 7) H3N2

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
  select(PID, baseline_titers = B_baseline, outcome_titers = B_outcome)

titerset_h1n1 <- titerset_combined %>% 
  select(PID, baseline_titers = H1_baseline, outcome_titers = H1_outcome)

titerset_h3n2 <- titerset_combined %>% 
  select(PID, baseline_titers = H3_baseline, outcome_titers = H3_outcome)

# 4)-------------------------------------------------------------------------------------------------------------------
#slopes and medians: 
# B: 0.3063  (from b_model_titers)
b_model_titers <- lm(log2(outcome_titers) ~ log2(baseline_titers), data = titerset_b)
b_slope <- coef(b_model_titers)["log2(baseline_titers)"] 
b_med <- median(titerset_b$baseline_titers)

# H1N1: 
h1n1_model_titers <- lm(log2(outcome_titers) ~ log2(baseline_titers), data = titerset_h1n1)
h1n1_slope <- coef(h1n1_model_titers)["log2(baseline_titers)"] 
h1n1_med <- median(titerset_h1n1$baseline_titers)

# H3N2: 
h3n2_model_titers <- lm(log2(outcome_titers) ~ log2(baseline_titers), data = titerset_h3n2)
h3n2_slope <- coef(h3n2_model_titers)["log2(baseline_titers)"] 
h3n2_med <- median(titerset_h3n2$baseline_titers)

# 5)-------------------------------------------------------------------------------------------------------------------
#### B Victoria Strain
slope <- b_slope            # slope of lin. regression
# Set constants for seronegative as baseline constant
baseline_constant_0 <- log2(39)   # Correction level, log2(39) for seronegative as 39 is entered for <40
baseline_constant_med <- log2(b_med) # Correction level, median titer for b-strain

# Apply the function with seronegative as baseline constant
titerset_b_corrected_seroneg <- titerset_b %>%
  mutate(
    outcome_titers_corrected = correct_T_post(outcome_titers, baseline_titers, baseline_constant_0, slope),
    outcome_titers_corrected_unlogged = 2^outcome_titers_corrected
  ) 

titerset_b_corrected_med <- titerset_b %>%
  mutate(
    outcome_titers_corrected = correct_T_post(outcome_titers, baseline_titers, baseline_constant_med, slope),
    outcome_titers_corrected_unlogged = 2^outcome_titers_corrected
  ) 

# 6)-------------------------------------------------------------------------------------------------------------------
#### H1N1 Victoria Strain
slope <- h1n1_slope            # slope of lin. regression
# Set constants for seronegative as baseline constant
baseline_constant_0 <- log2(39)   # Correction level, log2(39) for seronegative as 39 is entered for <40
baseline_constant_med <- log2(h1n1_med) # Correction level, median titer for b-strain

# Apply the function with seronegative as baseline constant
titerset_h1n1_corrected_seroneg <- titerset_h1n1 %>%
  mutate(
    outcome_titers_corrected = correct_T_post(outcome_titers, baseline_titers, baseline_constant_0, slope),
    outcome_titers_corrected_unlogged = 2^outcome_titers_corrected
  ) 

titerset_h1n1_corrected_med <- titerset_h1n1 %>%
  mutate(
    outcome_titers_corrected = correct_T_post(outcome_titers, baseline_titers, baseline_constant_med, slope),
    outcome_titers_corrected_unlogged = 2^outcome_titers_corrected
  ) 

# look at results
# seronegative-corrected
h1n1_titers_uncorrected_seroneg <- ggplot(data = titerset_h1n1_corrected_seroneg, aes(y = log2(outcome_titers), x = log2(baseline_titers))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5.5, 16)) +
  scale_y_continuous(limits = c(4, 17.5)) +
  labs(x = "Baseline Titers (log2)",
       y = "Raw Outcome Titers (log2)",
       title = "H1N1")  +
  theme_classic()

h1n1_titers_corrected_seroneg <- ggplot(data = titerset_h1n1_corrected_seroneg, aes(y = outcome_titers_corrected, x = log2(baseline_titers))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5.5, 16)) +
  scale_y_continuous(limits = c(4, 17.5)) +
  labs(x = "Baseline Titers (log2)",
       y = "Corrected Outcome Titers (log2)",
       title = "H1N1") +
  theme_classic()

# 7)-------------------------------------------------------------------------------------------------------------------
#### H3N2 Victoria Strain
slope <- h3n2_slope            # slope of lin. regression
# Set constants for seronegative as baseline constant
baseline_constant_0 <- log2(39)   # Correction level, log2(39) for seronegative as 39 is entered for <40
baseline_constant_med <- log2(h3n2_med) # Correction level, median titer for b-strain

# Apply the function with seronegative as baseline constant
titerset_h3n2_corrected_seroneg <- titerset_h3n2 %>%
  mutate(
    outcome_titers_corrected = correct_T_post(outcome_titers, baseline_titers, baseline_constant_0, slope),
    outcome_titers_corrected_unlogged = 2^outcome_titers_corrected
  ) 

titerset_h3n2_corrected_med <- titerset_h3n2 %>%
  mutate(
    outcome_titers_corrected = correct_T_post(outcome_titers, baseline_titers, baseline_constant_med, slope),
    outcome_titers_corrected_unlogged = 2^outcome_titers_corrected
  ) 

# look at results
# seronegative-corrected
h3n2_titers_uncorrected_seroneg <- ggplot(data = titerset_h3n2_corrected_seroneg, aes(y = log2(outcome_titers), x = log2(baseline_titers))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5, 13)) +
  scale_y_continuous(limits = c(4, 16)) +
  labs(x = "Baseline Titers (log2)",
       y = "Raw Outcome Titers (log2)",
       title = "H3N2")  +
  theme_classic()

h3n2_titers_corrected_seroneg <- ggplot(data = titerset_h3n2_corrected_seroneg, aes(y = outcome_titers_corrected, x = log2(baseline_titers))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5, 13)) +
  scale_y_continuous(limits = c(4, 16)) +
  labs(x = "Baseline Titers (log2)",
       y = "Corrected Outcome Titers (log2)",
       title = "H3N2") +
  theme_classic()

# 8)-------------------------------------------------------------------------------------------------------------------
# Combining dataset
titerset_comb_b_corrected_seroneg <- titerset_b_corrected_med %>% 
  #mutate(baseline_titers_seroneg = 39.000) %>% 
  mutate(fold_change = outcome_titers / baseline_titers) %>% 
  mutate(fold_change_corrbase_meyer = outcome_titers_corrected_unlogged / baseline_titers) %>% 
  select(PID, baseline_titers, outcome_titers, outcome_titers_corrected_unlogged, fold_change, fold_change_corrbase_meyer) 
 # FC <1 are still in -> problematic!

