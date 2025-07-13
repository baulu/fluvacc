#Index
# 1) Load libraries and data 
# 2) HAI Analysis
# 3) Microneutralisation Titers

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
library(ggpubr)
library(car)

#load data
microneut_analysis_raw  <- read.xlsx(here::here("processed", "microneut_analysis_raw.xlsx"))
hai_analysis_raw  <- read.xlsx(here::here("processed", "hai_analysis_raw.xlsx"))
influenza_antibody_results <- read.xlsx(here::here("processed", "influenza_antibody_results.xlsx"))
basefile_withPID <- read.xlsx(here::here("processed", "FluVac_basefile_withPID.xlsx"))basecorrect_categories <- read_csv("data/summary_classification.csv")
basecorrect_categories <- read_csv("data/summary_classification.csv")

# 2) -------------------------------------------------------------------------------------------------------------------
### First look at HAI data
#####facet_grid by patient group
#preparation
hai_analysis_raw_p <- hai_analysis_raw %>% 
  select(pat_group, Sampling_number, H1N1, H3N2, BVic, BYam, PID) %>% 
  pivot_longer(cols = 3:6, names_to = "strain", values_to = "result")  

plot_hai_h1n1_wide <- hai_analysis_raw_p %>% 
  select(pat_group, strain, Sampling_number, result) %>% 
  filter(strain == "H1N1") %>% 
  ggplot(aes(x = factor(Sampling_number), y = result, color = pat_group)) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "HAI - H1N1", x = "Sampling Number") +
  scale_y_log10(limits = c(10^0.5, 10^3.5), labels = scales::comma) +
  facet_wrap(~ pat_group, nrow = 1) +  # Facet by patgroup
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    strip.text = element_text(face = "bold"),  # Make facet labels bold
    legend.position = "none"  # Remove legend
  )

plot_hai_h3n2_wide <- hai_analysis_raw_p %>% 
  select(pat_group, strain, Sampling_number, result) %>% 
  filter(strain == "H3N2") %>% 
  ggplot(aes(x = factor(Sampling_number), y = result, color = pat_group)) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "HAI - H3N2", x = "Sampling Number") +
  scale_y_log10(limits = c(10^0.5, 10^3.5), labels = scales::comma) +
  facet_wrap(~ pat_group, nrow = 1) +  # Facet by patgroup
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    strip.text = element_blank(), 
    legend.position = "none"  # Remove legend
  )

plot_hai_bvic_wide <- hai_analysis_raw_p %>% 
  select(pat_group, strain, Sampling_number, result) %>% 
  filter(strain == "BVic") %>% 
  ggplot(aes(x = factor(Sampling_number), y = result, color = pat_group)) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "HAI - BVic", x = "Sampling Number") +
  scale_y_log10(limits = c(10^0.5, 10^3.5), labels = scales::comma) +
  facet_wrap(~ pat_group, nrow = 1) +  # Facet by patgroup
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    strip.text = element_blank(), 
    legend.position = "none"  # Remove legend
  )


plot_hai_byam_wide <- hai_analysis_raw_p %>% 
  select(pat_group, strain, Sampling_number, result) %>% 
  filter(strain == "BYam") %>% 
  ggplot(aes(x = factor(Sampling_number), y = result, color = pat_group)) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "HAI - BYam", x = "Sampling Number") +
  scale_y_log10(limits = c(10^0.5, 10^3.5), labels = scales::comma) +
  facet_wrap(~ pat_group, nrow = 1) +  # Facet by patgroup
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.text.x = element_text(face = "bold"),  # Make y-axis title bold
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    strip.text = element_blank(), 
    legend.position = "none"  # Remove legend
  )

plot_hai_comb_group_wide <- plot_hai_h1n1_wide / plot_hai_h3n2_wide / plot_hai_bvic_wide / plot_hai_byam_wide

# 3) -------------------------------------------------------------------------------------------------------------------
##### First look at microneutralisation data
#####facet_grid by patient group
#preparation
microneut_analysis_raw_p <- microneut_analysis_raw %>% #CAVE: USED IN OTHER SCRIPTS - e.g. 04_hai_vs_mic...
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50, PID) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50")  

plot_h1_wide <- microneut_analysis_raw_p %>% 
  select(pat_group, strain, Sampling_number, ic50) %>% 
  filter(strain == "H1") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, color = pat_group)) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - H1", x = "Sampling Number") +
  scale_y_log10(limits = c(10^1.5, 10^5.2), labels = scales::comma) +
  facet_wrap(~ pat_group, nrow = 1) +  # Facet by patgroup
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    strip.text = element_text(face = "bold"),  # Make facet labels bold
    legend.position = "none"  # Remove legend
  )

plot_h3_wide <- microneut_analysis_raw_p %>% 
  select(pat_group, strain, Sampling_number, ic50) %>% 
  filter(strain == "H3") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, color = pat_group)) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "pink", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - H3", x = "Sampling Number") +
  scale_y_log10(limits = c(10^1.5, 10^5.2), labels = scales::comma) +
  facet_wrap(~ pat_group, nrow = 1) +  # Facet by patgroup
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    strip.text = element_blank(),
    legend.position = "none"  # Remove legend
  )


plot_vic_wide <- microneut_analysis_raw_p %>% 
  select(pat_group, strain, Sampling_number, ic50) %>% 
  filter(strain == "Vic") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, color = pat_group)) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightgreen", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - Vic", x = "Sampling Number") +
  scale_y_log10(limits = c(10^1.5, 10^5.2), labels = scales::comma) +
  facet_wrap(~ pat_group, nrow = 1) +  # Facet by patgroup
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    axis.text.x = element_text(face = "bold"),  # Make y-axis title bold
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"  # Remove legend
  )

plot_ic50_comb_group_wide <- plot_h1_wide / plot_h3_wide / plot_vic_wide


#facet_grid by strain
plot_co <- microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50") %>%
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
  filter(pat_group == "Control") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - Control", x = "Sampling Number") +
  scale_y_log10(limits = c(10^1.5, 10^5.2), labels = scales::comma) +
  facet_wrap(~ strain, scales = "free_y") +  # Facet by strain
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    strip.text = element_text(face = "bold"),  # Make facet labels bold
    legend.position = "none"  # Remove legend
  )


plot_hiv <- microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50") %>%
  filter(pat_group == "HIV") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightgreen", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - HIV", x = "Sampling Number") +
  scale_y_log10(limits = c(10^1.5, 10^5.2), labels = scales::comma) +
  facet_wrap(~ strain, scales = "free_y") +  # Facet by strain
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    legend.position = "none"  # Remove legend
  )


plot_rh <- microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50") %>%
  filter(pat_group == "RH") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "brown1", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - Rheuma", x = "Sampling Number") +
  scale_y_log10(limits = c(10^1.5, 10^5.2), labels = scales::comma) +
  facet_wrap(~ strain, scales = "free_y") +  # Facet by strain
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    legend.position = "none"  # Remove legend
  )

plot_ms <- microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50") %>%
  filter(pat_group == "MS") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "deepskyblue3", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - MS", x = "Sampling Number") +
  scale_y_log10(limits = c(10^1.5, 10^5), labels = scales::comma) +
  facet_wrap(~ strain, scales = "free_y") +  # Facet by strain
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),  # Make y-axis title bold
    legend.position = "none"  # Remove legend
  )

plot_onc <- microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50") %>%
  filter(pat_group == "ONK") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "grey7", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - CART", x = "Sampling Number") +
  scale_y_log10(limits = c(10^1.5, 10^5.2), labels = scales::comma) +
  facet_wrap(~ strain, scales = "free_y") +  # Facet by strain
  scale_x_discrete(limits = c("1", "2"), labels = c("preVac", "1m")) +  # Ensure the x-axis shows 1 and 2
  theme_minimal() +
  guides(color = "none") +
  theme(
    axis.text.x = element_text(face = "bold"),  
    axis.title.y = element_text(face = "bold"),  
    axis.title.x = element_text(face = "bold"),  
    strip.text = element_blank(),  # Optionally remove facet labels
    legend.position = "none"  # Remove legend
  )

plot_ic50_comb_strain_wide <- plot_co / plot_hiv / plot_rh / plot_ms / plot_onc

