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
library(scales)
library(patchwork)

#load data
microneut_analysis_raw  <- read.xlsx("processed/microneut_analysis_raw.xlsx")

### Some look at the data

#facet_grid by patient group
microneut_analysis_raw %>% 
  #replace Non-inhibition values (=NA) with 39
  mutate(across(c(FluA_H1_ic50, FluA_H3_ic50, FluA_Vic_ic50), ~ ifelse(is.na(.), 39, .))) %>% 
  ggplot(aes(x = Sampling_number, y = FluA_Vic_ic50, color = pat_group)) +
  geom_jitter(width = 0.1, height = 0)+
  scale_y_log10(labels = comma)+
  facet_wrap(~ pat_group, scales = "free_y") +
  theme_classic()

#

plot_co <- microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50") %>%
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
  filter(pat_group == "Control") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - Control", x = "Sampling Number") +
  scale_y_log10(labels = scales::comma) +
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
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
  filter(pat_group == "HIV") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "lightgreen", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - HIV", x = "Sampling Number") +
  scale_y_log10(labels = scales::comma) +
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
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
  filter(pat_group == "RH") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "brown1", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - Rheuma", x = "Sampling Number") +
  scale_y_log10(labels = scales::comma) +
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
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
  filter(pat_group == "MS") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "deepskyblue3", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - MS", x = "Sampling Number") +
  scale_y_log10(labels = scales::comma) +
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
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
  filter(pat_group == "ONK") %>% 
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "grey7", size = 2, alpha = 0.5) +  # Add jittered points
  labs(y = "IC50 - CART", x = "Sampling Number") +
  scale_y_log10(labels = scales::comma) +
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

plot_comb <- plot_co / plot_hiv / plot_rh / plot_ms / plot_onc


#fold increase analysis
microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluV_H1_fold, H3 = FluV_H3_fold, Vic = FluV_Vic_fold) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "fold") %>%  
  filter(Sampling_number == 2) %>% 
  select(pat_group, strain, fold) %>% 
  filter(pat_group == "MS") %>% 
  ggplot(aes(x = factor(strain), y = fold, fill = factor(strain))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "grey7", size = 2, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 10))


#+  # Add jittered points
labs(y = "IC50 - CART", x = "Strain") +
  scale_y_log10(labels = scales::comma) +
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