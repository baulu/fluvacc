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

#####facet_grid by patient group
#preparation
microneut_analysis_raw_p <- microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50") %>%
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .)))  

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

plot_comb_group_wide <- plot_h1_wide / plot_h3_wide / plot_vic_wide



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
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
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
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
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
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
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
  mutate(across(c(ic50), ~ ifelse(is.na(.), 39, .))) %>%
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

plot_comb_strain_wide <- plot_co / plot_hiv / plot_rh / plot_ms / plot_onc

### plots with unifrom fixes y-axis
df <- microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50) %>% 
  pivot_longer(cols = c(H1, H3, Vic), names_to = "strain", values_to = "ic50") %>%
  mutate(across(ic50, ~ ifelse(is.na(.), 39, .))) %>%
  filter(pat_group == "Control")

# Define individual plots
plot_H1 <- df %>% filter(strain == "H1") %>%
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +
  geom_boxplot(outlier.shape = NA, alpha = 0) +
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +
  scale_y_log10(limits = c(10^1, 10^4), labels = scales::comma) +
  labs(title = "H1", y = "IC50 - Control", x = "Sampling Number") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

plot_H3 <- df %>% filter(strain == "H3") %>%
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +
  geom_boxplot(outlier.shape = NA, alpha = 0) +
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +
  scale_y_log10(limits = c(10^1, 10^4), labels = scales::comma) +
  labs(title = "H3", y = "IC50 - Control", x = "Sampling Number") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

plot_Vic <- df %>% filter(strain == "Vic") %>%
  ggplot(aes(x = factor(Sampling_number), y = ic50, fill = factor(Sampling_number))) +
  geom_boxplot(outlier.shape = NA, alpha = 0) +
  geom_jitter(width = 0.2, height = 0, color = "lightblue", size = 2, alpha = 0.5) +
  scale_y_log10(limits = c(10^1.5, 10^4), labels = scales::comma) +
  labs(title = "Vic", y = "IC50 - Control", x = "Sampling Number") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# Combine plots using patchwork
plot_H1 + plot_H3 + plot_Vic + plot_layout(ncol = 3)

###

#fold increase analysis
microneut_analysis_raw %>% 
  select(pat_group, Sampling_number, H1 = FluV_H1_fold, H3 = FluV_H3_fold, Vic = FluV_Vic_fold) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "fold") %>%  
  filter(Sampling_number == 2) %>% 
  select(pat_group, strain, fold) %>% print()
  filter(pat_group == "ONK") %>% 
  ggplot(aes(x = factor(strain), y = fold, fill = factor(strain))) +  # Factorize Sampling_number for distinct boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0) +  # Boxplot with transparency
  geom_jitter(width = 0.2, height = 0, color = "grey7", size = 2, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 10))
### Zusätzlich linie bei 1, 
### regression nach invers. baseline titer (fold change ln2 transformieren vorher und anschl, 2summieren, vgl.pnk-paper)
### adjusten unf auf MEdian baselinetiter (vgl. hirzel neurofilament paper)

h1_base <-  microneut_analysis_raw %>%
    filter(Sampling_number == 1) %>% 
    select(PID, H1_baseline = FluA_H1_ic50)

h1_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, H1_fold = FluV_H1_fold)
  
hi1_foldbase <- h1_base %>% 
  left_join(h1_fold) %>% 
  mutate(H1_fold_ln2 = log2(H1_fold))  %>% 
  mutate(H1_fold_ln10 = log(H1_fold))  %>% 
  as_tibble()

hi1_foldbase %>% 
  ggplot(aes(x = H1_baseline, y = H1_fold_ln10)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(0, 13000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()


