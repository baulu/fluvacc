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
microneut_analysis_raw  <- read.xlsx("processed/microneut_analysis_raw.xlsx")
hai_analysis_raw  <- read.xlsx("processed/hai_analysis_raw.xlsx")
influenza_antibody_results <- read.xlsx("processed/influenza_antibody_results.xlsx")

### Firs look at HAI data
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



##### First look at microneutralisation data
#####facet_grid by patient group
#preparation
microneut_analysis_raw_p <- microneut_analysis_raw %>% 
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

plot_comb_strain_wide <- plot_co / plot_hiv / plot_rh / plot_ms / plot_onc

######put together dataset with fold changes and baseline titer for all strains
###h1
h1_base <-  microneut_analysis_raw %>%
    filter(Sampling_number == 1) %>% 
    select(PID, H1_baseline = FluA_H1_ic50)

h1_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, H1_fold = FluV_H1_fold)
  
hi1_foldbase <- h1_base %>% 
  left_join(h1_fold) %>% 
  mutate(H1_fold_log2 = log2(H1_fold))  %>% 
  mutate(H1_fold_log10 = log10(H1_fold))  %>% 
  mutate(H1_base_log2 = log2(H1_baseline))  %>% 
  mutate(H1_base_log10 = log10(H1_baseline))  %>% 
  as_tibble()

hi1_linreg_plot_baselog2 <- hi1_foldbase %>% 
  ggplot(aes(x = H1_base_log2, y = H1_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  #scale_x_continuous(limits = c(0, 13000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4)  + # Display Pearson's r
  labs(y = "H1 - log2 Fold-Change", x = NULL) 

hi1_linreg_plot_basecont <- hi1_foldbase %>% 
  ggplot(aes(x = H1_baseline, y = H1_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(0, 13000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = NULL, x = NULL) 


# Perform Pearson correlation tests
cor_result1 <- cor.test(hi1_foldbase$H1_baseline, hi1_foldbase$H1_fold_log2, method = "pearson")
cor_result2 <- cor.test(hi1_foldbase$H1_base_log2, hi1_foldbase$H1_fold_log2, method = "pearson")

# Extract values and save in a dataframe
cor_df <- data.frame(
  Variable1 = "H1_baseline",
  Variable2 = "H1_fold_log2",
  Correlation_Coefficient = cor_result1$estimate,
  P_Value = cor_result1$p.value
)

cor_df2 <- data.frame(
  Variable1 = "H1_base_log2",
  Variable2 = "H1_fold_log2",
  Correlation_Coefficient = cor_result2$estimate,
  P_Value = cor_result1$p.value
)

cor_comb <- cor_df %>% 
  rows_append(cor_df2) %>% 
  print()
 

###h3
h3_base <-  microneut_analysis_raw %>%
  filter(Sampling_number == 1) %>% 
  select(PID, H3_baseline = FluA_H3_ic50)

h3_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, H3_fold = FluV_H3_fold)

hi3_foldbase <- h3_base %>% 
  left_join(h3_fold) %>% 
  mutate(H3_fold_log2 = log2(H3_fold))  %>% 
  mutate(H3_fold_log10 = log10(H3_fold))  %>% 
  mutate(H3_base_log2 = log2(H3_baseline))  %>% 
  mutate(H3_base_log10 = log10(H3_baseline))  %>% 
  as_tibble()

hi3_linreg_plot_baselog2 <- hi3_foldbase %>% 
  ggplot(aes(x = H3_base_log2, y = H3_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  #scale_x_continuous(limits = c(0, 13000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = "H3 - log2 Fold-Change", x = "log2 Transformed Baseline Titers") 
  
hi3_linreg_plot_basecont <- hi3_foldbase %>% 
  ggplot(aes(x = H3_baseline, y = H3_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = NULL, x = "Continuous Baseline Titers") 

# Perform Pearson correlation tests
corh3_result1 <- cor.test(hi3_foldbase$H3_baseline, hi3_foldbase$H3_fold_log2, method = "pearson")
corh3_result2 <- cor.test(hi3_foldbase$H3_base_log2, hi3_foldbase$H3_fold_log2, method = "pearson")

# Extract values and save in a dataframe
corh3_df <- data.frame(
  Variable1 = "H3_baseline",
  Variable2 = "H3_fold_log2",
  Correlation_Coefficient = corh3_result1$estimate,
  P_Value = corh3_result1$p.value
)

corh3_df2 <- data.frame(
  Variable1 = "H3_base_log2",
  Variable2 = "H3_fold_log2",
  Correlation_Coefficient = corh3_result2$estimate,
  P_Value = corh3_result2$p.value
)

corh3_comb <- corh3_df %>% 
  rows_append(corh3_df2) %>% 
  print()

###B
b_base <-  microneut_analysis_raw %>%
  filter(Sampling_number == 1) %>% 
  select(PID, Vic_baseline = FluA_Vic_ic50)

b_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, Vic_fold = FluV_Vic_fold)

b_foldbase <- b_base %>% 
  left_join(b_fold) %>% 
  mutate(Vic_fold_log2 = log2(Vic_fold))  %>% 
  mutate(Vic_fold_log10 = log10(Vic_fold))  %>% 
  mutate(Vic_base_log2 = log2(Vic_baseline))  %>% 
  mutate(Vic_base_log10 = log10(Vic_baseline))  %>% 
  as_tibble()

b_linreg_plot_baselog2 <- b_foldbase %>% 
  ggplot(aes(x = Vic_base_log2, y = Vic_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = "B - log2 Fold-Change", x = NULL) 


### multivariables lineares model
model1 <- lm(Vic_fold_log2 ~ Vic_base_log2, data = b_foldbase)
summary(model1). # +Age, sex (als factor), gruppe (als factor)
qqnorm(model1)

#Vergleiche --> adjustet R2, Residuals, AIC
#model2 <- lm(#titer_outc_log# ~ Vic_base_log2, data = b_foldbase)
#summary(model1). # +Age, sex (als factor), gruppe (als factor)
#qqnorm(model1) #qqplot vgl.

## repeated measures model


# residuals plot, qq plot
b_linreg_plot_basecont <- b_foldbase %>% 
  ggplot(aes(x = Vic_baseline, y = Vic_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(0, 20000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = NULL, x = NULL) 

# Perform Pearson correlation tests
corb_result1 <- cor.test(b_foldbase$Vic_baseline, b_foldbase$Vic_fold_log2, method = "pearson")
corb_result2 <- cor.test(b_foldbase$Vic_base_log2, b_foldbase$Vic_fold_log2, method = "pearson")

# Extract values and save in a dataframe
corb_df <- data.frame(
  Variable1 = "B_baseline",
  Variable2 = "B_fold_log2",
  Correlation_Coefficient = corb_result1$estimate,
  P_Value = corb_result1$p.value
)

corb_df2 <- data.frame(
  Variable1 = "B_base_log2",
  Variable2 = "B_fold_log2",
  Correlation_Coefficient = corb_result2$estimate,
  P_Value = corb_result2$p.value
)

corb_comb <- corb_df %>% 
  rows_append(corb_df2) %>% 
  print()

###combining lin. regression plots
linreg_plot_comb = (b_linreg_plot_baselog2 + b_linreg_plot_basecont) / (hi1_linreg_plot_baselog2 + hi1_linreg_plot_basecont) / (hi3_linreg_plot_baselog2 + hi3_linreg_plot_basecont)

# linear regression of baseline titer only in controls - change baseline-dataset for analysis before
#microneut_analysis_raw <- microneut_analysis_raw %>% 
#  filter(pat_group == "Control")


#### Checking if homoscedasticity is given (constant scattering of residuals independant of baseline)
###H3-Strain
# Remove NAs to ensure matching rows
h3_foldbase_clean <- hi3_foldbase %>% drop_na(H3_base_log2, H3_fold_log2)

# Fit linear model
adj_model_h3 <- lm(H3_fold_log2 ~ H3_base_log2, data = h3_foldbase_clean)

# Check for homoscedasticity using residuals vs fitted values plot
h3_residuals_plot <- ggplot(data = h3_foldbase_clean, aes(x = fitted(adj_model_h3), y = resid(adj_model_h3))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = "Fitted Values", y = "Residuals") +
  #ggtitle("H3 Residuals vs Fitted Values") +
  theme_minimal()

# Perform Breusch-Pagan test to check for homoscedasticity
durbinWatsonTest(adj_model_h3)
bp_test_h3 <- ncvTest(adj_model_h3)
print(bp_test_h3)

# QQ-Plot of residuals from a regression model
h3_qq_plot <- ggplot(data = data.frame(resid = residuals(adj_model_h3)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

###H3-Strain
# Remove NAs to ensure matching rows
h1_foldbase_clean <- hi1_foldbase %>% drop_na(H1_base_log2, H1_fold_log2)

# Fit linear model
adj_model_h1 <- lm(H1_fold_log2 ~ H1_base_log2, data = h1_foldbase_clean)

# Check for homoscedasticity using residuals vs fitted values plot
h1_residuals_plot <- ggplot(data = h1_foldbase_clean, aes(x = fitted(adj_model_h1), y = resid(adj_model_h1))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = "Fitted Values", y = "Residuals") +
  #ggtitle("H1 Residuals vs Fitted Values") +
  theme_minimal()

# Perform Breusch-Pagan test to check for homoscedasticity
durbinWatsonTest(adj_model_h1)
bp_test_h1 <- ncvTest(adj_model_h1)
print(bp_test_h1)

# QQ-Plot of residuals from a regression model
h1_qq_plot <- ggplot(data = data.frame(resid = residuals(adj_model_h1)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

###B-Strain
# Remove NAs to ensure matching rows
b_foldbase_clean <- b_foldbase %>% drop_na(Vic_base_log2, Vic_fold_log2)

# Fit linear model
adj_model_b <- lm(Vic_fold_log2 ~ Vic_base_log2, data = b_foldbase_clean)

# Check for homoscedasticity using residuals vs fitted values plot
b_residuals_plot <- ggplot(data = b_foldbase_clean, aes(x = fitted(adj_model_b), y = resid(adj_model_b))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = "Fitted Values", y = "Residuals") +
  ggtitle("Residuals vs Fitted Values") +
  theme_minimal()

# Perform Breusch-Pagan test to check for homoscedasticity
durbinWatsonTest(adj_model_b)
bp_test_b <- ncvTest(adj_model_b)
print(bp_test_b)


# QQ-Plot of residuals from a regression model
b_qq_plot <- ggplot(data = data.frame(resid = residuals(adj_model_b)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(title = "QQ-Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles")


###combining regression plots
linreg_plot_test_combined = ((b_linreg_plot_baselog2 + b_residuals_plot + b_qq_plot) /
                      (hi1_linreg_plot_baselog2 + h1_residuals_plot + h1_qq_plot) / 
                      (hi3_linreg_plot_baselog2 + h3_residuals_plot + h3_qq_plot))

#### Comparing HAI and Microneutralisation assay
#harmonize dataset for analysis - same vraiables, row_append hai
influenza_antibody_harmonised <- microneut_analysis_raw_p %>% 
  mutate(strain = case_when( strain == "H1" ~ "H1N1",
                             strain == "H3" ~ "H3N2",
                             strain == "Vic" ~ "BVic")) %>% 
  mutate(type = "ic50_microneutr") %>% 
  select(pat_group, Sampling_number, strain, result = ic50, type, PID) %>%  
  rows_append(hai_analysis_raw_p %>% mutate(type = "hai")) %>% print()

  
#rank values by strain and type(Hai vs. Microneut) for correlation - ranking only within the same strain
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
  
ggplot(data=df3, aes(x=hai_rank, y=ic50_rank)) +
  geom_point(size = 2, alpha = 0.5)+
  geom_smooth(method = "lm", se = TRUE)   # Add regression line with confidence interval

spearman <- cor.test(df3$hai_rank, df3$ic50_rank, method = "spearman")
rho <- round(spearman$estimate, 2)
pval <- format.pval(spearman$p.value, digits = 3, eps = .001)


ggplot(data=df3, aes(x = hai_rank, y = ic50_rank)) +  # Factorize Sampling_number for distinct boxplots
  geom_jitter( size = 2, alpha = 0.5) +  # Add jittered points
  geom_smooth(method = "lm", se = TRUE) +   # Add regression line with confidence interval
  annotate("text",
           x = Inf, y = -Inf,
           hjust = 1.1, vjust = -0.7,
           label = paste0("Spearman ρ = ", rho, "\nP = ", pval),
           size = 4) +
  labs(title = "Rangsummenkorrelation: HAI vs IC50",
       x = "HAI Rank",
       y = "IC50 Rank") +
  theme_minimal()
  
#alternative approach using corr (ranking across strains)
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
  
corr <- cor.test(x=tf3$haires, tf3$ic50res, method = 'spearman')
  
