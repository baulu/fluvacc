#Index
# 1) Load libraries and data 
# 2) Correlation of baseline titers with fold-changes of microneutralisation titers
    # 2b) lineare regression, modelle
# 3) Correlation of HAI and microneutralisation assays

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
library(knitr)
library(kableExtra)

#load data
microneut_analysis_raw  <- read.xlsx(here::here("processed", "microneut_analysis_raw.xlsx"))
hai_analysis_raw  <- read.xlsx(here::here("processed", "hai_analysis_raw.xlsx"))
influenza_antibody_results <- read.xlsx(here::here("processed", "influenza_antibody_results.xlsx"))
basefile_withPID <- read.xlsx(here::here("processed", "FluVac_basefile_withPID.xlsx"))

# 2) -------------------------------------------------------------------------------------------------------------------
##### assess correlation of microneutralisation fold changes with baseline titer
#put together dataset with fold changes and baseline titer for all strains
#h1
h1_base <-  microneut_analysis_raw %>%
    filter(Sampling_number == 1) %>% 
    select(PID, H1_baseline = FluA_H1_ic50)

h1_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, H1_fold = FluV_H1_fold, H1_fold_log2first = FluV_H1_fold_log2first, H1_fold_log10first = FluV_H1_fold_log10first)
  
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
  scale_x_continuous(limits = c(5, 16)) +
  scale_y_continuous(limits = c(-5, 8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4)  + # Display Pearson's r
  labs(y = "H1 - log2 FC", x = NULL) 

hi1_linreg_plot_basecont <- hi1_foldbase %>% 
  ggplot(aes(x = H1_baseline, y = H1_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(0, 13000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = NULL, x = NULL) 

hi1_linreg_plot_log2first <- hi1_foldbase %>% 
  ggplot(aes(x = H1_base_log2, y = H1_fold_log2first)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5, 16)) +
  scale_y_continuous(limits = c(0, 3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = 0.25) + # Display Pearson's r
  labs(y = "H1 - log2 first FC", x = NULL) 

# Perform Pearson correlation tests
cor_result1 <- cor.test(hi1_foldbase$H1_baseline, hi1_foldbase$H1_fold_log2, method = "pearson")
cor_result2 <- cor.test(hi1_foldbase$H1_base_log2, hi1_foldbase$H1_fold_log2, method = "pearson")
cor_result3_logfirst <- cor.test(hi1_foldbase$H1_base_log2, hi1_foldbase$H1_fold_log2first, method = "pearson")
cor_result4_log10first <- cor.test(hi1_foldbase$H1_base_log10, hi1_foldbase$H1_fold_log10first, method = "pearson")

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

cor_df3_log2first <- data.frame(
  Variable1 = "H1_base_log2",
  Variable2 = "H1_fold_log2first",
  Correlation_Coefficient = cor_result3_logfirst$estimate,
  P_Value = cor_result3_logfirst$p.value
)

cor_df4_log10first <- data.frame(
  Variable1 = "H1_base_log10",
  Variable2 = "H1_fold_log10first",
  Correlation_Coefficient = cor_result4_log10first$estimate,
  P_Value = cor_result4_log10first$p.value
)

cor_comb1 <- cor_df %>% 
  rows_append(cor_df2) %>% 
  rows_append(cor_df3_log2first) %>% 
  rows_append(cor_df4_log10first) %>% 
  print()
 

###h3
h3_base <-  microneut_analysis_raw %>%
  filter(Sampling_number == 1) %>% 
  select(PID, H3_baseline = FluA_H3_ic50)

h3_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, H3_fold = FluV_H3_fold, H3_fold_log2first = FluV_H3_fold_log2first, H3_fold_log10first = FluV_H3_fold_log10first)

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
  scale_x_continuous(limits = c(5, 16)) +
  scale_y_continuous(limits = c(-5, 8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = "H3 - log2 FC", x = "log2 Transformed Baseline Titers") 
  
hi3_linreg_plot_basecont <- hi3_foldbase %>% 
  ggplot(aes(x = H3_baseline, y = H3_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = NULL, x = "Continuous Baseline Titers") 

hi3_linreg_plot_log2first <- hi3_foldbase %>% 
  ggplot(aes(x = H3_base_log2, y = H3_fold_log2first)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5, 16)) +
  scale_y_continuous(limits = c(0, 3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = 0.25) + # Display Pearson's r
  labs(y = "H3 - log2 first FC", x = "log2 Transformed Baseline Titers") 

# Perform Pearson correlation tests
corh3_result1 <- cor.test(hi3_foldbase$H3_baseline, hi3_foldbase$H3_fold_log2, method = "pearson")
corh3_result2 <- cor.test(hi3_foldbase$H3_base_log2, hi3_foldbase$H3_fold_log2, method = "pearson")
corh3_result3_log2first <- cor.test(hi3_foldbase$H3_base_log2, hi3_foldbase$H3_fold_log2first, method = "pearson")
corh3_result4_log10first <- cor.test(hi3_foldbase$H3_base_log10, hi3_foldbase$H3_fold_log10first, method = "pearson")


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

corh3_df2_log2first <- data.frame(
  Variable1 = "H3_base_log2",
  Variable2 = "H3_fold_log2first",
  Correlation_Coefficient = corh3_result3_log2first$estimate,
  P_Value = corh3_result3_log2first$p.value
)

corh3_df2_log10first <- data.frame(
  Variable1 = "H3_base_log10",
  Variable2 = "H3_fold_log10first",
  Correlation_Coefficient = corh3_result4_log10first$estimate,
  P_Value = corh3_result4_log10first$p.value
)

corh3_comb <- corh3_df %>% 
  rows_append(corh3_df2) %>% 
  rows_append(corh3_df2_log2first) %>% 
  rows_append(corh3_df2_log10first) %>% 
  print()

###B
b_base <-  microneut_analysis_raw %>%
  filter(Sampling_number == 1) %>% 
  select(PID, Vic_baseline = FluA_Vic_ic50)

b_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, Vic_fold = FluV_Vic_fold, Vic_fold_log2first = FluV_Vic_fold_log2first, Vic_fold_log10first = FluV_Vic_fold_log10first)

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
  scale_x_continuous(limits = c(5, 16)) +
  scale_y_continuous(limits = c(-5, 8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = "B - log2 FC", x = NULL) 

b_linreg_plot_basecont <- b_foldbase %>% 
  ggplot(aes(x = Vic_baseline, y = Vic_fold_log2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(0, 20000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = NULL, x = NULL) 

b_linreg_plot_log2first <- b_foldbase %>% 
  ggplot(aes(x = Vic_base_log2, y = Vic_fold_log2first)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  scale_x_continuous(limits = c(5, 16)) +
  scale_y_continuous(limits = c(0, 3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = 0.25) + # Display Pearson's r
  labs(y = "B - log2 first FC", x = NULL) 

# Perform Pearson correlation tests
corb_result1 <- cor.test(b_foldbase$Vic_baseline, b_foldbase$Vic_fold_log2, method = "pearson")
corb_result2 <- cor.test(b_foldbase$Vic_base_log2, b_foldbase$Vic_fold_log2, method = "pearson")
corb_result3_log2first <- cor.test(b_foldbase$Vic_base_log2, b_foldbase$Vic_fold_log2first, method = "pearson")
corb_result4_log2first <- cor.test(b_foldbase$Vic_base_log10, b_foldbase$Vic_fold_log10first, method = "pearson")

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

corb_df3_log2first <- data.frame(
  Variable1 = "B_base_log2",
  Variable2 = "B_fold_log2first",
  Correlation_Coefficient = corb_result3_log2first$estimate,
  P_Value = corb_result3_log2first$p.value
)

corb_df3_log10first <- data.frame(
  Variable1 = "B_base_log10",
  Variable2 = "B_fold_log10first",
  Correlation_Coefficient = corb_result4_log2first$estimate,
  P_Value = corb_result4_log2first$p.value
)

corb_comb <- corb_df %>% 
  rows_append(corb_df2) %>%
  rows_append(corb_df3_log2first) %>% 
  rows_append(corb_df3_log10first) %>% 
  print()

#combining all correlation tests
corb_comb_all <- corb_comb %>% 
  rows_append(corh3_comb) %>% 
  rows_append(cor_comb1) %>% 
  print()

###combining lin. regression plots
linreg_plot_comb = (b_linreg_plot_baselog2 + b_linreg_plot_basecont) / (hi1_linreg_plot_baselog2 + hi1_linreg_plot_basecont) / (hi3_linreg_plot_baselog2 + hi3_linreg_plot_basecont)
linreg_plot_comb2 = (b_linreg_plot_baselog2 + b_linreg_plot_log2first) / (hi1_linreg_plot_baselog2 + hi1_linreg_plot_log2first) / (hi3_linreg_plot_baselog2 + hi3_linreg_plot_log2first)  


# 2b)) -------------------------------------------------------------------------------------------------------------------
### Linere regression log2 baseline titer vs. log2 fold-change and qq-plots
## univariables lineares model / B-Victoria
b_model1 <- lm(Vic_fold_log2 ~ Vic_base_log2, data = b_foldbase)
sum_b_model1 <- summary(b_model1)

#qq plot
b_uni_qq_plot <- ggplot(data = data.frame(resid = residuals(b_model1)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "B univar Q-Q Plot",
          x = "Theoretical Quantiles", 
          y = "Sample Quantiles")+
  theme(plot.title = element_text(hjust = 0.5))

# log first
b_model1_lf <- lm(Vic_fold_log2first ~ Vic_base_log2, data = b_foldbase)
sum_b_model1_lf <- summary(b_model1_lf)

#qq plot
b_uni_lf_qq_plot <- ggplot(data = data.frame(resid = residuals(b_model1_lf)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "B univar Q-Q Plot / LogF",
          x = "Theoretical Quantiles", 
          y = NULL)+
  theme(plot.title = element_text(hjust = 0.5))

## multivariables lineares model / B-Victoria
b_foldbase_meta <- b_foldbase %>% 
  left_join(basefile_withPID)

b_model2 <- lm(Vic_fold_log2 ~ Vic_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = b_foldbase_meta)
summary(b_model2) # +Age, sex (als factor), gruppe (als factor)

#qq plot
b_multi_qq_plot <- ggplot(data = data.frame(resid = residuals(b_model2)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "B multivariant Q-Q Plot",
          x = "Theoretical Quantiles", 
          y = "Sample Quantiles")+
  theme(plot.title = element_text(hjust = 0.5))

# log first
b_model2_lf <- lm(Vic_fold_log2first ~ Vic_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = b_foldbase_meta)
summary(b_model2_lf) # +Age, sex (als factor), gruppe (als factor)

#qq plot
b_multi_lf_qq_plot <- ggplot(data = data.frame(resid = residuals(b_model2_lf)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "B multivariant Q-Q Plot / LogF",
          x = "Theoretical Quantiles", 
          y = NULL)+
  theme(plot.title = element_text(hjust = 0.5))

## univeriables lineares model / H1N1
h1_model1 <- lm(H1_fold_log2 ~ H1_base_log2, data = hi1_foldbase)
sum_h1_model1 <- summary(h1_model1)

#qq plot
h1_uni_qq_plot <- ggplot(data = data.frame(resid = residuals(h1_model1)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "H1N1 univar Q-Q Plot",
          x = NULL, 
          y = "Sample Quantiles")+
  theme(plot.title = element_text(hjust = 0.5))

# log first
h1_model1_lf <- lm(H1_fold_log2first ~ H1_base_log2, data = hi1_foldbase)
sum_h1_model1_lf <- summary(h1_model1)

#qq plot
h1_uni_lf_qq_plot <- ggplot(data = data.frame(resid = residuals(h1_model1_lf)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "H1N1 univar Q-Q Plot / logF",
          x = NULL, 
          y = NULL)+
  theme(plot.title = element_text(hjust = 0.5))

## multivariables lineares model / H1N1
hi1_foldbase_meta <- hi1_foldbase %>% 
  left_join(basefile_withPID)

h1_model2 <- lm(H1_fold_log2 ~ H1_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = hi1_foldbase_meta)
sum_h1_model2 <- summary(h1_model2) # +Age, sex (als factor), gruppe (als factor)

#qq plot
h1_multi_qq_plot <- ggplot(data = data.frame(resid = residuals(h1_model2)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "H1N1 multivar Q-Q Plot",
          x = NULL, 
          y = "Sample Quantiles")+
  theme(plot.title = element_text(hjust = 0.5))

# log first
h1_model2_lf <- lm(H1_fold_log2first ~ H1_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = hi1_foldbase_meta)
sum_h1_model2_lf <- summary(h1_model2_lf) # +Age, sex (als factor), gruppe (als factor)

#qq plot
h1_multi_lf_qq_plot <- ggplot(data = data.frame(resid = residuals(h1_model2_lf)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "H1N1 multivar Q-Q Plot / logF",
          x = NULL, 
          y = NULL)+
  theme(plot.title = element_text(hjust = 0.5))

## univeriables lineares model / H3N2
h3_model5 <- lm(H3_fold_log2 ~ H3_base_log2, data = hi3_foldbase)
sum_h3_model5 <- summary(h3_model5)

#qq plot
h3_uni_qq_plot <- ggplot(data = data.frame(resid = residuals(h3_model5)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "H3N2 univar Q-Q Plot",
          x = NULL, 
          y = "Sample Quantiles")+
  theme(plot.title = element_text(hjust = 0.5))

#Log first
h3_model5_lf <- lm(H3_fold_log2first ~ H3_base_log2, data = hi3_foldbase)
sum_h3_model5_lf <- summary(h3_model5_lf)

#qq plot
h3_uni_lf_qq_plot <- ggplot(data = data.frame(resid = residuals(h3_model5_lf)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "H3N2 univar Q-Q Plot / logF",
          x = NULL, 
          y = NULL)+
  theme(plot.title = element_text(hjust = 0.5))

## multivariables lineares model / H3N2
hi3_foldbase_meta <- hi3_foldbase %>% 
  left_join(basefile_withPID)

h3_model6 <- lm(H3_fold_log2 ~ H3_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = hi3_foldbase_meta)
sum_h3_model6 <- summary(h3_model6) # +Age, sex (als factor), gruppe (als factor)

#qq plot
h3_multi_qq_plot <- ggplot(data = data.frame(resid = residuals(h3_model6)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "H3N2 multi Q-Q Plot",
          x = NULL, 
          y = "Sample Quantiles")+
  theme(plot.title = element_text(hjust = 0.5))

# log first
h3_model6_lf <- lm(H3_fold_log2 ~ H3_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = hi3_foldbase_meta)
sum_h3_model6_lf <- summary(h3_model6_lf) # +Age, sex (als factor), gruppe (als factor)

#qq plot
h3_multi_lf_qq_plot <- ggplot(data = data.frame(resid = residuals(h3_model6_lf)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(   title = "H3N2 multi Q-Q Plot / logF",
          x = NULL, 
          y = NULL)+
  theme(plot.title = element_text(hjust = 0.5))

sumplot_qq_uni <- (h3_uni_qq_plot + h3_uni_lf_qq_plot) / (h1_uni_qq_plot + h1_uni_lf_qq_plot) / ((b_uni_qq_plot + b_uni_lf_qq_plot))
sumplot_qq_multi <- (h3_multi_qq_plot + h3_multi_lf_qq_plot) / (h1_multi_qq_plot + h1_multi_lf_qq_plot) / ((b_multi_qq_plot + b_multi_lf_qq_plot))


#### Checking if homoscedasticity is given (constant scattering of residuals independant of baseline)
###H3N1-Strain
# Remove NAs to ensure matching rows
h3_foldbase_clean <- hi3_foldbase %>% drop_na(H3_base_log2, H3_fold_log2)

# Fit linear model
h3_model5 <- lm(H3_fold_log2 ~ H3_base_log2, data = h3_foldbase_clean)
h3_model5_lf <- lm(H3_fold_log2first ~ H3_base_log2, data = h3_foldbase_clean)

## Check for homoscedasticity using residuals vs fitted values plot
h3_residuals_plot <- ggplot(data = h3_foldbase_clean, aes(x = fitted(h3_model5), y = resid(h3_model5))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = "Fitted Values", y = "Residuals",
       title = "H3N1 univariable model") +
  #ggtitle("H3 Residuals vs Fitted Values") +
  theme_minimal()

#log first
h3_residuals_plot_lf <- ggplot(data = h3_foldbase_clean, aes(x = fitted(h3_model5_lf), y = resid(h3_model5_lf))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = "Fitted Values", y = NULL,
       title = "H3N1 univariable model - LogF") +
  #ggtitle("H3 Residuals vs Fitted Values") +
  theme_minimal()

# Perform Breusch-Pagan test to check for homoscedasticity
durbinWatsonTest(h3_model5)
bp_test_h3 <- ncvTest(h3_model5)
print(bp_test_h3)

durbinWatsonTest(h3_model5_lf)
bp_test_h3_lf <- ncvTest(h3_model5_lf)
print(bp_test_h3_lf)

###H1N1-Strain
# Remove NAs to ensure matching rows
h1_foldbase_clean <- hi1_foldbase %>% drop_na(H1_base_log2, H1_fold_log2)

# Check for homoscedasticity using residuals vs fitted values plot
h1_residuals_plot <- ggplot(data = h1_foldbase_clean, aes(x = fitted(h1_model1), y = resid(h1_model1))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = NULL, y = "Residuals",
       title = "H1N1 univariable model") +
  #ggtitle("H1 Residuals vs Fitted Values") +
  theme_minimal()

#log first
h1_residuals_plot_lf <- ggplot(data = h1_foldbase_clean, aes(x = fitted(h1_model1_lf), y = resid(h1_model1_lf))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = NULL, y = NULL,
       title = "H1N1 univariable model - LogF") +
  #ggtitle("H1 Residuals vs Fitted Values") +
  theme_minimal()

# Perform Breusch-Pagan test to check for homoscedasticity
durbinWatsonTest(h1_model1)
bp_test_h1 <- ncvTest(h1_model1)
print(bp_test_h1)

durbinWatsonTest(h1_model1_lf)
bp_test_h1_lf <- ncvTest(h1_model1_lf)
print(bp_test_h1_lf)

###B-Strain
# Remove NAs to ensure matching rows
b_foldbase_clean <- b_foldbase %>% drop_na(Vic_base_log2, Vic_fold_log2)

# Check for homoscedasticity using residuals vs fitted values plot
b_residuals_plot <- ggplot(data = b_foldbase_clean, aes(x = fitted(b_model1), y = resid(b_model1))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = "Fitted Values", y = "Residuals",
       title = "B Victoria univariable model") +
  theme_minimal()

#Lof first
b_residuals_plot_lf <- ggplot(data = b_foldbase_clean, aes(x = fitted(b_model1_lf), y = resid(b_model1_lf))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(x = NULL, y = NULL,
       title = "B Victoria univariable model - LogF") +
  theme_minimal()

# Perform Breusch-Pagan test to check for homoscedasticity
dw_b <- durbinWatsonTest(b_model1)
bp_test_b <- ncvTest(b_model1)
print(bp_test_b)

durbinWatsonTest(b_model1)
bp_test_b_lf <- ncvTest(b_model1)
print(bp_test_b_lf)

# Combining RR-Plots
residplot_combined <- (b_residuals_plot + b_residuals_plot_lf) / (h1_residuals_plot + h1_residuals_plot_lf) / (h3_residuals_plot + h3_residuals_plot_lf)

###combining regression plots
#linreg_plot_test_combined = ((b_linreg_plot_baselog2 + b_residuals_plot + b_qq_plot) /
     #                 (hi1_linreg_plot_baselog2 + h1_residuals_plot + h1_qq_plot) / 
     #                 (hi3_linreg_plot_baselog2 + h3_residuals_plot + h3_qq_plot))



# Create a dataframe with Breusch-Pagan and Durbin-Watson test results
# Extract Breusch-Pagan p-values
bp_pvalues <- c(
  bp_test_b$p, 
  bp_test_b_lf$p, 
  bp_test_h1$p, 
  bp_test_h1_lf$p, 
  bp_test_h3$p, 
  bp_test_h3_lf$p
)

# Extract Durbin-Watson statistics
dw_stat <- c(as.numeric(durbinWatsonTest(b_model1)[2]),
             as.numeric(durbinWatsonTest(b_model1_lf)[2]),
             as.numeric(durbinWatsonTest(h1_model1)[2]),
             as.numeric(durbinWatsonTest(h1_model1_lf)[2]),
             as.numeric(durbinWatsonTest(h3_model5)[2]),
             as.numeric(durbinWatsonTest(h3_model5_lf)[2]))


# Create summary table
bp_dw_results <- data.frame(
  Strain = c("B", "B", "H1N1", "H1N1", "H3N2", "H3N2"),
  Model = c("Univariable", "Univariable log2 first", 
            "Univariable", "Univariable log2 first",
            "Univariable", "Univariable log2 first"),
  Breusch_Pagan_p_value = round(bp_pvalues, 4),
  Durbin_Watson_statistic = round(dw_stat, 3)
)

# View the table
print(bp_dw_results)


##### Compare univariant models (Adjustet R2, Estimates, P values)
#Put  models in a list:
models_uni_comb <- list(
  model_b = b_model1,
  b_logfirst = b_model1_lf,
  model_h1n1 = h1_model1,
  h1n1_logfirst = h1_model1_lf,
  model_h3n2 = h3_model5,
  h3n2_logfirst = h3_model5_lf
)
#Loop through models and extract info:
result_list <- lapply(names(models_uni_comb), function(name) {
  model <- models_uni_comb[[name]]
  sum_model <- summary(model)
  
  # Find the coefficient name containing "_base_log2"
  coef_name <- grep("_base_log2", rownames(sum_model$coefficients), value = TRUE)
  
  data.frame(
    Model = name,
    Base_Titer_Estimate = sum_model$coefficients[coef_name, "Estimate"],
    P_value = sum_model$coefficients[coef_name, "Pr(>|t|)"],
    Adjusted_R2 = sum_model$adj.r.squared,
    Multiple_R2 = sum_model$r.squared,
    AIC = AIC(model)
  )
})

model_results <- do.call(rbind, result_list)
model_results$P_value <- formatC(model_results$P_value, digits = 6, format = "f")

# Create the table (nice looking)
model_results_nice <- model_results %>%
  kable("html", digits = 4, align = rep('c', ncol(model_results))) %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(seq(2, nrow(model_results), by = 2), extra_css = "border-bottom: 2px solid black;")


## repeated measures model
####
###

