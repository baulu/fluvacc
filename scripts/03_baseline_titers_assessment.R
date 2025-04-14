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

#load data
microneut_analysis_raw  <- read.xlsx("processed/microneut_analysis_raw.xlsx")
hai_analysis_raw  <- read.xlsx("processed/hai_analysis_raw.xlsx")
influenza_antibody_results <- read.xlsx("processed/influenza_antibody_results.xlsx")
basefile_withPID <- read.xlsx("processed/FluVac_basefile_withPID.xlsx")

# 2) -------------------------------------------------------------------------------------------------------------------
##### assess correlation of microneutralisation fold changes with baseline titer
#put together dataset with fold changes and baseline titer for all strains
#h1
h1_base <-  microneut_analysis_raw %>%
    filter(Sampling_number == 1) %>% 
    select(PID, H1_baseline = FluA_H1_ic50)

h1_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, H1_fold = FluV_H1_fold, H1_fold_log2first = FluV_H1_fold_log2first)
  
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

hi1_linreg_plot_log2first <- hi1_foldbase %>% 
  ggplot(aes(x = H1_base_log2, y = H1_fold_log2first)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  #scale_x_continuous(limits = c(0, 13000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -0.25) + # Display Pearson's r
  labs(y = NULL, x = NULL) 

# Perform Pearson correlation tests
cor_result1 <- cor.test(hi1_foldbase$H1_baseline, hi1_foldbase$H1_fold_log2, method = "pearson")
cor_result2 <- cor.test(hi1_foldbase$H1_base_log2, hi1_foldbase$H1_fold_log2, method = "pearson")
cor_result3_logfirst <- cor.test(hi1_foldbase$H1_base_log2, hi1_foldbase$H1_fold_log2first, method = "pearson")

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

cor_comb1 <- cor_df %>% 
  rows_append(cor_df2) %>% 
  rows_append(cor_df3_log2first) %>% 
  print()
 

###h3
h3_base <-  microneut_analysis_raw %>%
  filter(Sampling_number == 1) %>% 
  select(PID, H3_baseline = FluA_H3_ic50)

h3_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, H3_fold = FluV_H3_fold, H3_fold_log2first = FluV_H3_fold_log2first)

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

hi3_linreg_plot_log2first <- hi3_foldbase %>% 
  ggplot(aes(x = H3_base_log2, y = H3_fold_log2first)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +  # Add regression line with confidence interval
  #scale_x_continuous(limits = c(0, 13000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -0.25) + # Display Pearson's r
  labs(y = NULL, x = NULL) 

# Perform Pearson correlation tests
corh3_result1 <- cor.test(hi3_foldbase$H3_baseline, hi3_foldbase$H3_fold_log2, method = "pearson")
corh3_result2 <- cor.test(hi3_foldbase$H3_base_log2, hi3_foldbase$H3_fold_log2, method = "pearson")
corh3_result3_log2first <- cor.test(hi3_foldbase$H3_base_log2, hi3_foldbase$H3_fold_log2first, method = "pearson")


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

corh3_comb <- corh3_df %>% 
  rows_append(corh3_df2) %>% 
  rows_append(corh3_df2_log2first) %>% 
  print()

###B
b_base <-  microneut_analysis_raw %>%
  filter(Sampling_number == 1) %>% 
  select(PID, Vic_baseline = FluA_Vic_ic50)

b_fold <-  microneut_analysis_raw %>%
  filter(Sampling_number == 2) %>% 
  select(PID, Vic_fold = FluV_Vic_fold, Vic_fold_log2first = FluV_Vic_fold_log2first)

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
  #stat_cor(method = "pearson", label.x = 6, label.y = -4) + # Display Pearson's r
  labs(y = "B - log2 Fold-Change", x = NULL) 

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
  #scale_x_continuous(limits = c(0, 13000)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add horizontal line
  theme_minimal()+
  stat_cor(method = "pearson", label.x = 6, label.y = -0.25) + # Display Pearson's r
  labs(y = NULL, x = NULL) 

# Perform Pearson correlation tests
corb_result1 <- cor.test(b_foldbase$Vic_baseline, b_foldbase$Vic_fold_log2, method = "pearson")
corb_result2 <- cor.test(b_foldbase$Vic_base_log2, b_foldbase$Vic_fold_log2, method = "pearson")
corb_result3_log2first <- cor.test(b_foldbase$Vic_base_log2, b_foldbase$Vic_fold_log2first, method = "pearson")


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

corb_comb <- corb_df %>% 
  rows_append(corb_df2) %>%
  rows_append(corb_df3_log2first) %>% 
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
## univeriables lineares model / B-Victoria
model1 <- lm(Vic_fold_log2 ~ Vic_base_log2, data = b_foldbase)
summary(model1)
qqnorm(residuals(model1))
qqline(residuals(model1), col = "red")

model1_lf <- lm(Vic_fold_log2first ~ Vic_base_log2, data = b_foldbase)
summary(model1_lf )
qqnorm(residuals(model1_lf ))
qqline(residuals(model1_lf ), col = "red")

## multivariables lineares model / B-Victoria
b_foldbase_meta <- b_foldbase %>% 
  left_join(basefile_withPID)

model2 <- lm(Vic_fold_log2 ~ Vic_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = b_foldbase_meta)
summary(model2) # +Age, sex (als factor), gruppe (als factor)
qqnorm(residuals(model2))
qqline(residuals(model2), col = "red")

model2_lf <- lm(Vic_fold_log2first ~ Vic_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = b_foldbase_meta)
summary(model2_lf) # +Age, sex (als factor), gruppe (als factor)
qqnorm(residuals(model2_lf))
qqline(residuals(model2_lf), col = "red")

## univeriables lineares model / H1N1
model3 <- lm(H1_fold_log2 ~ H1_base_log2, data = hi1_foldbase)
summary(model3)
qqnorm(residuals(model3))
qqline(residuals(model3), col = "red")

## multivariables lineares model / H1N1
hi1_foldbase_meta <- hi1_foldbase %>% 
  left_join(basefile_withPID)

model4 <- lm(H1_fold_log2 ~ H1_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = hi1_foldbase_meta)
summary(model4) # +Age, sex (als factor), gruppe (als factor)
qqnorm(residuals(model4))
qqline(residuals(model4), col = "red")

model4_lf <- lm(H1_fold_log2first ~ H1_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = hi1_foldbase_meta)
summary(model4_lf) # +Age, sex (als factor), gruppe (als factor)
qqnorm(residuals(model4_lf))
qqline(residuals(model4_lf), col = "red")

## univeriables lineares model / H3N2
model5 <- lm(H3_fold_log2 ~ H3_base_log2, data = hi3_foldbase)
summary(model5)
qqnorm(residuals(model5))
qqline(residuals(model5), col = "red")

## multivariables lineares model / H3N2
hi3_foldbase_meta <- hi3_foldbase %>% 
  left_join(basefile_withPID)

model6 <- lm(H3_fold_log2 ~ H3_base_log2 + gen_age + as.factor(gen_sex) + as.factor(study_group), data = hi3_foldbase_meta)
summary(model6) # +Age, sex (als factor), gruppe (als factor)
qqnorm(residuals(model6))
qqline(residuals(model6), col = "red")

#Vergleiche --> adjustet R2, Residuals, AIC

## repeated measures model
####
###

# residuals plot, qq plot

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

