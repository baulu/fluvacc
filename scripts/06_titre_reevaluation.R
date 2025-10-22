## Setup

###Load Libraries
library(broom)
library(dplyr)
library(tidyr)
library(haven)
library(here)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(readr)
library(stringr)
###Load Dataset
df <- read.xlsx(here::here("processed", "influenza_antibody_results_ngs.xlsx"))

### Reshape: long to wide and exclude microneutralisation assay
# Pivot longer to isolate strain + metric info
df_long <- df %>%
  pivot_longer(
    cols = starts_with("hai_") | starts_with("microneut_") | starts_with("mirconeut_"),
    names_to = "key",
    values_to = "value"
  ) %>% 
  mutate(
    measurement = case_when(
      str_starts(key, "hai_") ~ "hai",
      str_starts(key, "microneut_") ~ "microneut",
      str_starts(key, "mirconeut_") ~ "microneut"
    ),
    strain = str_extract(key, "(H1N1|H3N2|BVic|BYam)"),
    metric = case_when(
      str_detect(key, "ic50") ~ "ic50",
      TRUE ~ "titer"
    )
  ) %>%
  select(Pat_ID, pat_group, Titer_type, strain, measurement, metric, value)

# Pivot wider based on Titer_type (baseline/outcome)
df_wide <- df_long %>%
  pivot_wider(
    names_from = Titer_type,
    values_from = value
  ) %>%
  rename(
    baseline = baseline_titer,
    outcome  = outcome_titer
  ) %>% 
  filter(measurement == "hai")

### Compute fold change
df_wide <- df_wide %>%
  mutate(fold_change = outcome / baseline)

###drop NA
df_wide <- df_wide %>%
  drop_na()


## Unadjusted HAI Fold-Changes / baseline values
###basic plot

ggplot(df_wide, aes(x = baseline, y = fold_change, color = fold_change > 3.9)) +
  geom_jitter(
    width = 10,    # increase for more horizontal jitter
    height = 0,  # increase for more vertical jitter
    seed = 5        # optional: makes jitter reproducible
  ) +
  scale_color_manual(values = c(`TRUE` = "green", `FALSE` = "red"))+
  theme_minimal()

#--> Most unadjusted responders with low baseline-titers. Need for and impact of adjustment for baseline titers likely not big.

### Add linear regression line 
ggplot(df_wide, aes(x = baseline, y = fold_change, color = fold_change > 3.9)) +
  geom_jitter(
    width = 10,    # increase for more horizontal jitter
    height = 0,  # increase for more vertical jitter
    seed = 5        # optional: makes jitter reproducible
  ) +
  scale_color_manual(values = c(`TRUE` = "green", `FALSE` = "red"))+
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red")+
  theme_minimal()

#--> linearity assumptions probably does not hold and using this model for Fold change correction seems problematic.


model_base <- lm(fold_change~baseline, data = df_wide)
summary(model_base) %>% print()


### plot log(FC)/log(baseline) + corresponding liner model log(FC) ~ log(baseline titers)
ggplot(data = df_wide %>% filter(strain == "BYam"), aes(
  x = log(baseline),
  y = log(fold_change),
  color = fold_change > 3.9   # TRUE if >4, FALSE otherwise
)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red")+
  scale_color_manual(
    values = c("TRUE" = "green", "FALSE" = "red")
  )+
  theme_minimal()
# looking somewhat more linear


model_loglog <- lm(log(fold_change)~log(baseline), data = df_wide)
summary(model_loglog) %>% print()
#checking model for linearity assumptions
plot(model_loglog)


##### Model max-titer as by marco
# Calculate fold change and standardize. Note that standardisation is done within measurement, strain and metric
df_fc <- df_wide %>%
  group_by(strain) %>%
  mutate(
    fc_median_g   = median(fold_change, na.rm = TRUE),
    fc_sd_g       = sd(fold_change, na.rm = TRUE),
    base_median_g = median(baseline, na.rm = TRUE),
    base_sd_g     = sd(baseline, na.rm = TRUE),
    # standardized (hybrid: median + SD)
    fc_std_norm       = (fold_change - fc_median_g) / fc_sd_g,
    baseline_std_norm = (baseline    - base_median_g) / base_sd_g
  ) %>%
  ungroup()

model_max <- lm(fc_std_norm ~ baseline_std_norm, data = df_fc)

df_fc_std <- df_fc %>%
  mutate(padjMFCstd = augment(model_max, data = .)$`.resid`)

discretize <- c(0.2, 0.3)
responseLabels <- c("lowResponder", "moderateResponder", "highResponder")

for (dis in discretize) {
  q_low <- quantile(df_fc_std$padjMFCstd, dis, na.rm = TRUE)
  q_high <- quantile(df_fc_std$padjMFCstd, 1 - dis, na.rm = TRUE)
  
  df_fc_std <- df_fc_std %>%
    mutate(!!paste0("padjMFC_d", dis * 100) := factor(case_when(
      is.na(padjMFCstd) ~ NA_character_,
      padjMFCstd <= q_low ~ responseLabels[1],
      padjMFCstd >= q_high ~ responseLabels[3],
      TRUE ~ responseLabels[2]
    ), levels = responseLabels))
}

#Backtransform Fold changes  (to foldchange of 1 as reference! - not as done by marco!)
# 1) Adjusted fold change in original units at reference baseline (x_std = 0)
b0_max <- coef(model_max)[["(Intercept)"]]            # α
b1_max <- coef(model_max)[["baseline_std_norm"]]         # β

df_fc_std <- df_fc_std %>%
  mutate(
    fc_std_adj0      = padjMFCstd + 0,          # ε + α
    fc_adj_backt = fc_std_adj0 * fc_sd_g + fc_median_g
  ) %>% 
  select(Pat_ID, pat_group, strain, baseline, outcome, fold_change,fc_std_adj0, fc_adj_backt, responsed20_std = padjMFC_d20, responsed30_std = padjMFC_d30)



### Show in plot
df_fc_std %>% 
  select(baseline, fold_change, fc_adj_backt) %>%
  pivot_longer(cols = 2:3, names_to = "fc_type", values_to = "fc_value") %>%
  ggplot()+
  geom_jitter(aes(x = baseline, y = fc_value, color = fc_type))+
  theme_minimal()

###### use a loglog transformed model instead (OVER ALL STRAINS!!! needs to be corrected!)
df_responder_loglog <- df_fc %>%
  filter(measurement == "hai", metric == "titer") %>%
  mutate(
    log_fc = log(fold_change),
    log_baseline = log(baseline),
  )

model_max_loglog <- lm(log_fc ~ log_baseline, data = df_responder_loglog)

df_responder_loglog <- df_responder_loglog %>%
  mutate(padjMFClog = augment(model_max_loglog, data = .)$`.resid`)

discretize <- c(0.2, 0.3)
responseLabels <- c("lowResponder", "moderateResponder", "highResponder")

for (dis in discretize) {
  q_low <- quantile(df_responder_loglog$padjMFClog, dis, na.rm = TRUE)
  q_high <- quantile(df_responder_loglog$padjMFClog, 1 - dis, na.rm = TRUE)
  
  df_responder_loglog <- df_responder_loglog %>%
    mutate(!!paste0("padjMFC_d", dis * 100) := factor(case_when(
      is.na(padjMFClog) ~ NA_character_,
      padjMFClog <= q_low ~ responseLabels[1],
      padjMFClog >= q_high ~ responseLabels[3],
      TRUE ~ responseLabels[2]
    ), levels = responseLabels))
}


#Backtransform Fold changes (to foldchange of 1 as reference! - not as done by marco!)
# 1) Adjusted fold change in original units at reference baseline (x_std = 0)
b0_loglog <- coef(model_max_loglog)[["(Intercept)"]]            # α
b1_loglog <- coef(model_max_loglog)[["log_baseline"]]         # β

df_responder_loglog <- df_responder_loglog %>%
  mutate(
    fc_loglog_adj0      = padjMFClog + 0,          # ε + α
    fc_loglog_backt = exp(fc_loglog_adj0)
  ) %>% 
  select(Pat_ID, pat_group, strain, baseline, outcome, fold_change,fc_loglog_adj0, fc_loglog_backt, responsed20_log = padjMFC_d20, responsed30_log = padjMFC_d30)


### Show in plot 
df_fc_std %>% 
  select(baseline, fold_change, fc_adj_backt) %>%
  pivot_longer(cols = 2:3, names_to = "fc_type", values_to = "fc_value") %>%
  ggplot()+
  geom_jitter(aes(x = baseline, y = fc_value, color = fc_type))+
  theme_minimal()

df_responder_loglog %>% 
  select(Pat_ID, strain, baseline, fold_change, fc_loglog_backt) %>% 
  left_join(
    df_fc_std %>% select(Pat_ID, strain, fc_adj_backt),  # keep join keys
    by = join_by(Pat_ID, strain)
  ) %>%
  pivot_longer(cols = 4:6, names_to = "fc_type", values_to = "fc_value") %>% 
  ggplot()+
  geom_jitter(aes(x = baseline, y = fc_value, color = fc_type),
              width = 5,    # increase for more horizontal jitter
              height = 0.05)+  # increase for more vertical jitter
  theme_minimal()

df_responder_loglog %>% 
  select(Pat_ID, strain, baseline, fold_change, fc_loglog_backt) %>% 
  left_join(
    df_fc_std %>% select(Pat_ID, strain, fc_adj_backt),  # keep join keys
    by = join_by(Pat_ID, strain)
  ) %>%
  pivot_longer(cols = 4:6, names_to = "fc_type", values_to = "fc_value") %>% 
ggplot(aes(x = baseline, y = fc_value, color = fc_value > 3.9)) +
  geom_jitter(width = 5, height = 0, seed = 1) +
  facet_wrap(~ fc_type, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Baseline", y = "Fold change")

# compare d30 categorisation betwen std and loglog models
df_fc_std %>% 
  left_join(df_responder_loglog) %>% 
  #pivot_longer(cols = contains("responsed30"), names_to = "response_type", values_to = "response_d30") %>% 
  print()
  
df_fc_std %>%
  left_join(df_responder_loglog) %>% 
  select(baseline, fold_change,fc_loglog_backt, fc_adj_backt, responsed30_std, responsed30_log) %>%
  #pivot_longer(cols = 2:3, names_to = "fc_type", values_to = "fc_value") %>%
  ggplot()+
  geom_jitter(aes(x = baseline, y = fc_loglog_backt, color = responsed30_log))+
  theme_minimal()

df_fc_std %>%
  left_join(df_responder_loglog) %>% 
  select(baseline, fold_change, fc_loglog_backt, fc_adj_backt, responsed30_std, responsed30_log) %>%
  #pivot_longer(cols = 2:3, names_to = "fc_type", values_to = "fc_value") %>%
  ggplot()+
  geom_jitter(aes(x = baseline, y = fc_adj_backt, color = responsed30_std))+
  theme_minimal()


df_wide %>% 
  filter(baseline <=250) %>% 
  select(Pat_ID, strain, baseline) %>% 
  pivot_wider(names_from = strain, values_from = baseline) %>% 
  drop_na() %>% 
  print()



############## 2. Updated Models (GLMM / GEE / Beyer) Analysis of raw data

influenza_folds2 <- read_csv(here::here("data/titer_classification_summary_2ndround.csv"))

inlfzenza_raw_updatedmodels_rank <- df_wide %>% 
  left_join(influenza_folds2 %>% select(Pat_ID, beyer_sum_rank, beyer_sum_d20, beyer_sum_d30, GLMM_rank)) %>% 
  print()

inlfzenza_raw_updatedmodels <- inlfzenza_raw_updatedmodels_rank %>%   
  select(Pat_ID, pat_group, strain, baseline, fold_change,beyer_sum_d20, beyer_sum_d30) %>% 
  print()

inlfzenza_raw_updatedmodels_rank %>% distinct(Pat_ID, .keep_all = TRUE) %>% print()

#look at droupouts by high baseline titer
inlfzenza_raw_updatedmodels %>% filter(baseline >=250) %>% left_join(inlfzenza_raw_updatedmodels) %>% 
  select(Pat_ID, pat_group, strain) %>% 
  print()



# Assuming you already have inlfzenza_raw_updatedmodels from your previous code
# Define color bins for 'baseline'
inlfzenza_raw_updatedmodels <- inlfzenza_raw_updatedmodels %>%
  mutate(
    baseline_group = case_when(
      baseline >= 0 & baseline <= 9 ~ "0-9",
      baseline >= 10 & baseline <= 39 ~ "10-39",
      baseline >= 40 ~ "40+"
    )
  )

# Create ggplot
ggplot(inlfzenza_raw_updatedmodels %>% filter(beyer_sum_d30 == "highResponder"), aes(x = strain, y = fold_change, color = baseline_group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 4, color = "lightgreen") +
  facet_wrap(~ Pat_ID) +
  scale_color_manual(
    values = c("0-9" = "grey", "10-39" = "orange", "40+" = "green"),
    name = "Baseline Titer"
  ) +
  scale_y_continuous(limits = c(0, 40))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "d30 High Responder (Beyer)",
    x = "Strain",
    y = "Fold Change"
  )


