library(broom)
library(dplyr)
library(tidyr)
library(haven)
library(here)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(stringr)

#
df <- read.xlsx(here::here("processed", "influenza_antibody_results_ngs.xlsx"))

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
  )
#compute FC
df_wide <- df_wide %>%
  mutate(fold_change = outcome / baseline)

## Standardization 
df_fc <- df_wide %>%
  group_by(measurement, strain, metric) %>%
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

to_exclude <- df_fc %>%
  filter(measurement == "hai") %>%
  filter(is.na(baseline) | is.na(outcome) | !is.finite(baseline) | !is.finite(outcome)) %>% 
  pull(Pat_ID) %>% unique()

df_fc <- df_fc %>% filter(!Pat_ID %in% to_exclude)

df_responder <- df_fc %>%
  filter(measurement == "hai", metric == "titer") %>%
  group_by(Pat_ID) %>%
  mutate(
    fc_norm_max_ivt = max(fc_std_norm, na.rm = TRUE),
    d0_norm_paired  = baseline_std_norm[which.max(fc_std_norm)]
  ) %>%
  slice_max(fc_std_norm, n = 1, with_ties = FALSE) %>%
  ungroup()

model <- lm(fc_norm_max_ivt ~ d0_norm_paired, data = df_responder)

df_responder <- df_responder %>%
  mutate(padjMFC = augment(model, data = .)$`.resid`)


# 1) De-standardize a standardized fold change (simple inversion)
df_responder <- df_responder %>%
  mutate(
    fold_change_from_std = fc_norm_max_ivt * fc_sd_g + fc_median_g
  )

# 2) Adjusted fold change in original units at reference baseline (x_std = 0)
b0 <- coef(model)[["(Intercept)"]]            # α
b1 <- coef(model)[["d0_norm_paired"]]         # β

df_responder <- df_responder %>%
  mutate(
    fc_std_adj0      = padjMFC + b0,          # ε + α
    fold_change_adj0 = fc_std_adj0 * fc_sd_g + fc_median_g
  )

model_max <- lm(fc_norm_max_ivt ~ d0_norm_paired, data = df_responder)
summary(model_max)

plot(model_max)


ggplot(data = df_responder, aes(x = d0_norm_paired, y = fc_norm_max_ivt)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  theme_minimal()

##### sum of responses

df_responder_sum <- df_fc %>%
  filter(measurement == "hai", metric == "titer") %>%
  group_by(Pat_ID) %>%
  summarise(
    pat_group = first(pat_group),
    fc_sum = sum(fold_change, na.rm = TRUE),
    fc_norm_sum_ivt = sum(fc_std_norm, na.rm = TRUE),
    d0_norm_sum = sum(baseline_std_norm, na.rm = TRUE),
    fc_median_g   = median(fold_change, na.rm = TRUE),
    fc_sd_g       = sd(fold_change, na.rm = TRUE),
    base_median_g = median(baseline, na.rm = TRUE),
    base_sd_g     = sd(baseline, na.rm = TRUE),
    .groups = "drop"
  )

model_sum <- lm(fc_norm_sum_ivt ~ d0_norm_sum, data = df_responder_sum)
summary(model_sum)
plot(model_sum)

ggplot(data = df_responder_sum, aes(x = d0_norm_sum, y = fc_norm_sum_ivt)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  theme_minimal()

df_responder_sum <- df_responder_sum %>%
  mutate(padjMFC = augment(model_sum, data = .)$`.resid`)

discretize <- c(0.2, 0.3)
responseLabels <- c("lowResponder", "moderateResponder", "highResponder")

for (dis in discretize) {
  q_low <- quantile(df_responder_sum$padjMFC, dis, na.rm = TRUE)
  q_high <- quantile(df_responder_sum$padjMFC, 1 - dis, na.rm = TRUE)
  
  df_responder_sum <- df_responder_sum %>%
    mutate(!!paste0("padjMFC_d", dis * 100) := factor(case_when(
      is.na(padjMFC) ~ NA_character_,
      padjMFC <= q_low ~ responseLabels[1],
      padjMFC >= q_high ~ responseLabels[3],
      TRUE ~ responseLabels[2]
    ), levels = responseLabels))
}


### look at single fold_changes and what happens with them when correcting by linear regression
#compute FC

## Standardization 
df_fc <- df_wide %>%
  group_by(measurement, strain, metric) %>%
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

to_exclude <- df_fc %>%
  filter(measurement == "hai") %>%
  filter(is.na(baseline) | is.na(outcome) | !is.finite(baseline) | !is.finite(outcome)) %>% 
  pull(Pat_ID) %>% unique()

df_fc <- df_fc %>% filter(!Pat_ID %in% to_exclude)

model_fcsd <- lm(fc_std_norm ~ baseline_std_norm, data = df_fc)
summary(model_fcsd)

model_logfc <- lm(log(fold_change) ~ log(baseline), data = df_fc)
summary(model_logfc)
plot(model_logfc)

df_fc <- df_fc %>%
  mutate(padjMFC = augment(model_logfc, data = .)$`.resid`)

# 2) Adjusted fold change in original units at reference baseline (x_std = 0)
b0 <- coef(model_logfc)[["(Intercept)"]]            # α
b1 <- coef(model_logfc)[["log(baseline)"]]         # β

df_fc <- df_fc %>%
  mutate(
    fc_std_adj0      = padjMFC + b0,          # ε + α
    fold_change_adj0 = fc_std_adj0 * fc_sd_g + fc_median_g
  )

ggplot(data = (df_fc %>% filter(strain == "H3N2")), aes(x = log(baseline), y = log(fold_change))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  theme_minimal()


# plotting raw fc / baseline-titer
ggplot(data = df_wide,aes(
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

model_loglog <- lm(log(fold_change) ~log(baseline), data = df_wide)
summary(model_loglog)

#looking only at values >4

ggplot(data = df_wide %>% filter(fold_change >= 3.9),aes(
  x = log(baseline),
  y = log(fold_change)
)) +
  geom_point() +
  theme_minimal()
## looking pretty good for linear model

df_hai_ov4 <- df_wide %>% 
  filter(fold_change >= 3.9) %>% 
  filter(measurement == "hai")

model_haiov4 <- lm(log(fold_change) ~ log(baseline), data = df_hai_ov4)
summary(model_haiov4)

ggplot(data = df_responder, aes(x = d0_norm_paired, y = fc_norm_max_ivt)) +
  geom_point() +
  #geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, color = "red") +
  theme_minimal()

model <- lm(fc_norm_max_ivt ~ d0_norm_paired, data = df_responder)
model_sq <- lm(fc_norm_max_ivt ~ d0_norm_paired + I(d0_norm_paired^2), data = df_responder)

plot(model)

summary(model)
summary(model_sq)

ggplot(data = df_responder, aes(x = d0_norm_paired, y = log(fc_norm_max_ivt))) +
  geom_point() #+
  #geom_smooth(method = "lm", formula = log(y) ~ log(x), se = TRUE, color = "red") +
  theme_minimal()

  
df_responder_sum <- df_fc %>%
    filter(measurement == "hai", metric == "titer") %>%
    group_by(Pat_ID) %>%
    summarise(
      pat_group = first(pat_group),
      fc_norm_sum_ivt = sum(fc_std_norm, na.rm = TRUE),
      d0_norm_sum = sum(baseline_std_norm, na.rm = TRUE),
      .groups = "drop"
    )
  
model_sum <- lm(fc_norm_sum_ivt ~ d0_norm_sum, data = df_responder_sum)
summary(model_sum)
plot(model_sum)

ggplot(data = df_responder_sum, aes(x = d0_norm_sum, y = fc_norm_sum_ivt)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  theme_minimal()  

  df_responder_sum <- df_responder_sum %>%
    mutate(padjMFC = augment(model_sum, data = .)$`.resid`)
  
  discretize <- c(0.2, 0.3)
  responseLabels <- c("lowResponder", "moderateResponder", "highResponder")
  
  for (dis in discretize) {
    q_low <- quantile(df_responder_sum$padjMFC, dis, na.rm = TRUE)
    q_high <- quantile(df_responder_sum$padjMFC, 1 - dis, na.rm = TRUE)
    
    df_responder_sum <- df_responder_sum %>%
      mutate(!!paste0("padjMFC_d", dis * 100) := factor(case_when(
        is.na(padjMFC) ~ NA_character_,
        padjMFC <= q_low ~ responseLabels[1],
        padjMFC >= q_high ~ responseLabels[3],
        TRUE ~ responseLabels[2]
      ), levels = responseLabels))
  }
  
  
dif_fc <- df_responder %>% select(Pat_ID, baseline, outcome, fold_change, fc_norm_max_ivt, fold_change_adj0) %>% 
  print()

dif_fc %>%
  select(Pat_ID, baseline, fold_change, fold_change_adj0) %>%
  pivot_longer(
    cols = c(fold_change, fold_change_adj0),
    names_to = "fc_type",
    values_to = "FC_value"
  ) %>%
  ggplot(aes(x = baseline, y = FC_value, color = fc_type, group = Pat_ID)) +
  geom_jitter() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 4, color = "green", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Raw vs Adjusted Fold Change per Patient",
    x = "Baseline",
    y = "Fold change value"
  )



  
  