library(broom)
library(dplyr)
library(tidyr)
library(haven)
library(here)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(stringr)

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

# Pivot wider based on Titer_type (baseline/outcome). 
df_wide <- df_long %>%
  pivot_wider(
    names_from = Titer_type,
    values_from = value
  ) %>%
  rename(
    baseline = baseline_titer,
    outcome = outcome_titer
  )

# Calculate fold change and standardize. Note that standardisation is done within measurement, strain and metric
df_fc <- df_wide %>%
  mutate(fold_change = outcome / baseline) %>%
  group_by(measurement, strain, metric) %>%
  mutate(
    fc_median_g = median(fold_change, na.rm = TRUE),
    fc_sd_g     = sd(fold_change, na.rm = TRUE),
    base_median_g = median(baseline, na.rm = TRUE),
    base_sd_g     = sd(baseline, na.rm = TRUE),
    fc_std_norm     = (fold_change - fc_median_g) / fc_sd_g,
    baseline_std_norm = (baseline    - base_median_g) / base_sd_g
  ) %>%
  ungroup()

# !!! ## Mix of mean + SD a but unconventional, to discuss

to_exclude <- df_fc %>%
  filter(measurement == "hai") %>%
  filter(is.na(baseline) | is.na(outcome) | !is.finite(baseline) | !is.finite(outcome)) %>% 
  pull(Pat_ID) %>%  unique()
print(to_exclude)

df_fc <- df_fc %>% 
  filter(!Pat_ID %in% to_exclude)

df_responder <- df_fc %>%
  filter(measurement == "hai", metric == "titer") %>%
  group_by(Pat_ID) %>%
  mutate(
    fc_norm_max_ivt = max(fc_std_norm, na.rm = TRUE),
    d0_norm_paired = baseline_std_norm[which.max(fc_std_norm)]
  ) %>%
  slice_max(fc_std_norm, n = 1, with_ties = FALSE) %>%
  ungroup()

model <- lm(fc_norm_max_ivt ~ d0_norm_paired, data = df_responder)

df_responder <- df_responder %>%
  mutate(padjMFC = augment(model, data = .)$`.resid`)

discretize <- c(0.2, 0.3)
responseLabels <- c("lowResponder", "moderateResponder", "highResponder")

for (dis in discretize) {
  q_low <- quantile(df_responder$padjMFC, dis, na.rm = TRUE)
  q_high <- quantile(df_responder$padjMFC, 1 - dis, na.rm = TRUE)
  
  df_responder <- df_responder %>%
    mutate(!!paste0("padjMFC_d", dis * 100) := factor(case_when(
      is.na(padjMFC) ~ NA_character_,
      padjMFC <= q_low ~ responseLabels[1],
      padjMFC >= q_high ~ responseLabels[3],
      TRUE ~ responseLabels[2]
    ), levels = responseLabels))
}




## 1) compute group stats once and keep them
df_fc <- df_wide %>%
  mutate(fold_change = outcome / baseline) %>%
  group_by(measurement, strain, metric) %>%
  mutate(
    fc_median_g = median(fold_change, na.rm = TRUE),
    fc_sd_g     = sd(fold_change, na.rm = TRUE),
    base_median_g = median(baseline, na.rm = TRUE),
    base_sd_g     = sd(baseline, na.rm = TRUE),
    fc_std_norm     = (fold_change - fc_median_g) / fc_sd_g,
    baseline_std_norm = (baseline    - base_median_g) / base_sd_g
  ) %>%
  ungroup()

## ---- everything else as in your script up to df_responder/model ----

model <- lm(fc_norm_max_ivt ~ d0_norm_paired, data = df_responder)

df_responder <- df_responder %>%
  mutate(padjMFC = augment(model, data = .)$`.resid`)

## 2) De-standardize a standardized fold change (simple case)
## If df_responder still carries the group's fc_median_g/fc_sd_g (it should if you didn’t drop them),
## you can invert directly:
df_responder <- df_responder %>%
  mutate(
    fold_change_from_std = fc_norm_max_ivt * fc_sd_g + fc_median_g
  )

## 3) OPTIONAL: get an *adjusted* fold change (original units) after removing baseline effect
## Your residuals are on the standardized scale of y.
## If you want the value at a reference baseline of 0 (i.e., x_std = 0),
## first compute y_std adjusted, then invert:
b0 <- coef(model)[["(Intercept)"]]                 # α
b1 <- coef(model)[["d0_norm_paired"]]              # β

df_responder <- df_responder %>%
  mutate(
    # y_std_adj at x_std = 0 (keep population intercept, remove effect of d0_norm_paired)
    fc_std_adj0 = padjMFC + b0,                    # ε + α
    # back to original fold-change units using the group’s scale/center
    fold_change_adj0 = fc_std_adj0 * fc_sd_g + fc_median_g
  )

