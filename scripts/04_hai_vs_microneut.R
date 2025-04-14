#Index
# 1) Load libraries and data 
# 2) Correlation of HAI and microneutralisation assays

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
#### Comparing HAI and Microneutralisation assay
#preparation
microneut_analysis_raw_p <- microneut_analysis_raw %>% #CAVE: USED IN OTHER SCRIPTS - e.g. 02_hai_micro...
  select(pat_group, Sampling_number, H1 = FluA_H1_ic50, H3 = FluA_H3_ic50, Vic = FluA_Vic_ic50, PID) %>% 
  pivot_longer(cols = 3:5, names_to = "strain", values_to = "ic50")  

#preparation
hai_analysis_raw_p <- hai_analysis_raw %>% #CAVE: USED IN OTHER SCRIPTS - e.g. 02_hai_micro...
  select(pat_group, Sampling_number, H1N1, H3N2, BVic, BYam, PID) %>% 
  pivot_longer(cols = 3:6, names_to = "strain", values_to = "result")  


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

