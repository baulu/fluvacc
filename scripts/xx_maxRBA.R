library(tibble)
library(tidyr)
library(dplyr)
library(here)

source(here::here("scripts","f01_maxRBA.R"))

titerset_long <- titerset_b %>% 
  mutate(strain = "B_victoria") %>% 
  rows_append(titerset_h1n1 %>% 
                mutate( strain = "H1N1")) %>% 
  rows_append(titerset_h3n2 %>% 
                mutate( strain = "H3N2")) %>% 
  select(SubjectID = PID, Pre = baseline_titers, Post = outcome_titers, Strain = strain)
  
titerset_format <- FormatTiters(titerset_long)

titerset_format_2 <- Calculate_maxRBA(titerset_format, method = "lm", )

# Assuming titerset_format_2$maxRBA is a named vector
long_df <- enframe(titerset_format_2$maxRBA, name = "PID", value = "maxRBA") %>% 
  mutate(PID = as.double(PID))

# View result
print(long_df)

titerset_comb_corrected_seroneg %>% 
  left_join(long_df) %>% view()




corrected_FC_matrix <- titerset_format_2$residualMatrix
corrected_FC_df <- corrected_FC_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "PID") %>% 
  mutate(PID = as.double(PID)) 


titerset5 <- titerset_comb_corrected_seroneg %>% 
  left_join(corrected_FC_df %>% select(PID, maxRBA_calcFC = B_victoria)) %>% 
  mutate(fold_change_log = log2(fold_change)) %>% 
  mutate(fold_change_corrbase_log = log2(fold_change_corrbase)) %>% 
  view()

