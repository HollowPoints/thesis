# Install packages
install.packages(c("anthro", "zscorer","haven","anthroplus", "tidyr", "dplyr", "gt"))

# Library packages
library(anthro)
library(zscorer)
library(haven)
library(anthroplus)
library(dplyr)
library(tidyr)
library(stringr)
library(gt)

# call functions script
source("functions.R")




# commands to push rscript changes to git hub
# git init
# git add MScTH_rscript.R
# git push origin main


# TO DO

# gt tables package for descriptive statistics  rmarkdown, quarto(how to add the tables to quarto)

# how mamy variables are missing from antrho datasets DONE!
# kind of rule of days difference in antrho to serology measurements DONE!



# write.csv(complete_data, "C:/Users/grigoris.kalampoukas/MScThesis/complete_data.csv", row.names = FALSE)




# data loading

# anthropometric data
obesity_BiB <- read_dta("C:/Users/grigoris.kalampoukas/MScThesis/Grigoris/marato_obesity_BiB.dta")
obesity_RHEA <-  read_dta("C:/Users/grigoris.kalampoukas/MScThesis/Grigoris/marato_obesity_RHEA.dta")
obesity_INMA <-  read_dta("C:/Users/grigoris.kalampoukas/MScThesis/Grigoris/marato_obesity_INMA.dta")

# serological data
serology_BiB <- read_dta("C:/Users/grigoris.kalampoukas/MScThesis/Grigoris/marato_serology_Bib.dta")
serology_RHEA <- read_dta("C:/Users/grigoris.kalampoukas/MScThesis/Grigoris/marato_serology_Rhea.dta")
serology_INMA <-  read_dta("C:/Users/grigoris.kalampoukas/MScThesis/Grigoris/marato_serology_INMA.dta")
#backups 

bckp_obesity_BiB <- obesity_BiB
bckp_obesity_INMA <- obesity_INMA
bckp_obesity_RHEA <- obesity_RHEA

bckp_serology_BiB <- serology_BiB
bckp_serology_INMA <- serology_INMA
bckp_serology_RHEA <- serology_RHEA


# transfor RHEA to long format

# find RHEA measurement columns.

measurement_cols <- names(obesity_RHEA) %>% 
  str_subset("\\d+$")


# extract variable roots   (e.g., "height", "height_age", "weight", ...)
variable_roots <- measurement_cols %>%
  str_extract("^[^\\d]+") %>%
  unique()


pattern_regex <- paste0("^(", paste(variable_roots, collapse = "|"), ")(\\d+)$")

# long format database

obesity_RHEA <- obesity_RHEA %>%
  pivot_longer(
    cols = all_of(measurement_cols),
    names_to = c(".value", "timepoint"),
    names_pattern = pattern_regex
  ) %>%
  # Drop rows where all values are NA
  filter(if_any(all_of(variable_roots), ~ !is.na(.)))




# make all columns of same variables have the same name
obesity_INMA <- obesity_INMA %>%
  rename(
    agecd_cgrowth = height_age,
    cheight = height_,
    agecm_cgrowth = age_months,
    agecy_cgrowth = age_years,
    cweight = weight_,
    cabdo = waistcirc_,
    csubscap = subscapsf_,
    ctriceps = tricepsf_
  )




obesity_RHEA <- obesity_RHEA %>%
  rename(
    agecd_cgrowth = weight_age,
    cheight = height_,
    cweight = weight_,
    cabdo = waistcirc_,
    csubscap = subscapsf_,
    ctriceps = tricepsf_
  )

# calculate RHEA measurment age in months and years
obesity_RHEA$agecm_cgrowth <- round(obesity_RHEA$agecd_cgrowth / 30.4375, 0)
obesity_RHEA$agecy_cgrowth <- round(obesity_RHEA$agecm_cgrowth /12, 0)
obesity_RHEA$cbmi <- obesity_RHEA$cweight / ((obesity_RHEA$cheight / 100) ^ 2 )





# initialize z score values
obesity_BiB[, c("bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore")] <- NA
obesity_INMA[, c("bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore")] <- NA
obesity_RHEA[, c("bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore")] <- NA






# check for NAs in agecd_growth(days) if yes calculate it with month *  30.4375 for INMA and RHEA
obesity_INMA$agecd_cgrowth[is.na(obesity_INMA$agecd_cgrowth)] <- round(obesity_INMA$agecm_cgrowth[is.na(obesity_INMA$agecd_cgrowth)] * 30.4375, 0)


# bmi calculation for INMA cohort
obesity_INMA$cbmi <- obesity_INMA$cweight / ((obesity_INMA$cheight / 100) ^ 2)

# calculate generate z-scores for weight, height, BMI, triceps skinfold and subscapular skinfold according to WHO 0-5 years old


obesity_BiB <- assign_anthro_scores(obesity_BiB)
obesity_INMA <- assign_anthro_scores(obesity_INMA)
obesity_RHEA <- assign_anthro_scores(obesity_RHEA)







# calculate generate z-scores for weight, height, BMI according to WHO 5-19 years old


obesity_BiB <- assign_anthroplus_scores(obesity_BiB)
obesity_INMA <- assign_anthroplus_scores(obesity_INMA)
obesity_RHEA <- assign_anthroplus_scores(obesity_RHEA)




# function that finds and removes implausible values (±5 SDs) removes the whole row if all z scores are implausible values. 
# Creates a exclusion log for each dataset with h_id and values excluded.





clean_zscores(obesity_BiB)
clean_zscores(obesity_INMA)
clean_zscores(obesity_RHEA)



# function to classify children as Normal, Overweight, Obese based on z score , depending on age < 5 and age > 5


obesity_BiB <- classify_bmi_category(obesity_BiB)
obesity_INMA <- classify_bmi_category(obesity_INMA)
obesity_RHEA <- classify_bmi_category(obesity_RHEA)


# indetify central obesity for children over 4 years of age. 



obesity_BiB <- classify_central_obesity(obesity_BiB)
obesity_INMA <- classify_central_obesity(obesity_INMA)
obesity_RHEA <- classify_central_obesity(obesity_RHEA)


# function to classify stunted using height for age z-score


obesity_BiB <- classify_stunting(obesity_BiB)
obesity_INMA <- classify_stunting(obesity_INMA)
obesity_RHEA <- classify_stunting(obesity_RHEA)

# function to make serology to long format

serology_BiB <- longify_serology_0_4(serology_BiB)
serology_INMA <- longify_serology_0_4(serology_INMA)
serology_RHEA <- longify_serology_0_4(serology_RHEA)


# calculate age in days for all cohorts, and place the new column next to measurment date

serology_BiB <- serology_BiB %>%
  mutate(age_days = age * 365) %>%
  relocate(age_days, .after = 10)

serology_INMA <- serology_INMA %>%
  mutate(age_days = age * 365) %>%
  relocate(age_days, .after = 10)

serology_RHEA <- serology_RHEA %>%
  mutate(age_days = age * 365) %>%
  relocate(age_days, .after = 10)



# function to match obesity data to serology data
# for each serology timepoint measurment the functions finds the best match
# based on BMI existence, then other antrho variables count("cabdo", "csubscap", "ctriceps") 
# and finally by measuremnt age difference with a cuttoff of 60 days.





data_BiB <- match_growth_to_serology(serology_BiB, obesity_BiB, cutoff_days = 60)
data_INMA <- match_growth_to_serology(serology_INMA, obesity_INMA, cutoff_days = 60)
data_RHEA <- match_growth_to_serology(serology_RHEA, obesity_RHEA, cutoff_days = 60)

complete_data <- bind_rows(data_BiB, data_INMA, data_RHEA)



# drop fup_ column as it does not add any information
complete_data <- complete_data %>% select(-fup_)



# summary table of unmatched serology with anthropometry. Number of NAs in each column 


main_summary <- complete_data %>%
  filter(!is.na(age_days)) %>%
  summarise(
    across(c("cbmi", "cabdo", "csubscap", "ctriceps"), ~sum(is.na(.)), .names = "na_in_{.col}")
  ) %>%
  pivot_longer(cols = everything(), names_to = "column", values_to = "n_NA")

total_sample_row <- tibble(
  column = "total_sample",
  n_NA = nrow(complete_data)
)

na_sum_table <- bind_rows(main_summary, total_sample_row)



# compute serology - anthropometry age differences by timepoint

age_diff_summary <- complete_data %>%
  filter(!is.na(age_days), !is.na(agecd_cgrowth)) %>%
  mutate(age_diff = age_days - agecd_cgrowth) %>%
  group_by(timepoint) %>%
  summarise(
    mean_diff = mean(age_diff),
    sd_diff = sd(age_diff),
    min_diff = min(age_diff),
    max_diff = max(age_diff),
    n = n(),
    .groups = "drop"
  )

age_diff_summary




#### DESCRIPTIVE STATISTICS #### 



gt_tabble <- complete_data %>%
  mutate(cohort_label = case_when(
    coh == 1 ~ "BiB",
    coh == 2 ~ "INMA",
    coh == 3 ~ "RHEA",
    TRUE ~ as.character(coh)
  )) %>%
  group_by(cohort_label) %>%
  summarise(
    N = n(),
    Male = sum(sex == 1, na.rm = TRUE),
    Female = sum(sex == 2, na.rm = TRUE),
    BMI_Mean = mean(cbmi, na.rm = TRUE),
    BMI_SD = sd(cbmi, na.rm = TRUE),
    BMIz_Mean = mean(bmi_zscore, na.rm = TRUE),
    BMIz_SD = sd(bmi_zscore, na.rm = TRUE),
    Heightz_Mean = mean(height_zscore, na.rm = TRUE),
    Heightz_SD = sd(height_zscore, na.rm = TRUE),
    WHtR_Mean = mean(whtr, na.rm = TRUE),
    WHtR_SD = sd(whtr, na.rm = TRUE)
  ) %>%
  gt() %>%
  tab_header(
    title = "Descriptive Statistics by Cohort",
    subtitle = "Mean and SD of Anthropometric Measures"
  ) %>%
  tab_spanner(label = "BMI", columns = c(BMI_Mean, BMI_SD)) %>%
  tab_spanner(label = "BMI z-score", columns = c(BMIz_Mean, BMIz_SD)) %>%
  tab_spanner(label = "Height z-score", columns = c(Heightz_Mean, Heightz_SD)) %>%
  tab_spanner(label = "Waist-to-Height Ratio", columns = c(WHtR_Mean, WHtR_SD)) %>%
  fmt_number(
    columns = c(BMI_Mean, BMI_SD, BMIz_Mean, BMIz_SD, Heightz_Mean, Heightz_SD, WHtR_Mean, WHtR_SD),
    decimals = 2
  ) %>%
  cols_label(cohort_label = "Cohort") %>%
  tab_source_note(source_note = "Source: BiB, INMA, RHEA Cohorts")
