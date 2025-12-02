
install.packages(c("anthro", "zscorer","haven","anthroplus", "tidyr", "dplyr", "gt", "rstatix", "janitor", "purrr", "ggplot2", 
                   "gtsummary","binom", "lme4", "broom.mixed", "lmerTest" , "mice", "randomForest"))

# Library packages
library(anthro)
library(zscorer)
library(haven)
library(anthroplus)
library(dplyr)
library(tidyr)
library(stringr)
library(gt)
library(rstatix)
library(janitor)
library(purrr)
library(ggplot2)
library(gtsummary)
library(binom)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(mice)
library(lmerTest)
library(randomForest)

source("functions.R")
# call functions script

pdf("everything.pdf", width = 8, height = 6)


# commands to push rscript changes to git hub
# git init
# git add MScTH_rscript.R
# git push origin main





# TO DO

# 6 datasets: 
# 3 wide split by coh with everything(maybe 1 pooled wide)
# 1 long anthropometry
# 1 long serostatus


# plot trajectories
# .pdf to export all plots







# dag, missing, gt tables (add categorical data) timepoint -> cohort , repeated measueres 
# https://www.dagitty.net/dags.html
# https://academic.oup.com/ije/article/50/2/620/6012812



# write.csv(complete_data, "C:/Users/grigoris.kalampoukas/MScThesis/complete_data.csv", row.names = FALSE)




# data loading

# anthropometric data

obesity_BiB <- read_dta("G:/Msc Thesis/Grigoris/marato_obesity_BiB.dta")
obesity_RHEA <-  read_dta("G:/Msc Thesis/Grigoris/marato_obesity_RHEA.dta")
obesity_INMA <-  read_dta("G:/Msc Thesis/Grigoris/marato_obesity_INMA.dta")

# serological data
serology_BiB <- read_dta("G:/Msc Thesis/Grigoris/marato_serology_Bib.dta")
serology_RHEA <- read_dta("G:/Msc Thesis/Grigoris/marato_serology_Rhea.dta")
serology_INMA <-  read_dta("G:/Msc Thesis/Grigoris/marato_serology_INMA.dta")

#cofounders

confounders <- read_dta("D:/Downloads Real/variables_grigoris.dta")

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

complete_data_bckp <- complete_data

complete_data <- complete_data[complete_data$timepoint != 0, ]


# drop fup_ column as it does not add any information
complete_data <- complete_data %>% select(-fup_)

# make sex as factor
complete_data$sex <- factor(
  complete_data$sex,
  levels = c(1, 2),
  labels = c("Male", "Female")
)

# Convert coh to factor with labeled levels
complete_data$coh <- factor(complete_data$coh,
                               levels = c(1, 2, 3),
                               labels = c("BIB", "INMA", "RHE"))



# summary table of unmatched serology with anthropometry. Number of NAs in each column 

total_sample_row <- tibble(
  column = "total_sample",
  n_NA = nrow(data_INMA %>% filter(timepoint == 4))
)




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




complete_data$timepoint <- as.factor(complete_data$timepoint)




### ERROR? ###

main_summary_INMA <- data_INMA %>%
  filter(timepoint == 4, !is.na(age_days)) %>%
  summarise(
    across(c("cbmi", "cabdo", "csubscap", "ctriceps"), ~sum(is.na(.)), .names = "na_in_{.col}")
  ) %>%
  pivot_longer(cols = everything(), names_to = "column", values_to = "n_NA")

total_sample_row_INMA <- tibble(
  column = "total_sample",
  n_NA = nrow(data_INMA %>% filter(timepoint == 4))
)

na_sum_table_INMA <- bind_rows(main_summary, total_sample_row)






#### DESCRIPTIVE STATISTICS #### 

anthro_long_sero_wide <- anthro_long_sero_wide %>%
  rename(sex = sex.y) %>%
  mutate(
    sex = factor(sex, levels = c(1,2), labels = c("Male", "Female")),
    # create timepoint based on age
    timepoint = case_when(
      age >= 1.5 & age < 3.5  ~ 1,
      age >= 4   & age < 6  ~ 2,
      age >= 6   & age < 10 ~ 3,
      age >= 10  & age <= 12 ~ 4,
      TRUE ~ NA_real_
    ),
    timepoint_label = case_when(
      timepoint == 1 ~ "2 Years",
      timepoint == 2 ~ "4–6 Years",
      timepoint == 3 ~ "6–9 Years",
      timepoint == 4 ~ "11 Years",
      TRUE ~ as.character(timepoint)
    )
  )


summary_stats <- serostatus_all_long %>%
  group_by(timepoint) %>%
  summarise(
    median_age = median(age, na.rm = TRUE),
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE)
  )

# extract unique h_ids
unique_ids <- unique(serostatus_all_long$h_id)


closest_points <- anthro_long_sero_wide %>%
  group_by(h_id) %>%
  group_modify(~ {
    df_sub <- .
    assigned <- lapply(1:nrow(summary_stats), function(i) {
      # Get timepoint info
      median_val <- summary_stats$median_age[i]
      min_val <- summary_stats$min_age[i]
      max_val <- summary_stats$max_age[i]
      
      # Filter rows within the allowed interval
      df_valid <- df_sub %>% filter(age >= min_val, age <= max_val)
      
      # If no valid measurement, skip
      if(nrow(df_valid) == 0) return(NULL)
      
      # Pick closest to median among valid rows
      idx <- which.min(abs(df_valid$age - median_val))
      df_row <- df_valid[idx, ]
      df_row$timepoint <- summary_stats$timepoint[i]
      df_row
    })
    bind_rows(assigned)
  }) %>%
  ungroup()




# anthropometry table
library(gt)
library(dplyr)

# anthropometry table using closest_points
anthro_table <- closest_points %>%
  filter(!is.na(timepoint)) %>%  # ensure we only include assigned timepoints
  mutate(
    cohort_label = case_when(
      coh == 1 ~ "BiB",
      coh == 2 ~ "INMA",
      coh == 3 ~ "RHEA",
      TRUE ~ as.character(coh)
    ),
    timepoint_label = case_when(
      timepoint == 1 ~ "2 Years",
      timepoint == 2 ~ "4–6 Years",
      timepoint == 3 ~ "6–9 Years",
      timepoint == 4 ~ "11 Years",
      TRUE ~ as.character(timepoint)
    )
  ) %>%
  group_by(timepoint, timepoint_label, cohort_label) %>%
  summarise(
    N = n(),
    Male = sum(sex == "Male", na.rm = TRUE),
    Female = sum(sex == "Female", na.rm = TRUE),
    
    # BMI
    BMI_Mean = mean(cbmi, na.rm = TRUE),
    BMI_SD = sd(cbmi, na.rm = TRUE),
    BMI_n = sum(!is.na(cbmi)),
    
    # BMI z-score
    BMIz_Mean = mean(bmi_zscore, na.rm = TRUE),
    BMIz_SD = sd(bmi_zscore, na.rm = TRUE),
    BMIz_n = sum(!is.na(bmi_zscore)),
    
    # Height z-score
    Heightz_Mean = mean(height_zscore, na.rm = TRUE),
    Heightz_SD = sd(height_zscore, na.rm = TRUE),
    Heightz_n = sum(!is.na(height_zscore)),
    
    # WHtR
    WHtR_Mean = mean(whtr, na.rm = TRUE),
    WHtR_SD = sd(whtr, na.rm = TRUE),
    WHtR_n = sum(!is.na(whtr)),
    
    # Categorical
    Stunted_n = sum(stunted == 1, na.rm = TRUE),
    Stunted_pct = mean(stunted == 1, na.rm = TRUE) * 100,
    
    Thin_n = sum(bmi_category == "Thin", na.rm = TRUE),
    Thin_pct = mean(bmi_category == "Thin", na.rm = TRUE) * 100,
    
    Normal_n = sum(bmi_category == "Normal", na.rm = TRUE),
    Normal_pct = mean(bmi_category == "Normal", na.rm = TRUE) * 100,
    
    Overweight_n = sum(bmi_category == "Overweight", na.rm = TRUE),
    Overweight_pct = mean(bmi_category == "Overweight", na.rm = TRUE) * 100,
    
    Obese_n = sum(bmi_category == "Obese", na.rm = TRUE),
    Obese_pct = mean(bmi_category == "Obese", na.rm = TRUE) * 100,
    
    CentralOb_n = sum(central_obesity == 1, na.rm = TRUE),
    CentralOb_pct = mean(central_obesity == 1, na.rm = TRUE) * 100,
    
    .groups = "drop"
  ) %>%
  mutate(
    cohort_label = sprintf("%s (n = %d)", cohort_label, N),
    
    BMI = sprintf("%.1f (%.1f)\nn=%d", BMI_Mean, BMI_SD, BMI_n),
    BMIz = sprintf("%.2f (%.2f)\nn=%d", BMIz_Mean, BMIz_SD, BMIz_n),
    Heightz = sprintf("%.2f (%.2f)\nn=%d", Heightz_Mean, Heightz_SD, Heightz_n),
    WHtR = sprintf("%.2f (%.2f)\nn=%d", WHtR_Mean, WHtR_SD, WHtR_n),
    
    Stunted = sprintf("%d (%.1f%%)", Stunted_n, Stunted_pct),
    Thin = sprintf("%d (%.1f%%)", Thin_n, Thin_pct),
    Normal = sprintf("%d (%.1f%%)", Normal_n, Normal_pct),
    Overweight = sprintf("%d (%.1f%%)", Overweight_n, Overweight_pct),
    Obese = sprintf("%d (%.1f%%)", Obese_n, Obese_pct),
    CentralOb = sprintf("%d (%.1f%%)", CentralOb_n, CentralOb_pct)
  ) %>%
  arrange(timepoint, cohort_label) %>%
  select(cohort_label, timepoint_label, Male, Female,
         BMI, BMIz, Heightz, WHtR,
         Stunted, Thin, Normal, Overweight, Obese, CentralOb) %>%
  gt(groupname_col = "timepoint_label") %>%
  tab_header(
    title = "Anthropometric Measures by Cohort and Timepoint",
    subtitle = "Stratified Summary per Follow-Up Period"
  ) %>%
  tab_spanner(label = "BMI", columns = BMI, id = "spanner_bmi") %>%
  tab_spanner(label = "BMI z-score", columns = BMIz, id = "spanner_bmiz") %>%
  tab_spanner(label = "Height z-score", columns = Heightz, id = "spanner_heightz") %>%
  tab_spanner(label = "Waist-to-Height Ratio", columns = WHtR, id = "spanner_whtr") %>%
  cols_label(
    timepoint_label = "Timepoint",
    Male = "Male",
    Female = "Female",
    BMI = "Mean (SD)",
    BMIz = "Mean (SD)",
    Heightz = "Mean (SD)",
    WHtR = "Mean (SD)",
    Stunted = "Stunted",
    Thin = "Thin",
    Normal = "Normal",
    Overweight = "Overweight",
    Obese = "Obese",
    CentralOb = "Central Obesity"
  ) %>%
  tab_source_note(source_note = "Stratified by cohort and timepoint (BiB, INMA, RHEA)")



anthro_table







# Serostatus table
sero_table <- serostatus_all_long%>%
  filter(timepoint != 0) %>%
  mutate(
    cohort_label = case_when(
      coh == 1 ~ "BiB",
      coh == 2 ~ "INMA",
      coh == 3 ~ "RHEA",
      TRUE ~ as.character(coh)
    ),
    timepoint_label = case_when(
      timepoint == 1 ~ "2 Years",
      timepoint == 2 ~ "4–6 Years",
      timepoint == 3 ~ "6–9 Years",
      timepoint == 4 ~ "11 Years",
      TRUE ~ as.character(timepoint)
    )
  ) %>%
  group_by(timepoint, cohort_label, timepoint_label) %>%
  summarise(
    N = n(),
    CMV_Positive = mean(CMV_class == 1, na.rm = TRUE) * 100,
    VZV_Positive = mean(cut_VZV == 1, na.rm = TRUE) * 100,
    ADV36_Positive = mean(cut_Avd36 == 1, na.rm = TRUE) * 100,
    EBV_Positive = mean(EBV_class == 1, na.rm = TRUE) * 100,
    BK_Positive = mean(cut_BK == 1, na.rm = TRUE) * 100,
    JC_Positive = mean(cut_JC == 1, na.rm = TRUE) * 100,
    KI_Positive = mean(cut_KI == 1, na.rm = TRUE) * 100,
    WU_Positive = mean(cut_WU == 1, na.rm = TRUE) * 100,
    MCV_Positive = mean(cut_MCV == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(timepoint, cohort_label) %>%
  select(cohort_label, timepoint_label, N, everything(), -timepoint) %>%
  gt(groupname_col = "timepoint_label") %>%
  tab_header(
    title = "Seropositivity by Cohort and Timepoint",
    subtitle = "Stratified Summary per Follow-Up Period"
  ) %>%
  tab_spanner(label = "Seropositivity (%)",
              columns = c(CMV_Positive, VZV_Positive, ADV36_Positive, EBV_Positive,
                          BK_Positive, JC_Positive, KI_Positive, WU_Positive, MCV_Positive)) %>%
  fmt_number(columns = N, decimals = 0) %>%
  fmt_number(columns = c(CMV_Positive, VZV_Positive, ADV36_Positive, EBV_Positive,
                         BK_Positive, JC_Positive, KI_Positive, WU_Positive, MCV_Positive),
             decimals = 1, suffixing = FALSE) %>%
  fmt_percent(columns = c(CMV_Positive, VZV_Positive, ADV36_Positive, EBV_Positive,
                          BK_Positive, JC_Positive, KI_Positive, WU_Positive, MCV_Positive),
              decimals = 1, scale_values = FALSE) %>%
  cols_label(
    timepoint_label = "Timepoint",
    N = "N",
    CMV_Positive = "CMV",
    VZV_Positive = "VZV",
    ADV36_Positive = "ADV-36",
    EBV_Positive = "EBV",
    BK_Positive = "BK",
    JC_Positive = "JC",
    KI_Positive = "KI",
    WU_Positive = "WU",
    MCV_Positive = "MCV"
  ) %>%
  tab_source_note(source_note = "Stratified by cohort and timepoint (BiB, INMA, RHEA)")


sero_table





# histograms for normality check


norm_histograms(anthro_long_sero_wide, "cbmi")
norm_histograms(anthro_long_sero_wide, "bmi_zscore")
norm_histograms(anthro_long_sero_wide, "height_zscore")
norm_histograms(anthro_long_sero_wide, "weight_zscore")


complete_data %>%
  filter(!is.na(whtr)) %>%
  mutate(
    cohort_label = case_when(
      coh == 1 ~ "BiB",
      coh == 2 ~ "INMA",
      coh == 3 ~ "RHEA",
      TRUE ~ as.character(coh)
    )
  ) %>%
  ggplot(aes(x = whtr)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "skyblue") +
  facet_wrap(~ cohort_label, scales = "free_y") +
  labs(
    title = "Distribution of whtr by Cohort",
    x = "whtr",
    y = "Count"
  ) +
  theme_minimal()


# scatter plots. distribution of variable by cohort

scatter_by_cohort(anthro_long_sero_wide, "age", "cbmi")
scatter_by_cohort(anthro_long_sero_wide, "age", "bmi_zscore")
scatter_by_cohort(anthro_long_sero_wide, "age", "height_zscore")
scatter_by_cohort(anthro_long_sero_wide,"age", "weight_zscore")
scatter_by_cohort(anthro_long_sero_wide, "age", "whtr" )




# normality tests
norm_cbmi <- shapiro_by_group(anthro_long_sero_wide, "cbmi")
norm_bmiz <- shapiro_by_group(complete_data, "bmi_zscore")
norm_heightz <- shapiro_by_group(complete_data, "height_zscore")
norm_weight_z <- shapiro_by_group(complete_data, "weight_zscore")
norm_whtr <- shapiro_by_group(complete_data, "whtr")












# repeated measures vector

repeated_participants <- complete_data %>%
  group_by(h_id) %>%
  summarise(timepoint_count = n_distinct(timepoint)) %>%
  filter(timepoint_count > 1) %>%
  pull(h_id)



# plot participants repeated counts

timepoint_counts <- complete_data %>%
  group_by(h_id) %>%
  summarise(timepoint_count = n_distinct(timepoint))

ggplot(timepoint_counts, aes(x = factor(timepoint_count))) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Number of Participants by Count of Measured Timepoints",
    x = "Number of Timepoints Measured",
    y = "Count of Participants"
  ) +
  theme_minimal()







########## DATASETS TRANSFORMATION ########## 


# widen anthropometric datasets



# 3 wide datasets split by coh with everything(maybe 1 pooled wide)

data_BiB_wide <- widen2(data_BiB)
data_INMA_wide <- widen2(data_INMA)
data_RHEA_wide <- widen2(data_RHEA)


complete_data_wide <- bind_rows(data_BiB_wide, data_INMA_wide, data_RHEA_wide)



# all anthropomety in long format

cols_to_keep <- c(
  "h_id",  "sex",  "coh",
  "agecd_cgrowth", "cheight", "cweight", "cbmi", "cabdo", "csubscap", "ctriceps", 
  "bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore",
  "bmi_category", "whtr", "central_obesity", "stunted"
)
obesity_BiB$coh <- 1
obesity_INMA$coh <- 2
obesity_RHEA$coh <- 3


obesity_all_long <- bind_rows(
  obesity_BiB %>% select(all_of(cols_to_keep)),
  obesity_INMA %>% select(all_of(cols_to_keep)),
  obesity_RHEA %>% select(all_of(cols_to_keep))
)





serostatus_all_long <- bind_rows(serology_BiB, serology_INMA, serology_RHEA)






# 6 datasets: 
# 3 wide split by coh with everything(maybe 1 pooled wide) DONE!
# 1 long anthropometry DONE!
# 1 long serostatus done!


# Summary:
# complete_data = long format, all cohorts, serostatus with matched anthropometry
# data_$coh$_wide = wide format, by cohort, serostatus with matched anthropometry
# complete_data_wide = wide format, all cohorts, serostatus with matched anthropometry

# obesity_all = long format, all cohorts, all anthropometry
# serostatus_all_long = long format, all cohorts, all serology


# make sure timepoint is numeric or factor with defined levels
serostatus_all_long <- serostatus_all_long %>%
  mutate(timepoint = as.numeric(timepoint))




# define which columns to widen (everything except IDs & grouping vars)
value_vars <- setdiff(
  colnames(serostatus_all_long),
  c("h_id", "m_id", "c_id", "helixid", "cohort", "coh",
    "n_samples", "fup_", "date")   
)

serostatus_all_wide <- serostatus_all_long %>%
  pivot_wider(
    id_cols = h_id,
    names_from = timepoint,
    values_from = all_of(value_vars),
    names_glue = "{.value}_{timepoint}",
    names_sort = TRUE,                  # ensure 0,1,2 order
    names_vary = "fastest"                 # consistent column order
  )


anthro_long_sero_wide <- obesity_all_long %>%
  left_join(serostatus_all_wide, by = "h_id")


anthro_long_sero_wide <- anthro_long_sero_wide %>%
  left_join(confounders, by = "h_id") %>%
  relocate(all_of(c("pre_bmi_c", "urb_area_id",
                    "smk_p", "parity_m", "ethn3_m", "edu_m_0", "birth_weight","birth_length", "birth_head_circum",
                    "delivery_type", "m_age", "sex.y", "breastfed_ever", "preg_smk", "nursery_upto2years" )), .after = 3)



## to do:

# add obesity all add serostatus in wide format on each id of participant DONE!

# methods, serology , obesity

# read bibliography 


# quick fixes

summary(complete_data$age)

complete_data <- complete_data %>%
  mutate(
    agecd_cgrowth = if_else(is.na(agecd_cgrowth), round(age * 365.25), agecd_cgrowth),
    agecm_cgrowth = if_else(is.na(agecm_cgrowth), round(age * 12), agecm_cgrowth),
    agecy_cgrowth = if_else(is.na(agecy_cgrowth), round(age), agecy_cgrowth)
  )






# Model 1: BMI predicting viral infection(ADV36)
model1 <- glmer(
  cut_Avd36 ~ scale(bmi_zscore) * age + scale(cbmi) * cohort + (1 | h_id),
  data = complete_data,
  family = binomial(link = "logit")
)


summary(model1)

# Model 2: ADV36 predicting BMI
model2 <- lmer(
  bmi_zscore ~ cut_Avd36 * age + cut_Avd36 * cohort + (1 | h_id),
  data = complete_data
)

summary(model2)


# add interaction and cohort(drop age + sex)





# --- Virus variable names ---
virus_vars <- c("CMV_class", "cut_VZV", "cut_Avd36", "EBV_class",
                "cut_BK", "cut_JC", "cut_KI", "cut_WU", "cut_MCV")


# --- Function to run models and format results ---

# Run models
all_results <- run_bidirectional_models_table(complete_data, virus_vars)

# Format for display (drop raw numeric columns afterwards)
all_results_fmt <- all_results %>%
  mutate(
    p_value_fmt = case_when(
      p_value == 0 ~ "0",
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    # add stars for significance
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    p_value_fmt = paste0(p_value_fmt, sig),
    
    OR_or_Beta_fmt = sprintf("%.2f", OR_or_Beta),
    Estimate_fmt = sprintf("%.2f", Estimate),
    Std_Error_fmt = sprintf("%.2f", Std_Error),
    Random_SD_fmt = sprintf("%.2f", Random_SD)
  ) %>%
  select(Virus, Model, Predictor,
         Estimate_fmt, Std_Error_fmt, OR_or_Beta_fmt,
         p_value_fmt, Random_SD_fmt)

# --- GT table ---
results_gt <- all_results_fmt %>%
  gt(groupname_col = "Virus") %>%
  tab_header(
    title = "Bidirectional Mixed Models: Viral Infections and BMI (with Age Interaction)",
    subtitle = "Estimates, Odds Ratios/Beta, Standard Errors, p-values, Random Effects SD"
  ) %>%
  tab_spanner(
    label = "Fixed Effects",
    columns = c("Predictor", "Estimate_fmt", "Std_Error_fmt", "OR_or_Beta_fmt", "p_value_fmt")
  ) %>%
  cols_label(
    Predictor = "Predictor",
    Estimate_fmt = "Estimate",
    Std_Error_fmt = "Std. Error",
    OR_or_Beta_fmt = "OR/Beta",
    p_value_fmt = "p-value",
    Random_SD_fmt = "Random SD",
    Model = "Model Direction"
  ) %>%
  cols_merge(
    columns = c("Estimate_fmt", "Std_Error_fmt"),
    pattern = "{1} ({2})"
  ) %>%
  cols_move(columns = "Random_SD_fmt", after = "p_value_fmt") %>%
  opt_table_font(font = list(google_font("Roboto")))

results_gt






# group by timepoint and adv36 positive/negative
# check normality and based on normality check mean difference

data_filtered <- complete_data %>%
  filter(timepoint != 0)

test_results <- data_filtered %>%
  group_by(timepoint) %>%
  group_modify(~ {
    df <- .
    cbmi_pos <- df$cbmi[df$cut_Avd36 == 1]
    cbmi_neg <- df$cbmi[df$cut_Avd36 == 0]
    
    n_pos <- length(cbmi_pos)
    n_neg <- length(cbmi_neg)
    
    shapiro_pos <- if(n_pos >= 3) shapiro.test(cbmi_pos)$p.value else NA_real_
    shapiro_neg <- if(n_neg >= 3) shapiro.test(cbmi_neg)$p.value else NA_real_
    
    # Decide test
    test_method <- case_when(
      n_pos < 2 | n_neg < 2 ~ NA_character_,
      !is.na(shapiro_pos) & !is.na(shapiro_neg) & shapiro_pos > 0.05 & shapiro_neg > 0.05 ~ "t.test",
      TRUE ~ "wilcox.test"
    )
    
    # Compute mean difference
    mean_diff <- if(n_pos >= 1 & n_neg >= 1) mean(cbmi_pos, na.rm = TRUE) - mean(cbmi_neg, na.rm = TRUE) else NA_real_
    
    
    # Compute p-value
    p_val <- case_when(
      test_method == "t.test" ~ t.test(cbmi_pos, cbmi_neg, var.equal = FALSE)$p.value,
      test_method == "wilcox.test" ~ wilcox.test(cbmi_pos, cbmi_neg)$p.value,
      TRUE ~ NA_real_
    )
    
    tibble(
      n_pos = n_pos,
      n_neg = n_neg,
      shapiro_pos = shapiro_pos,
      shapiro_neg = shapiro_neg,
      test = test_method,
      mean_diff = mean_diff,
      p_value = p_val
    )
  }) %>%
  ungroup()

test_results




dev.off()

