# Install packages
install.packages(c("anthro", "zscorer","haven","anthroplus", "tidyr", "dplyr"))

# Library packages
library(anthro)
library(zscorer)
library(haven)
library(anthroplus)
library(dplyr)
library(tidyr)
library(stringr)

##################### TO DO #####################
#  tell Marianna about covariates mother smoking, bmi, socioeconomic status, activity level, 






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





# initialize z score values
obesity_BiB[, c("bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore")] <- NA
obesity_INMA[, c("bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore")] <- NA
obesity_RHEA[, c("bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore")] <- NA






# check for NAs in agecd_growth(days) if yes calculate it with month *  30.4375 for INMA and RHEA
obesity_INMA$agecd_cgrowth[is.na(obesity_INMA$agecd_cgrowth)] <- round(obesity_INMA$agecm_cgrowth[is.na(obesity_INMA$agecd_cgrowth)] * 30.4375, 0)


# bmi calculation for INMA cohort
obesity_INMA$cbmi <- obesity_INMA$cweight / ((obesity_INMA$cheight / 100) ^ 2)

# calculate generate z-scores for weight, height, BMI, triceps skinfold and subscapular skinfold according to WHO 0-5 years old


assign_anthro_scores <- function(df) {
  rows <- !is.na(df$agecd_cgrowth) & df$agecd_cgrowth <= 1826
  
  # Run anthroplus_zscores() once on subset rows
  z_scores <- with(
    df[rows, ],
    anthro_zscores(
      sex = sex,
      age = agecd_cgrowth,
      weight = cweight,
      lenhei = cheight,
      triskin = ctriceps,
      subskin = csubscap
    )
  )
  
  target_cols <- c("bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore")
  source_cols <- c("zbmi", "zwei", "zlen", "zts", "zss")
  
  for (i in seq_along(target_cols)) {
    df[[target_cols[i]]][rows] <- z_scores[[source_cols[i]]]
  }
  
  return(df)
}


obesity_BiB <- assign_anthro_scores(obesity_BiB)
obesity_INMA <- assign_anthro_scores(obesity_INMA)
obesity_RHEA <- assign_anthro_scores(obesity_RHEA)







# calculate generate z-scores for weight, height, BMI according to WHO 5-19 years old


assign_anthroplus_scores <- function(df) {
  rows <- !is.na(df$agecd_cgrowth) & df$agecd_cgrowth > 1826
  
  # Run anthroplus_zscores() once on subset rows
  z_scores <- with(
    df[rows, ],
    anthroplus_zscores(
      sex = sex,
      age_in_months = agecm_cgrowth,
      weight_in_kg =  cweight,
      height_in_cm =  cheight,
      
    )
  )
  
  target_cols <- c("bmi_zscore", "weight_zscore", "height_zscore" )
  source_cols <- c("zbfa", "zwfa", "zhfa" )
  
  for (i in seq_along(target_cols)) {
    df[[target_cols[i]]][rows] <- z_scores[[source_cols[i]]]
  }
  
  return(df)
}



obesity_BiB <- assign_anthroplus_scores(obesity_BiB)
obesity_INMA <- assign_anthroplus_scores(obesity_INMA)
obesity_RHEA <- assign_anthroplus_scores(obesity_RHEA)




# function that finds and removes implausible values (±5 SDs) removes the whole row if all z scores are implausible values. 
# Creates a exclusion log for each dataset with h_id and values excluded.


clean_zscores <- function(df) {
  # Capture original object name as string
  df_name <- deparse(substitute(df))
  
  # Define z-score columns
  z_cols_under5 <- c("bmi_zscore", "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore")
  z_cols_5to19  <- c("bmi_zscore", "weight_zscore", "height_zscore")
  
  log_list <- list()
  remove_rows <- logical(nrow(df))  # track rows to remove
  
  for (i in seq_len(nrow(df))) {
    age <- df$agecd_cgrowth[i]
    id  <- df$h_id[i]
    
    if (is.na(age)) next
    
    z_cols <- if (age <= 1826) z_cols_under5 else z_cols_5to19
    z_values <- unlist(df[i, z_cols])
    
    outlier_flags <- abs(z_values) >= 5 & !is.na(z_values)
    
    if (all(outlier_flags)) {
      remove_rows[i] <- TRUE
      log_list[[length(log_list) + 1]] <- data.frame(
        h_id = id,
        agecd_cgrowth = age,
        removed = paste(z_cols, collapse = ","),
        values = paste(round(z_values, 2), collapse = ","),
        action = "row_removed"
      )
    } else if (any(outlier_flags)) {
      outlier_cols <- z_cols[outlier_flags]
      outlier_vals <- z_values[outlier_flags]
      df[i, outlier_cols] <- NA
      log_list[[length(log_list) + 1]] <- data.frame(
        h_id = id,
        agecd_cgrowth = age,
        removed = paste(outlier_cols, collapse = ","),
        values = paste(round(outlier_vals, 2), collapse = ","),
        action = "values_set_to_NA"
      )
    }
  }
  
  # Remove flagged rows
  df <- df[!remove_rows, ]
  
  # Create exclusion log
  exclusion_log <- if (length(log_list) > 0) do.call(rbind, log_list) else data.frame()
  
  # Overwrite the original dataset with cleaned version
  assign(df_name, df, envir = .GlobalEnv)
  
  # Save exclusion log with dynamic name
  assign(paste0("exclusionlog_", df_name), exclusion_log, envir = .GlobalEnv)
}



clean_zscores(obesity_BiB)
clean_zscores(obesity_INMA)
clean_zscores(obesity_RHEA)



# function to classify children as Normal, Overweight, Obese based on z score , depending on age < 5 and age > 5



classify_bmi_category <- function(df) {
  df <- df %>%
    mutate(
      bmi_category = case_when(
        is.na(agecd_cgrowth) | is.na(bmi_zscore) ~ NA_character_,
        
        # Under 5 years
        agecd_cgrowth <= 1826 & bmi_zscore < 2 ~ "Normal",
        agecd_cgrowth <= 1826 & bmi_zscore >= 2 & bmi_zscore < 3 ~ "Overweight",
        agecd_cgrowth <= 1826 & bmi_zscore >= 3 ~ "Obese",
        
        # 5 years and older
        agecd_cgrowth > 1826 & bmi_zscore < -3 ~ "Severely thin",
        agecd_cgrowth > 1826 & bmi_zscore >= -3 & bmi_zscore < -2 ~ "Thin",
        agecd_cgrowth > 1826 & bmi_zscore >= -2 & bmi_zscore <= 1 ~ "Normal",
        agecd_cgrowth > 1826 & bmi_zscore > 1 & bmi_zscore <= 2 ~ "Overweight",
        agecd_cgrowth > 1826 & bmi_zscore > 2 ~ "Obese",
        
        TRUE ~ NA_character_
      )
    )
  
  return(df)
}


obesity_BiB <- classify_bmi_category(obesity_BiB)
obesity_INMA <- classify_bmi_category(obesity_INMA)
obesity_RHEA <- classify_bmi_category(obesity_RHEA)





# indetify central obesity for children over 4 years of age.  

classify_central_obesity <- function(df) {
  df <- df %>%
    mutate(
      whtr = cabdo / cheight,  # assumes waist and height are in cm
      central_obesity = case_when(
        is.na(agecd_cgrowth) | is.na(cabdo) | is.na(cheight) ~ NA,
        agecd_cgrowth < 1460 ~ NA,  # under 4 years, not applicable
        whtr >= 0.5 ~ TRUE,
        whtr < 0.5 ~ FALSE,
        TRUE ~ NA
      )
    )
  
  return(df)
}


obesity_BiB <- classify_central_obesity(obesity_BiB)
obesity_INMA <- classify_central_obesity(obesity_INMA)
obesity_RHEA <- classify_central_obesity(obesity_RHEA)




# function to classify stunted using height for age z-score

classify_stunting <- function(df) {
  df <- df %>%
    mutate(
      stunted = case_when(
        is.na(height_zscore) ~ NA,
        height_zscore < -2 ~ TRUE,
        height_zscore >= -2 ~ FALSE
      )
    )
  
  return(df)
}

obesity_BiB <- classify_stunting(obesity_BiB)
obesity_INMA <- classify_stunting(obesity_INMA)
obesity_RHEA <- classify_stunting(obesity_RHEA)






