library(anthro)
library(zscorer)
library(haven)
library(anthroplus)
library(dplyr)
library(tidyr)
library(stringr)
library(gt)




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



# function to make serology to long format

longify_serology_0_4 <- function(df) {
  # match variables ending in a digit (0–4)
  measurement_cols <- names(df) %>% str_subset("\\d$")
  
  # extract roots (everything except last digit)
  variable_roots <- measurement_cols %>%
    str_remove("\\d$") %>%
    unique()
  
  # build regex pattern to split root and timepoint
  pattern_regex <- paste0("^(", paste(variable_roots, collapse = "|"), ")([0-4])$")
  
  df_long <- df %>%
    pivot_longer(
      cols = all_of(measurement_cols),
      names_to = c(".value", "timepoint"),
      names_pattern = pattern_regex
    ) %>%
    {
      existing_vars <- intersect(variable_roots, names(.))
      filter(., !if_all(all_of(existing_vars), is.na))
    } %>%
    filter(fup_ == 1)
  
  return(df_long)
}







# function to match obesity data to serology data
# for each serology timepoint measurment the functions finds the best match
# based on BMI existence, then other antrho variables count("cabdo", "csubscap", "ctriceps") 
# and finally by measuremnt age difference with a cuttoff of 60 days.
match_growth_to_serology <- function(
    serology,
    obesity,
    anthropo_vars = c("cbmi", "cabdo", "csubscap", "ctriceps"),
    extra_vars = c(
      "sex",
      "agecd_cgrowth",
      "cheight",
      "agecm_cgrowth",
      "agecy_cgrowth",
      "cweight",
      "cabdo",
      "csubscap",
      "ctriceps",
      "bmi_zscore",
      "weight_zscore",
      "height_zscore",
      "tris_zscore",
      "subscap_zscore",
      "bmi_category",
      "whtr",
      "central_obesity",
      "stunted"
    ),
    cutoff_days = 60
) {
  
  obesity <- obesity %>%
    mutate(orig_row = row_number())
  
  count_anthro <- function(df) {
    rowSums(!is.na(df[, anthropo_vars, drop=FALSE]))
  }
  
  obesity <- obesity %>%
    mutate(
      anthro_count = count_anthro(.)
    ) %>%
    filter(anthro_count > 0)
  
  matched_indices <- integer(nrow(serology))
  matched_indices[] <- NA_integer_
  
  debug_matches <- list()
  
  for (i in seq_len(nrow(serology))) {
    ser_row <- serology[i, ]
    ser_age <- ser_row$age_days
    ser_cid <- ser_row$h_id
    
    if (is.na(ser_age)) next
    
    candidates <- obesity %>%
      filter(h_id == ser_cid) %>%
      mutate(age_diff = ceiling(abs(agecd_cgrowth - ser_age))) %>%
      filter(age_diff <= cutoff_days)
    
    if (nrow(candidates) == 0) next
    
    candidates <- candidates %>%
      mutate(
        cbmi_present = !is.na(cbmi),
        other_anthro_count = rowSums(!is.na(select(., setdiff(anthropo_vars, "cbmi"))))
      ) %>%
      arrange(
        desc(cbmi_present),
        desc(other_anthro_count),
        age_diff
      )
    
    best_match <- candidates[1, ]
    matched_indices[i] <- best_match$orig_row
    
    debug_matches[[length(debug_matches) + 1]] <- candidates %>%
      mutate(serology_row = i, matched = (orig_row == best_match$orig_row))
  }
  
  serology$matched_growth_row <- matched_indices
  
  # Select all required variables for merging
  vars_to_select <- unique(c("orig_row", extra_vars, anthropo_vars))
  
  matched_growth_data <- obesity %>%
    select(all_of(vars_to_select))
  
  # Merge with serology
  serology_merged <- serology %>%
    left_join(matched_growth_data, by = c("matched_growth_row" = "orig_row"))
  
  # Insert extra_vars then anthropo_vars after 12th column
  insert_after <- 12
  original_cols <- names(serology)
  new_extra <- setdiff(extra_vars, original_cols)
  new_anthro <- setdiff(anthropo_vars, original_cols)
  
  insert_order <- append(
    names(serology_merged)[1:insert_after],
    c(new_extra, new_anthro),
    after = insert_after
  )
  
  final_cols <- unique(c(insert_order, names(serology_merged)))  # ensure all columns are included
  
  serology_merged <- serology_merged[, final_cols]
  
  return(serology_merged)
}
