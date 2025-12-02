obesity_all_long$age <- obesity_all_long$agecd_cgrowth / 365

obesity_all_long <- obesity_all_long %>%
  filter(agecd_cgrowth != 0)


serostatus_all_long <- serostatus_all_long %>%
  filter(timepoint != 0)






propagate_sero_nn_multi <- function(anthro, sero, virus_vec) {
  
  anthro_out <- anthro  # start with the input dataset
  
  for (virus in virus_vec) {
    
    sero_col <- rep(NA_integer_, nrow(anthro))  # empty vector for this virus
    sero_split <- split(sero, sero$h_id)
    
    for(i in seq_len(nrow(anthro))) {
      
      h <- anthro$h_id[i]
      age_i <- anthro$age[i]
      
      if(!h %in% names(sero_split)) next
      s <- sero_split[[as.character(h)]]
      s_status <- s[[virus]]
      
      prev_idx <- which(s$age <= age_i)
      if(length(prev_idx) > 0) {
        prev_row <- prev_idx[which.max(s$age[prev_idx])]
        prev_age <- s$age[prev_row]
        prev_status <- s_status[prev_row]
      } else {
        prev_status <- NA
      }
      
      next_idx <- which(s$age >= age_i)
      if(length(next_idx) > 0) {
        next_row <- next_idx[which.min(s$age[next_idx])]
        next_age <- s$age[next_row]
        next_status <- s_status[next_row]
      } else {
        next_status <- NA
      }
      
      if(!is.na(prev_status) & !is.na(next_status)) {
        if(prev_status == next_status) {
          sero_col[i] <- prev_status
        } else if(prev_status == 0 & next_status == 1) {
          midpoint <- (prev_age + next_age)/2
          sero_col[i] <- ifelse(age_i < midpoint, 0L, 1L)
        }
      } else if(!is.na(prev_status) & is.na(next_status)) {
        sero_col[i] <- ifelse(prev_status == 1, 1L, NA_integer_)
      } else if(is.na(prev_status) & !is.na(next_status)) {
        sero_col[i] <- ifelse(next_status == 0, 0L, NA_integer_)
      }
    }
    
    # Attach the new column to the output dataset
    anthro_out[[paste0(virus, "_sero")]] <- sero_col
  }
  
  return(anthro_out)  # just one dataset with 9 new columns
}




virus_vec <- c("CMV_class", "cut_VZV", "cut_Avd36", "EBV_class",
               "cut_BK", "cut_JC", "cut_KI", "cut_WU", "cut_MCV")




# turn imputed into list of completed datasets
imputed_list <- lapply(1:imputed$m, function(i) complete(imputed, i))


# propagate for each imputed dataset
processed_list <- lapply(imputed_list, function(df) {
  propagate_sero_nn_multi(df, serostatus_all_long, virus_vec)
})











# 1. Identify log10 columns from the first dataframe
log_cols <- names(processed_list[[1]])[str_detect(names(processed_list[[1]]), "^log10_")]

# 2. Apply the rule to all 20 datasets
processed_list <- lapply(processed_list, function(df) {
  df %>%
    mutate(
      across(
        all_of(log_cols),
        ~ ifelse(is.na(sero_age), NA_real_, .)
      )
    )
})


#



closest_points <-propagate_sero_nn_multi(closest_points, serostatus_all_long, virus_vec)





# Example for Adv36; can be repeated for other viruses
bmi_sero_summary <- closest_points %>%
  select(h_id, age, timepoint, sex, coh, bmi_zscore, cbmi, bmi_category, cut_Avd36_sero) %>%
  filter(!is.na(cut_Avd36_sero)) %>%   # only children with serostatus
  group_by(timepoint, cut_Avd36_sero) %>%
  summarise(
    N = n(),
    BMI_Mean = mean(cbmi, na.rm = TRUE),
    BMI_SD = sd(cbmi, na.rm = TRUE),
    BMIz_Mean = mean(bmi_zscore, na.rm = TRUE),
    BMIz_SD = sd(bmi_zscore, na.rm = TRUE),
    Overweight_n = sum(bmi_category %in% c("Overweight", "Obese"), na.rm = TRUE),
    Overweight_pct = mean(bmi_category %in% c("Overweight", "Obese"), na.rm = TRUE)*100,
    .groups = "drop"
  ) %>%
  mutate(
    BMI = sprintf("%.1f (%.1f)", BMI_Mean, BMI_SD),
    BMIz = sprintf("%.2f (%.2f)", BMIz_Mean, BMIz_SD),
    Overweight = sprintf("%d (%.1f%%)", Overweight_n, Overweight_pct)
  ) %>%
  arrange(timepoint, cut_Avd36_sero)






bmi_sero_summary_cohort <- closest_points %>%
  select(h_id, age, timepoint, sex, coh, bmi_zscore, cbmi, bmi_category, cut_Avd36_sero) %>%
  filter(!is.na(cut_Avd36_sero)) %>%
  group_by(timepoint, coh, cut_Avd36_sero) %>%
  summarise(
    N = n(),
    BMI_Mean = mean(cbmi, na.rm = TRUE),
    BMI_SD = sd(cbmi, na.rm = TRUE),
    BMIz_Mean = mean(bmi_zscore, na.rm = TRUE),
    BMIz_SD = sd(bmi_zscore, na.rm = TRUE),
    Overweight_n = sum(bmi_category %in% c("Overweight", "Obese"), na.rm = TRUE),
    Overweight_pct = mean(bmi_category %in% c("Overweight", "Obese"), na.rm = TRUE)*100,
    .groups = "drop"
  ) %>%
  mutate(
    BMI = sprintf("%.1f (%.1f)", BMI_Mean, BMI_SD),
    BMIz = sprintf("%.2f (%.2f)", BMIz_Mean, BMIz_SD),
    Overweight = sprintf("%d (%.1f%%)", Overweight_n, Overweight_pct),
    cohort_label = case_when(
      coh == 1 ~ "BiB",
      coh == 2 ~ "INMA",
      coh == 3 ~ "RHEA",
      TRUE ~ as.character(coh)
    )
  ) %>%
  arrange(timepoint, cohort_label, cut_Avd36_sero)




library(ggplot2)

# Prepare data for plotting (mean per timepoint)
plot_data <- closest_points %>%
  select(h_id, age, timepoint, bmi_zscore, cut_Avd36_sero, coh) %>%
  filter(!is.na(cut_Avd36_sero)) %>%
  group_by(timepoint, cut_Avd36_sero, coh) %>%
  summarise(
    BMIz_mean = mean(bmi_zscore, na.rm = TRUE),
    BMIz_se = sd(bmi_zscore, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    cohort_label = case_when(
      coh == 1 ~ "BiB",
      coh == 2 ~ "INMA",
      coh == 3 ~ "RHEA",
      TRUE ~ as.character(coh)
    )
  )

ggplot(plot_data, aes(x = timepoint, y = BMIz_mean, color = factor(cut_Avd36_sero), group = cut_Avd36_sero)) +
  geom_line(linewidth = 1) +      # use linewidth instead of size
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = BMIz_mean - 1.96*BMIz_se,
                  ymax = BMIz_mean + 1.96*BMIz_se,
                  fill = factor(cut_Avd36_sero)), alpha = 0.2, color = NA) +
  facet_wrap(~cohort_label) +
  scale_x_continuous(breaks = 1:4, labels = c("2y", "4–6y", "6–9y", "11y")) +
  labs(x = "Timepoint", y = "BMI z-score", color = "Seropositive", fill = "Seropositive") +
  theme_minimal(base_size = 14)









sero_vars <- c("CMV_class_sero", "cut_VZV_sero", "cut_Avd36_sero",
               "EBV_class_sero", "cut_BK_sero", "cut_JC_sero",
               "cut_KI_sero", "cut_WU_sero", "cut_MCV_sero")

# Create summary tables for each serovar
bmi_sero_summaries <- map(sero_vars, function(sero) {
  closest_points %>%
    select(h_id, age, timepoint, sex, coh, bmi_zscore, cbmi, bmi_category, all_of(sero)) %>%
    filter(!is.na(.data[[sero]])) %>%
    group_by(timepoint, coh, .data[[sero]]) %>%
    summarise(
      N = n(),
      BMI_Mean = mean(cbmi, na.rm = TRUE),
      BMI_SD = sd(cbmi, na.rm = TRUE),
      BMIz_Mean = mean(bmi_zscore, na.rm = TRUE),
      BMIz_SD = sd(bmi_zscore, na.rm = TRUE),
      Overweight_n = sum(bmi_category %in% c("Overweight", "Obese"), na.rm = TRUE),
      Overweight_pct = mean(bmi_category %in% c("Overweight", "Obese"), na.rm = TRUE)*100,
      .groups = "drop"
    ) %>%
    mutate(
      BMI = sprintf("%.1f (%.1f)", BMI_Mean, BMI_SD),
      BMIz = sprintf("%.2f (%.2f)", BMIz_Mean, BMIz_SD),
      Overweight = sprintf("%d (%.1f%%)", Overweight_n, Overweight_pct),
      cohort_label = case_when(
        coh == 1 ~ "BiB",
        coh == 2 ~ "INMA",
        coh == 3 ~ "RHEA",
        TRUE ~ as.character(coh)
      ),
      sero_var = sero
    ) %>%
    arrange(timepoint, cohort_label, .data[[sero]])
})

# Combine into one table if needed
bmi_sero_summary_all <- bind_rows(bmi_sero_summaries)


# Create plots for each serovar
bmi_sero_plots <- map(sero_vars, function(sero) {
  plot_data <- closest_points %>%
    select(h_id, age, timepoint, bmi_zscore, coh, all_of(sero)) %>%
    filter(!is.na(.data[[sero]])) %>%
    group_by(timepoint, .data[[sero]], coh) %>%
    summarise(
      BMIz_mean = mean(bmi_zscore, na.rm = TRUE),
      BMIz_se = sd(bmi_zscore, na.rm = TRUE)/sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      cohort_label = case_when(
        coh == 1 ~ "BiB",
        coh == 2 ~ "INMA",
        coh == 3 ~ "RHEA",
        TRUE ~ as.character(coh)
      )
    )
  
  ggplot(plot_data, aes(x = timepoint, y = BMIz_mean, color = factor(.data[[sero]]), group = .data[[sero]])) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = BMIz_mean - 1.96*BMIz_se,
                    ymax = BMIz_mean + 1.96*BMIz_se,
                    fill = factor(.data[[sero]])), alpha = 0.2, color = NA) +
    facet_wrap(~cohort_label) +
    scale_x_continuous(breaks = 1:4, labels = c("2y", "4–6y", "6–9y", "11y")) +
    labs(x = "Timepoint", y = "BMI z-score", color = "Seropositive", fill = "Seropositive",
         title = paste("BMI z-score by", sero)) +
    theme_minimal(base_size = 14)
})

# Access the first plot as an example:
bmi_sero_plots[[1]]

# Or save all plots with a loop
# plotssss<- walk2(bmi_sero_plots, sero_vars, ~ ggsave(paste0(.y, "_bmi_z_plot.png"), .x, width = 8, height = 5))



library(dplyr)
library(gt)
library(tidyr)
library(stringr)

bmi_sero_gt_compact <- bmi_sero_summary_all %>%
  mutate(
    cohort_label = case_when(
      coh == 1 ~ "BiB",
      coh == 2 ~ "INMA",
      coh == 3 ~ "RHEA",
      TRUE ~ as.character(coh)
    ),
    serostatus_label = ifelse(.data[[names(.)[3]]] == 1, "Seropositive", "Seronegative")
  ) %>%
  select(sero_var, cohort_label, timepoint, serostatus_label, BMI, BMIz, Overweight) %>%
  group_by(sero_var, cohort_label, timepoint, serostatus_label) %>%
  summarise(
    BMI = first(BMI),
    BMIz = first(BMIz),
    Overweight = first(Overweight),
    .groups = "drop"
  ) %>%
  filter(!is.na(serostatus_label)) %>%

  pivot_wider(
    names_from = serostatus_label,
    values_from = c(BMI, BMIz, Overweight),
    names_glue = "{.value} ({serostatus_label})"
  ) %>%
  arrange(sero_var, cohort_label, timepoint) %>%
  gt(groupname_col = "sero_var") %>%
  tab_header(
    title = md("**BMI and BMI z-score by Serostatus, Cohort, and Timepoint**"),
    subtitle = md("Mean (SD) values and overweight prevalence (%) for seropositive and seronegative children")
  ) %>%
  cols_label(
    cohort_label = md("**Cohort**"),
    timepoint = md("**Timepoint**"),
    `BMI (Seronegative)` = md("**BMI (Mean (SD)) – Seronegative**"),
    `BMI (Seropositive)` = md("**BMI (Mean (SD)) – Seropositive**"),
    `BMIz (Seronegative)` = md("**BMI z-score (Mean (SD)) – Seronegative**"),
    `BMIz (Seropositive)` = md("**BMI z-score (Mean (SD)) – Seropositive**"),
    `Overweight (Seronegative)` = md("**Overweight/Obese (%) – Seronegative**"),
    `Overweight (Seropositive)` = md("**Overweight/Obese (%) – Seropositive**")
  ) %>%
  tab_spanner(
    label = md("**BMI (kg/m²)**"),
    columns = starts_with("BMI (")
  ) %>%
  tab_spanner(
    label = md("**BMI z-score**"),
    columns = starts_with("BMIz (")
  ) %>%
  tab_spanner(
    label = md("**Overweight/Obese**"),
    columns = starts_with("Overweight (")
  ) %>%
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(13),
    data_row.padding = px(3),
    table.width = pct(100)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

bmi_sero_gt_compact













