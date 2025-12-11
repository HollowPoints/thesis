anthro_long_sero_wide$age <- anthro_long_sero_wide$agecd_cgrowth / 365
densityplot(anthro_long_sero_wide,  ~ age)

anthro_long_sero_wide <- anthro_long_sero_wide %>%
  mutate(bmi_scaled = (cbmi - min(cbmi, na.rm = TRUE)) / (max(cbmi, na.rm = TRUE) - min(cbmi, na.rm = TRUE)))




database <- anthro_long_sero_wide





database  <- database  %>%
  mutate(
    sex.x = as_factor(sex.x),
    delivery_type = as_factor(delivery_type),
    across(where(is.labelled) & !c("sex.x", "delivery_type"), ~ as.numeric(.))
  )

database$delivery_type <- as.numeric(database$delivery_type)

database <- classify_bmi_category(database)





# Parameters
viruses <- c("CMV", "VZV", "Avd36", "EBV", "BKPyV", "JCPyV", "KIPyV", "WUPyV", "MCPyV", "HSV" ,"EVB")
virus_pattern <- paste0("^log10_(", paste(viruses, collapse = "|"), ").*norm_\\d+$")
day_range <- 60 / 365.25  # 60 days in years

# 1) Add row id
database2 <- database %>% mutate(.rowid = row_number())

# 2) Pivot long all virus columns
long_df <- database2 %>%
  select(.rowid, age, starts_with("age_"), matches(virus_pattern)) %>%
  pivot_longer(
    cols = matches(virus_pattern),
    names_to = c("virus_full", "tp"),
    names_pattern = "(log10_.*)norm_(\\d+)$",
    values_to = "log_value"
  ) %>%
  mutate(
    tp = as.integer(tp),
    virus = str_remove(virus_full, "^log10_"),
    age_tp = case_when(
      tp == 0 ~ age_0,
      tp == 1 ~ age_1,
      tp == 2 ~ age_2,
      tp == 3 ~ age_3,
      tp == 4 ~ age_4,
      TRUE ~ NA_real_
    ),
    age_diff = abs(age - age_tp)
  )

# 3) Keep only candidates within ±60 days
candidates <- long_df %>%
  filter(!is.na(log_value), !is.na(age_tp), age_diff <= day_range)

# 4) For each original row and virus, select the closest
closest <- candidates %>%
  group_by(.rowid, virus) %>%
  slice_min(order_by = age_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(.rowid, virus, log_value, sero_age = age_tp)  # add age_sero


# 5) Pivot wider to get one column per virus
wide_df <- closest %>%
  pivot_wider(
    names_from = virus,
    values_from = log_value,
    names_prefix = "log10_"
  )

# 6) Merge back into original dataset
database <- database2 %>%
  left_join(wide_df, by = ".rowid") %>%
  select(-.rowid)











# 1) Detect the virus columns that we actually want to clean:
#    all "log10_" columns that are NOT the original ..._norm_* intensity columns.
sero_vars <- names(database) %>%
  .[grepl("^log10_", .) & !grepl("norm_", .)]

# Check you really have a single serology age column called "sero_age"
# If it has a different name, change this line accordingly.
sero_age_col <- "sero_age"

# 2) Add stable row id
database2 <- database %>%
  mutate(.rowid = dplyr::row_number())

# 3) Loop over each virus column and keep only the closest row per h_id × timepoint_label × sero_age
for (v in sero_vars) {
 
  keep_ids <- database2 %>%
    filter(
      !is.na(.data[[v]]),
      !is.na(.data[[sero_age_col]]),
      !is.na(timepoint_label)
    ) %>%
    mutate(sero_age = .data[[sero_age_col]]) %>%
    group_by(h_id, timepoint_label, sero_age) %>%
    # distance between anthropometry age and serology age
    mutate(dist = abs(age - sero_age)) %>%
    # keep only the closest anthropometry row in this group
    slice_min(order_by = dist, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull(.rowid)
 
  # mask for rows that should be cleared for this virus
  mask <- !(database2$.rowid %in% keep_ids) & !is.na(database2[[v]])
 
  # set virus to NA on all "extra" rows
  database2[[v]][mask] <- NA_real_
}

# 4) Drop helper id and overwrite database
database <- database2 %>%
  select(-.rowid)





library(dplyr)
library(tidyr)
library(stringr)

# ---- 1. Parameters ----
day_range <- 60 / 365.25   # ±60 days in years

# identify all serology markers in serostatus_all_long: log10_*_norm
sero_markers <- names(serostatus_all_long) %>%
  .[grepl("^log10_.*_norm$", .)]

# ---- 2. Build serology-long with age in years ----
sero_long <- serostatus_all_long %>%
  filter(timepoint != 0) %>%                             # ignore timepoint 0 as you did
  select(h_id, age_days, all_of(sero_markers)) %>%       # adjust age_days name if needed
  mutate(age_sero = age_days / 365.25) %>%
  pivot_longer(
    cols = all_of(sero_markers),
    names_to   = "marker_norm",
    values_to  = "value"
  ) %>%
  filter(!is.na(value))                                  # keep only non-missing serology

# marker name in database (without "_norm")
sero_long <- sero_long %>%
  mutate(marker_db = sub("_norm$", "", marker_norm))

# ---- 3. Anthro ages available in database ----
anthro_times <- database %>%
  filter(!is.na(age)) %>%
  select(h_id, age) %>%
  distinct()

# ---- 4. Find serology measurements that have ANY database row within ±60 days ----
# join per h_id, compute age diff, keep matches within range
cand_matches <- sero_long %>%
  inner_join(anthro_times, by = "h_id") %>%
  mutate(age_diff = abs(age_sero - age)) %>%
  filter(age_diff <= day_range)

# each serology record is uniquely (h_id, age_sero, marker_db, value)
matched_keys <- cand_matches %>%
  distinct(h_id, age_sero, marker_db, value)

# ---- 5. Serology measurements with NO matching anthropometry ----
unmatched <- anti_join(
  sero_long,
  matched_keys,
  by = c("h_id", "age_sero", "marker_db", "value")
)

# If there are none, you're done
# if (nrow(unmatched) == 0) { database_expanded <- database }
library(dplyr)
library(tidyr)

# ---- 6. Collapse unmatched to one row per h_id × age_sero (multiple viruses in same row) ----
# First, ensure one numeric value per (h_id, age_sero, marker_db)
unmatched_wide <- unmatched %>%
  group_by(h_id, age_sero, marker_db) %>%
  summarise(value = first(na.omit(value)), .groups = "drop") %>%
  pivot_wider(
    names_from  = marker_db,
    values_from = value
  )

# ---- 7. Create "blank" rows with same structure as database ----
n_new <- nrow(unmatched_wide)

if (n_new > 0) {
  # empty frame with same columns
  template <- database[0, , drop = FALSE]
 
  # create NA-filled df with same columns
  new_rows <- as.data.frame(
    matrix(NA, nrow = n_new, ncol = ncol(template)),
    stringsAsFactors = FALSE
  )
  names(new_rows) <- names(template)
 
  # fill key fields
  new_rows$h_id <- unmatched_wide$h_id
  new_rows$age  <- unmatched_wide$age_sero
 
  # copy over virus columns (and any other matching names)
  for (col in names(unmatched_wide)) {
    if (col %in% names(new_rows)) {
      new_rows[[col]] <- unmatched_wide[[col]]
    }
  }
 
  # ---- 8. Bind new serology-only rows onto the database ----
  database_expanded <- bind_rows(database, new_rows)
} else {
  database_expanded <- database
}




fixed_vars <- c(
  "sex.x", "coh", "pre_bmi_c", "urb_area_id", "smk_p",
  "ethn3_m", "edu_m_0", "m_age", "breastfed_ever", "birth_head_circum", "birth_weight" ,  
  "birth_length" ,  "preg_smk", "nursery_upto2years", "sex" ,"parity_m"
)


lookup <- database_expanded %>%
  # keep rows where at least one fixed var is non-missing
  filter(if_any(all_of(fixed_vars), ~ !is.na(.))) %>%
  group_by(h_id) %>%
  slice(1) %>%          # one representative row per child
  ungroup() %>%
  select(h_id, all_of(fixed_vars))


database_expanded <- database_expanded %>%
  left_join(lookup, by = "h_id", suffix = c("", ".ref"))

for (v in fixed_vars) {
  ref <- paste0(v, ".ref")
  database_expanded[[v]] <- coalesce(database_expanded[[v]], database_expanded[[ref]])
}

database_expanded <- database_expanded %>%
  select(-ends_with(".ref"))


database_bckp <- database
database <- database_expanded












# Correlation


cor_matrix <- cor(select_if(database, is.numeric), use = "pairwise.complete.obs")


 
cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column(var = "var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "cor") %>%
  filter(var1 != var2) %>%              # exclude self-correlations
  mutate(abs_cor = abs(cor)) %>%
  filter(abs_cor >= 0.5) %>%            # set threshold for correlation
  arrange(desc(abs_cor))

cor_long

vars_to_impute <- c("cheight", "cweight", "cabdo", "csubscap", "ctriceps")


cor_selected <- cor_long %>%
  filter(abs_cor > 0.5 & (var1 %in% vars_to_impute | var2 %in% vars_to_impute))



# Extract unique variable names from these selected pairs
selected_vars <- unique(c(cor_selected$var1, cor_selected$var2))

 


 
 
 
 
 
  # --- 0. Prepare data ---
 
  # Assume obesity_all_long is your original dataset with missing values
  # Convert h_id to factor for modeling later
  database$h_id <- as.factor(database$h_id)
 
 
 
  # Identify column positions in `database`
  start_col <- match("timepoint_0", names(database))
  end_col   <- match("cut_Toxo_4", names(database))
 
  # Optional sanity check
  start_col
  end_col
  names(database)[start_col:end_col]
 
  # Drop the entire block of columns
  database <- database[ , -c(start_col:end_col), drop = FALSE]
 
 
 
  # Generate an integer version of h_id that preserves unique IDs (required by 2l.pan)
  imp_data <- database
  imp_data$h_id_int <- as.integer(factor(database$h_id))
 
 
 
 
 
 
  # map the connections for further reference if needed
  mapping <- data.frame(
    original_id = levels(database$h_id),
    int_id = as.integer(factor(levels(database$h_id)))
   
  )
 
 
 
  # --- 1. Variables to impute ---
  # Only raw continuous measurements
   vars_to_impute <- c("cheight", "cweight", "cabdo", "csubscap","ctriceps", "pre_bmi_c",
                       "urb_area_id", "smk_p", "edu_m_0", "ethn3_m", "m_age", "breastfed_ever", "preg_smk", "nursery_upto2years")
 
 
 
  # --- 2. Predictor variables ---
  # Exclude serology to avoid bias; include baseline/exogenous variables
  # Also include other anthropometrics to improve imputation
 
 
  predicting_vars <-  c("age", "weight_zscore", "height_zscore", "bmi_zscore", "coh",  "sex",  "cbmi", "cheight", "cweight", "cabdo", "csubscap","ctriceps", "pre_bmi_c",
                        "urb_area_id", "smk_p", "edu_m_0", "ethn3_m", "m_age", "breastfed_ever", "preg_smk", "nursery_upto2years")  
 
 
  sero_log_cols <- names(database)[
    str_detect(names(database), "^log10_.*_$") |
      names(database) == "log10_EBV_EAD"
  ]
 
 
 
 
  predicting_vars <- c(predicting_vars, sero_log_cols)
 
 
  # --- 3. Setup method vector ---
   meth <- make.method(imp_data)
   meth[vars_to_impute] <- "2l.pan"  # multilevel continuous imputation
   
   
   # All columns except vars_to_impute and predicting_vars
   derived_vars <- setdiff(colnames(imp_data), c(vars_to_impute, predicting_vars, "h_id", "h_id_int"))
   
   # Exclude these from imputation
   meth[derived_vars] <- ""
   
   
   # does the predictors take in account only the individual(h_id) values #sos question
   
   
   
 
 
   
   
   
   
   # --- 4. Build predictor matrix ---
   pred <- make.predictorMatrix(imp_data)
   pred[,] <- 0  # initialize
   
 
 
  # Include chosen predictors for each variable to impute
  for (var in vars_to_impute) {
    pred[var, predicting_vars] <- 1
  }
 
 
  # Random effects: clustering by household
  pred[vars_to_impute, "h_id_int"] <- -2
 
 
  # Longitudinal structure: timepoint as fixed effect
  pred[vars_to_impute, "age"] <- 1
 
  pred[vars_to_impute, sero_log_cols] <- 1
 
 
  # Prevent self-prediction
  diag(pred) <- 0
 
 
  # Exclude derived variables as predictors
  pred[, derived_vars] <- 0
 
 
 
  cat("\nPredictors for each target var (1 = used):\n")
  print(pred[vars_to_impute, c(predicting_vars, "h_id_int", "age")])
 
 
  # --- 5. Run multiple imputation ---
  # Increase number of imputations and iterations for stability
 
 
  # methods etc already defined:
  # meth, pred, imp_data, bounds_list
 
  post <- make.post(imp_data)
 
  post["cweight"] <- "imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], 1.5), 90)"
  post["cheight"] <- "imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], 30), 175)"
 
  post["weight_zscore"] <- "imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], -5), 5)"
  post["height_zscore"] <- "imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], -5), 5)"
  post["bmi_zscore"]    <- "imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], -5), 5)"
 
  set.seed(123)
  imputed <- mice(
    imp_data,
    m              = 20,
    method         = meth,
    predictorMatrix= pred,
    maxit          = 50,
    post           = post
  )
 
 
  imputed_bckp <- imputed
 
 
 
 
  ######### manipulation #########
  completed_data <- complete(imputed, "long") %>%
    group_by(.imp) %>%
    mutate(bmi_scaled = (cbmi - min(cbmi, na.rm = TRUE)) / (max(cbmi, na.rm = TRUE) - min(cbmi, na.rm = TRUE))) %>%
    ungroup()
 
 
  # extract all imputations in long form (including original .imp = 0)
  imp_long <- complete(imputed, action = "long", include = TRUE)
 
  # compute BMI per imputation
  imp_long <- imp_long %>%
    group_by(.imp) %>%
    mutate(cbmi = cweight / (cheight/100)^2) %>%   # if you need recalculation
    ungroup()
 
  # convert back to mids
  imputed <- as.mids(imp_long)
 
 
 
 
 
 
 
  imputed_long <- complete(imputed, action = "long", include = TRUE)
 
  # Compute bmi_scaled per imputation (including .imp = 0 to preserve structure)
  imputed_long <- imputed_long %>%
    group_by(.imp) %>%
    mutate(bmi_scaled = (cbmi - min(cbmi, na.rm = TRUE)) /
             (max(cbmi, na.rm = TRUE) - min(cbmi, na.rm = TRUE))) %>%
    ungroup()
 
  # Convert back to mids object safely
  imputed <- as.mids(imputed_long)




# Trace plot for weight zscore
plot(imputed, c("weight_zscore"))

# Trace plot for height zscore
plot(imputed, c("height_zscore"))

plot

# Trace plot for bmi zscore
plot(imputed, c("bmi_zscore"))

# Stripplot for height
stripplot(imputed, weight_zscore ~ .imp, pch = 20, cex = 1.2)

# Density plot
densityplot(imputed, ~ weight_zscore)
densityplot(imputed, ~ height_zscore)
densityplot(imputed, ~ bmi_zscore)
densityplot(imputed, ~ bmi_scaled)
densityplot(imputed, ~ cbmi)
densityplot(imputed, ~ cheight)
densityplot(imputed, ~ cweight)
#########  #########



# add z scroe and check if imputation zscore is same as normal zscore(weight, height) done


# imputation for serostatus bib timepoint 1-2 , rhea inma timepoint 2+






# Extract each completed dataset from the mids object and tag with .imp
imputed_long <- bind_rows(
  lapply(1:imputed$m, function(i) {
    complete(imputed, i) %>%
      mutate(.imp = i)
  })
)


ggplot(imputed_long, aes(x = weight_zscore, color = factor(.imp))) +
  geom_density(size = 1, alpha = 0.1) +
  labs(color = "Imputation") +
  theme_minimal()

ggplot(imputed_long, aes(x = height_zscore, color = factor(.imp))) +
  geom_density(size = 1, alpha = 0.1) +
  labs(color = "Imputation") +
  theme_minimal()

ggplot(imputed_long, aes(x = bmi_zscore, color = factor(.imp))) +
  geom_density(size = 1, alpha = 0.1) +
  labs(color = "Imputation") +
  theme_minimal()

