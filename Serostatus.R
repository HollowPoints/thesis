besity_all_long$age <- obesity_all_long$agecd_cgrowth / 365

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


## add cbmi
processed_list1 <- processed_list

processed_list <- lapply(
  processed_list,
  function(df) {
    cbmi_min <- min(df$cbmi, na.rm = TRUE)
    cbmi_max <- max(df$cbmi, na.rm = TRUE)
    denom    <- cbmi_max - cbmi_min
   
    df$bmi_scaled <- if (is.finite(denom) && denom > 0) {
      (df$cbmi - cbmi_min) / denom
    } else {
      NA_real_
    }
   
    df
  }
)

##
processed_list_20 <- lapply(processed_list, function(df) df[1:20, , drop = FALSE])

saveRDS(processed_list, "processed_list.rds")
saveRDS(anthro_long_sero_wide, "anthro_long_sero_wide.rds")
saveRDS(serostatus_all_long, "serostatus_all_long.rds")



# cox model

library(dplyr)
library(purrr)

processed_list <- lapply(processed_list, function(df) {
  df %>%
    mutate(
      bmi_scaled = (cbmi - min(cbmi, na.rm = TRUE)) /
        (max(cbmi, na.rm = TRUE) - min(cbmi, na.rm = TRUE))
    )
})


library(survival)

cox_list_bmi_scaled <- lapply(processed_list, function(df) {
  coxph(Surv(age, CMV_class_sero) ~ bmi_scaled + sex + coh, data = df)
})

m <- length(cox_list_bmi_scaled)

# extract beta and its variance for bmi_scaled from each fit
beta_vec <- sapply(cox_list_bmi_scaled, function(fit) {
  coef(fit)[["bmi_scaled"]]
})

var_vec <- sapply(cox_list_bmi_scaled, function(fit) {
  vcov(fit)[ "bmi_scaled", "bmi_scaled" ]
})

# Rubin's rules
Q_bar <- mean(beta_vec)           # pooled beta
W     <- mean(var_vec)            # within-imputation variance
B     <- var(beta_vec)            # between-imputation variance
T_var <- W + (1 + 1/m) * B        # total variance
se    <- sqrt(T_var)              # pooled SE

# degrees of freedom (Barnard–Rubin)
df <- (m - 1) * (1 + W / ((1 + 1/m) * B))^2

# test statistics
t_val <- Q_bar / se
p_val <- 2 * pt(-abs(t_val), df = df)

# Hazard ratio and CI
HR    <- exp(Q_bar)
lower <- exp(Q_bar - qt(0.975, df) * se)
upper <- exp(Q_bar + qt(0.975, df) * se)

pooled_result_bmi_scaled <- data.frame(
  term        = "bmi_scaled",
  beta        = Q_bar,
  se          = se,
  df          = df,
  HR          = HR,
  CI_lower    = lower,
  CI_upper    = upper,
  t_value     = t_val,
  p_value     = p_val
)

pooled_result_bmi_scaled


library(ggplot2)

ggplot(pooled_result_bmi_scaled,
       aes(x = term, y = HR)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Pooled effect of BMI (scaled) on CMV seropositivity",
    x = "",
    y = "Hazard Ratio (HR)"
  )

hr_each <- sapply(cox_list_bmi_scaled, function(fit) {
  exp(coef(fit)[["bmi_scaled"]])
})

df_hr <- data.frame(
  imp = 1:m,
  HR  = hr_each
)

ggplot(df_hr, aes(x = imp, y = HR)) +
  geom_point() +
  geom_line(alpha = 0.4) +
  geom_hline(yintercept = HR, linetype = "dashed", colour = "red") +
  theme_minimal() +
  labs(
    title = "HR(bmi_scaled) across imputations (CMV_class_sero)",
    x = "Imputation index",
    y = "Hazard Ratio"
  )


hr_each <- sapply(cox_list_bmi_scaled, function(fit) {
  exp(coef(fit)[["bmi_scaled"]])
})

df_hr <- data.frame(
  imp = 1:m,
  HR  = hr_each
)

ggplot(df_hr, aes(x = imp, y = HR)) +
  geom_point() +
  geom_line(alpha = 0.4) +
  geom_hline(yintercept = HR, linetype = "dashed", colour = "red") +
  theme_minimal() +
  labs(
    title = "HR(bmi_scaled) across imputations (CMV_class_sero)",
    x = "Imputation index",
    y = "Hazard Ratio"
  )


## first descriptive statistics imputed vs non-imputed (per cohort)


# maybe start with 4 year exposure, as one timepoint
## add the categorical variable of bmi(normal, obese etc), maybe skinfold, w circumference
## exposure bmi 3-5 age, (older child for duplicates)
## minimally adjusted for cofounders(sex, cohort, maybe exact age of exposure bmi)
## dag adjusted model with all relevant confounders



###

## --- 0. Setup -------------------------------------------------------------

library(dplyr)
library(purrr)
library(ggplot2)

## processed_list: list of 20 imputed + post-processed data frames
## Assumptions:
## - Each df has: h_id, age, bmi_scaled, sex, coh
## - Nine serology outcomes ending in "_sero", coded 0/1/NA


## --- 1. Make sure bmi_scaled exists and is consistent in all datasets -----

processed_list <- lapply(processed_list, function(df) {
  if (!"bmi_scaled" %in% names(df)) {
    df %>%
      mutate(
        bmi_scaled = (cbmi - min(cbmi, na.rm = TRUE)) /
          (max(cbmi, na.rm = TRUE) - min(cbmi, na.rm = TRUE))
      )
  } else {
    df
  }
})


## --- 2. Identify serology outcomes ---------------------------------------

outcome_vars <- names(processed_list[[1]])[grepl("_sero$", names(processed_list[[1]]))]
outcome_vars
# e.g. "CMV_class_sero", "cut_VZV_sero", ..., "cut_MCV_sero"


## --- 3. Fit logistic models for a single outcome across imputations ------

fit_logit_for_outcome <- function(outcome_name, data_list) {
  # For each imputed dataset, fit:
  #   outcome ~ bmi_scaled + sex + coh + age   (logistic)
  lapply(data_list, function(df) {
   
    # Drop rows with missing outcome or key covariates
    df2 <- df %>%
      filter(
        !is.na(.data[[outcome_name]]),
        !is.na(bmi_scaled),
        !is.na(sex),
        !is.na(coh),
        !is.na(age)
      )
   
    # If no variation in outcome, model is not estimable
    if (dplyr::n_distinct(df2[[outcome_name]]) < 2) {
      return(NULL)
    }
   
    # Fit logistic regression
    tryCatch(
      glm(
        formula = as.formula(paste0(outcome_name, " ~ bmi_scaled + sex + coh + age")),
        family  = binomial,
        data    = df2
      ),
      error = function(e) NULL
    )
  })
}


## --- 4. Rubin pooling for bmi_scaled from list of glm models --------------

pool_bmi_scaled_logit <- function(glm_list, outcome_name) {
  # Remove failed / NULL fits
  glm_list <- glm_list[!vapply(glm_list, is.null, logical(1))]
 
  if (length(glm_list) <= 1) {
    return(
      data.frame(
        outcome   = outcome_name,
        term      = "bmi_scaled",
        beta      = NA_real_,
        se        = NA_real_,
        df        = NA_real_,
        OR        = NA_real_,
        CI_lower  = NA_real_,
        CI_upper  = NA_real_,
        z_value   = NA_real_,
        p_value   = NA_real_
      )
    )
  }
 
  # Extract beta and variance of bmi_scaled
  beta_vec <- sapply(glm_list, function(fit) {
    coef(fit)[["bmi_scaled"]]
  })
 
  var_vec <- sapply(glm_list, function(fit) {
    vcov(fit)["bmi_scaled", "bmi_scaled"]
  })
 
  # Drop NAs
  keep <- !is.na(beta_vec) & !is.na(var_vec)
  beta_vec <- beta_vec[keep]
  var_vec  <- var_vec[keep]
  m <- length(beta_vec)
 
  if (m <= 1) {
    return(
      data.frame(
        outcome   = outcome_name,
        term      = "bmi_scaled",
        beta      = NA_real_,
        se        = NA_real_,
        df        = NA_real_,
        OR        = NA_real_,
        CI_lower  = NA_real_,
        CI_upper  = NA_real_,
        z_value   = NA_real_,
        p_value   = NA_real_
      )
    )
  }
 
  # Rubin's rules for logistic regression (same math as for linear/Cox beta)
  Q_bar <- mean(beta_vec)   # pooled beta
  W     <- mean(var_vec)    # within-imputation variance
  B     <- var(beta_vec)    # between-imputation variance
  T_var <- W + (1 + 1/m) * B
  se    <- sqrt(T_var)
 
  # Barnard–Rubin df
  df <- (m - 1) * (1 + W / ((1 + 1/m) * B))^2
 
  # z (or t) statistic and p-value
  z_val <- Q_bar / se
  # With large m, normal approx is fine, but keep df-based t for consistency
  p_val <- 2 * pt(-abs(z_val), df = df)
 
  # Odds ratios and 95% CI
  OR    <- exp(Q_bar)
  lower <- exp(Q_bar - qt(0.975, df) * se)
  upper <- exp(Q_bar + qt(0.975, df) * se)
 
  data.frame(
    outcome   = outcome_name,
    term      = "bmi_scaled",
    beta      = Q_bar,
    se        = se,
    df        = df,
    OR        = OR,
    CI_lower  = lower,
    CI_upper  = upper,
    z_value   = z_val,
    p_value   = p_val
  )
}


## --- 5. Fit and pool across all 9 serology outcomes ----------------------

# 5a. Fit all logistic models per outcome
glm_lists_all <- lapply(outcome_vars, fit_logit_for_outcome, data_list = processed_list)
names(glm_lists_all) <- outcome_vars

# 5b. Pool bmi_scaled effect for each outcome
pooled_results_bmi_scaled <- bind_rows(
  Map(pool_bmi_scaled_logit, glm_lists_all, outcome_vars)
)

pooled_results_bmi_scaled


## --- 6. Forest plot of pooled OR(bmi_scaled) by serology outcome ---------

# Clean nicer labels if you want (optional)
pooled_results_bmi_scaled$outcome_label <- pooled_results_bmi_scaled$outcome

ggplot(
  pooled_results_bmi_scaled %>% filter(!is.na(OR)),
  aes(x = outcome_label, y = OR)
) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Pooled odds ratio of BMI (scaled) for each serology outcome",
    x = "Serology outcome",
    y = "Odds Ratio (per 0–1 BMI scale)"
  )


## --- 7. Optional: show per-imputation OR(bmi_scaled) for each outcome ----

hr_long <- bind_rows(
  lapply(seq_along(outcome_vars), function(j) {
    out_name <- outcome_vars[j]
    glm_list <- glm_lists_all[[j]]
    glm_list <- glm_list[!vapply(glm_list, is.null, logical(1))]
   
    if (length(glm_list) == 0) {
      return(data.frame(
        outcome = character(0),
        imp     = integer(0),
        OR      = numeric(0)
      ))
    }
   
    OR_each <- sapply(glm_list, function(fit) {
      exp(coef(fit)[["bmi_scaled"]])
    })
   
    data.frame(
      outcome = out_name,
      imp     = seq_along(OR_each),
      OR      = OR_each
    )
  })
)

ggplot(hr_long %>% filter(!is.na(OR)),
       aes(x = imp, y = OR, group = outcome, colour = outcome)) +
  geom_line(alpha = 0.4) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "OR(bmi_scaled) across imputations for each serology outcome",
    x = "Imputation index",
    y = "Odds Ratio"
  )

 







##
