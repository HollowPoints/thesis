md.pattern(complete_data[,1:30])

pred <- make.predictorMatrix(complete_data)

impute_vars <- c("cheight", "cweight", "cabdo", "csubscap", "ctriceps", "bmi_zscore",
                 "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore",
                 "whtr", "central_obesity", "stunted", "cbmi")



pred <- make.predictorMatrix(complete_data)

# Set all entries to 0 to exclude all variables from prediction.
pred[,] <- 0

# For each variable you want to impute (the rows in the matrix)...
for (i in impute_vars) {
  # ...set its predictors to be the other variables from the same list.
  pred[i, impute_vars] <- 1
}

# A variable cannot predict itself, so set the diagonal to 0.
diag(pred) <- 0

# Create a methods vector for all columns in your data, defaulting to ""
meth <- make.method(complete_data)

# Specify methods for your continuous variables.
meth[c("cheight", "cweight", "cabdo", "csubscap", "ctriceps", "bmi_zscore",
       "weight_zscore", "height_zscore", "tris_zscore", "subscap_zscore",
       "whtr", "cbmi")] <- "pmm"

# Specify methods for your binary/categorical variables.
meth[c("central_obesity", "stunted")] <- "logreg"



imputed_data <- mice(complete_data, m = 5, method = meth, predictorMatrix = pred, maxit = 20, seed = 123)



