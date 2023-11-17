#' Fit a zero inflated negative binomial linear mixed model with covariate correction
#'
#' This function fits a zero inflated negative binomial linear mixed model with user specified Time variable, subjects as random effects, and covariates. The returned value are coutns per time point that have had the covariate effect removed from them. We recommend users to perform ZINBLMM model fitting and selection to determine the appropriate model and covariance structure before using this step. The function will fit a zero-inflated negative binomial linear mixed model with the covariates to correct by. The function will first convert the count table and sample data into a long formatted model table. The model will be fit with the outcome variable as "asvcounts" and the predictors as covariates to correct by with subjects modeled as random effects.
#'
#' @param Obj (Required) A LIMON Object
#' @param Time (Required) A column in "SampleData" for the Time variable. Must be numeric ordinal (1,2,3, etc)
#' @param Subject (Required) A column of subject IDs. Used as
#' @param Covariates (Required) Specified string of columns from SampleData to use as covariates in the model
#' @param model (Required) Specify the right side of the lmer equation (the covariates to correct by)
#' @param prop.data (Optional) Option to specify if these are proprotional data or not. [Default] = FALSE indicating these are raw counts. If using relative abundance or other proportion data, specify TRUE.
#' @param correlation (Optional) Option to specify the correlation structure for model fitting ([Default] of NULL specifies unstructured). Follows the same covariance structure names as in NBZIMM
#'
#' @returns LIMON Object with the addition; (1) Model_data: contains the Counts and SampleData in Long Format for the model. (2) equation: the formula for ZINBLMM. (3) Model_res: The results from running glmm.zinb(). (4) Corrected_Counts: corrected count values with samples as rows and taxa as columns. (5) Corrected_Counts_Time: corrected count values per time point with samples as rows and taxa as columns.
#' @export
#'
LIMON_DistrFit <- function(Obj, Time, Subject, Covariates, model, prop.data = FALSE,
                           correlation = NULL){

  # Required packages
  library(NBZIMM)

  # Step 1: Make data frame for running the model
  ##################################################

  # Create an empty list to store replicated covariates
  rep_cov <- list()
  # Loop to replicate covariate columns and store them in the list
  for (col_name in Covariates) {
    original_col <- Obj[["SampleData"]][[col_name]]
    replicated_col <- rep(original_col, times = ncol(Obj[["Counts"]]))
    rep_cov[[col_name]] <- replicated_col
  }

  # Replicate the sampleid, timepoint, and asvcounts
  timepoint <- rep(Obj[["SampleData"]][[Time]], ncol(Obj[["Counts"]]))
  subjectid <- rep(Obj[["SampleData"]][[Subject]], ncol(Obj[["Counts"]]))
  asvcounts <- unlist(Obj[["Counts"]])

  # Merge into one data frame & add to the Object
  model_data <- data.frame(subjectid, timepoint, asvcounts, rep_cov)
  Obj[["Model_data"]] <- model_data

  # Step 2: Run the appropriate model based on prop.data
  ##################################################

  # Generate the linear model equation
  eq <- formula(paste("asvcounts", "~", model))

  if (prop.data) {
    # Run lme.zig() when prop.data is TRUE
    zinblmm_res <- lme.zig(fixed = eq,
                           random = ~ 1 | subjectid,
                           data = model_data,
                           correlation = correlation,
                           zi_fixed = ~ 1)

    # Add model results to the object
    Obj[["equation"]] <- eq
    Obj[["Model_res"]] <- zinblmm_res
  }
  else {
    # Run glmm.zinb() when prop.data is FALSE
    zinblmm_res <- glmm.zinb(fixed = eq,
                             offset = log(sum("asvcounts")),
                             random = ~ 1 | subjectid,
                             data = model_data,
                             correlation = correlation,
                             family = zi.negative.binomial(),
                             zi_fixed = ~ 1)

    # Add model results to the object
    Obj[["equation"]] <- eq
    Obj[["Model_res"]] <- zinblmm_res
  }

  # Step 3: Extract the corrected counts and give as output
  ########################################################

  # Get corrected counts
  data_residuals <- resid(zinblmm_res)
  data_fitted <- fitted(zinblmm_res)


  lmm_corrected_predictors <- data_fitted - data_residuals
  lmmcor_pred_data <- data.frame(lmm_corrected_predictors)
  lmmcor_pred_data <- data.frame(split(lmmcor_pred_data,
                                       rep(1:ncol(Obj[["Counts"]]),
                                           each=nrow(Obj[["Counts"]]))))
  names(lmmcor_pred_data) <- names(Obj[["Counts"]])
  rownames(lmmcor_pred_data) <- rownames(Obj[["Counts"]])

  # Assign final all corrected counts to object
  Obj[["Corrected_Counts"]] <- lmmcor_pred_data

  # Step 4: Separate the corrected counts by visit
  ########################################################
  # Initialize an empty list to store count tables for each timepoint
  Corrected_Counts_Timepoint <- list()

  # Get unique timepoints from the sample data
  sample_data <- Obj[["SampleData"]]
  timepoints <- unique(Obj[["SampleData"]][[Time]])

  # Loop through each unique time point and create a separate count table
  for (tp in timepoints) {
    # Subset sample data and count table for the current timepoint
    subset_sample_data <- sample_data[sample_data[[Time]] == tp, ]
    timept_counts <- lmmcor_pred_data[row.names(lmmcor_pred_data) %in% row.names(subset_sample_data), ]
    # Time point specific count tables to the list
    timept_counts_name <- paste0("Corrected_Counts_", tp)
    Corrected_Counts_Timepoint[[timept_counts_name]] <- timept_counts
  }

  # Store the list of count tables in the original object
  Obj[["Corrected_Counts_Time"]] <- Corrected_Counts_Timepoint
  Obj
}
