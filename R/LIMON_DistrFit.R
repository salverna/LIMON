#' Fit a zero inflated negative binomial linear mixed model with covariate correction
#'
#' This function fits a zero inflated negative binomial linear mixed model with user specified Time variable, subjects as random effects, and covariates. The returned value are coutns per time point that have had the covariate effect removed from them. We recommend users to perform ZINBLMM model fitting and selection to determine the appropriate model and covariance structure before using this step. The function will fit a zero-inflated negative binomial linear mixed model with the covariates to correct by. The function will first convert the count table and sample data into a long formatted model table. The model will be fit with the outcome variable as "asvcounts" and the predictors as covariates to correct by with subjects modeled as random effects.
#'
#' @param Obj (Required) A LIMON Object
#' @param Time (Required) A column in "SampleData" for the Time variable. Must be numeric ordinal (1,2,3, etc)
#' @param Subject (Required) A column of subject IDs. Used as
#' @param Covariates (Required) Specified string of columns from SampleData to use as covariates in the model
#' @param model (Required) Specify the right side of the lmer equation (the covariates to correct by)
#' @param prop.data (Optional) Option to specify if these are proportional data or not. [Default] = FALSE indicating these are raw counts. If using relative abundance or other proportion data, specify TRUE. lme.zig() equation from the NBZIMM package.
#' @param correlation (Optional) Option to specify the correlation structure for model fitting ([Default] of NULL specifies unstructured). Follows the same covariance structure names as in NBZIMM
#'@param use_zinb (Optional) Species if zero inflated data. [Default] = TRUE. If true and prop.data is FALSE, then the model will run a Zero-Inflated Negative Binomial LMM, glmm.zinb() function in NBZIMM. If this is set to FALSE, it will run a Negative Binomial LMM, glmm.nb() in NBZIMM.
#'
#' @returns LIMON Object with the addition; (1) Model_data: contains the Counts and SampleData in Long Format for the model. (2) equation: the formula for ZINBLMM. (3) Model_res: The results from running the linear mixed model. (4) Corrected_Counts: corrected count values with samples as rows and taxa as columns. (5) Corrected_Counts_Time: corrected count values per time point with samples as rows and taxa as columns.
#' @export
#'
LIMON_DistrFit <- function(Obj, Time, Subject, Covariates, model, prop.data = FALSE, correlation = NULL, use_zinb = TRUE) {
  library(NBZIMM)

  # Check if row names are identical and in the same order in Counts and SampleData
  if (!identical(rownames(Obj$Counts), rownames(Obj$SampleData))) {
    stop("Row names of Counts and SampleData tables are not identical. Make sure the names and order match between the two")
  }

  # Initialize a list to store corrected counts for each column
  corrected_counts_list <- base::list()

  for (col_name in colnames(Obj$Counts)) {
    # Prep        are the data for the current column
    model_data <- cbind(Obj$SampleData[, c(Time, Subject)], Obj$Counts[, col_name])
    colnames(model_data) <- c("timepoint", "subjectid", "asvcounts")

    for (covariate in Covariates) {
      model_data[[covariate]] <- Obj$SampleData[[covariate]]
    }

    # Fit the model
    eq <- formula(base::paste("asvcounts", "~", model))

    # Choose the appropriate function based on user's choice
    if (use_zinb) {
      if (prop.data) {
        model_res <- NBZIMM::lme.zig(fixed = eq, random = ~1 | subjectid,
                                     data = model_data, correlation = correlation, zi_fixed = ~1)
      } else {
        model_res <- NBZIMM::glmm.zinb(fixed = eq, offset = log(sum(asvcounts)),
                                       random = ~1 | subjectid, data = model_data, correlation = correlation,
                                       family = zi.negative.binomial(), zi_fixed = ~1)
      }
    } else {
      model_res <- NBZIMM::glmm.nb(fixed = eq, offset = log(sum(asvcounts)),
                                   random = ~1 | subjectid, data = model_data, correlation = correlation)
    }

    Obj[["equation"]] <- eq
    Obj[["Model_res"]] <- model_res

    # Store corrected counts
    corrected_counts_list[[col_name]] <- resid(model_res)

  }
  Corrected_Counts <- as.data.frame(corrected_counts_list)


  # Remove "label" attribute from each column
  for (col in names(Corrected_Counts)) {
    base::attr(Corrected_Counts[[col]], "label") <- NULL
  }

  # Set column names and row names
  colnames(Corrected_Counts) <- colnames(Obj[["Counts"]])
  rownames(Corrected_Counts) <- rownames(Obj[["Counts"]])

  # Add one to scale the data for the appropriate distribution and save the object
  #Corrected_Counts <- Corrected_Counts + 1

  # scale to positive and make larger
  Corrected_Counts <- Corrected_Counts + abs(min(Corrected_Counts))
  Corrected_Counts <- round(Corrected_Counts*10)

  Obj$Corrected_Counts <- Corrected_Counts

  # Store per timepoint
  timepoints <- unique(Obj[["SampleData"]][[Time]])

  Corrected_Counts_Timepoint <- base::list()

  for (tp in timepoints) {
    subset_sample_data <- Obj$SampleData[Obj$SampleData[[Time]] == tp, ]
    timept_counts <- Obj$Corrected_Counts[row.names(Corrected_Counts) %in% row.names(subset_sample_data), ]
    timept_counts_name <- base::paste0("Corrected_Counts_", tp)
    Corrected_Counts_Timepoint[[timept_counts_name]] <- timept_counts
  }
  Obj$Corrected_Counts_Time <- Corrected_Counts_Timepoint

  return(Obj)
}
