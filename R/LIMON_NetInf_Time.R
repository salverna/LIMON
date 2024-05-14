#' Infer Network per time point
#'
#' Uses SPIEC-EASI to fit a network per Time point and return that network object, the optimal covariance matrix per time, and differences in covariance matrix per time. See the \pkg{SpiecEasi} documentation for more details.
#'
#' @param Obj LIMON Object after performing LIMON_DistrFit()
#' @param method Default "glasso". "mb" other option from SpiecEasi
#' @param sel.criterion Method for model selection. Either 'stars' [default] or 'bstars' (Bounded StARS)
#' @param lambda.min.ratio Scaling factor that will set the minimum lambda (sparsity) parameter
#' @param nlambda sparsity parameter (see SpiecEasi for more details)
#' @param pulsar.select Default TRUE
#' @param pulsar.params Default empty list()
#' @param icov.select Default "pulsar.select"
#' @param icov.select.params Default "pulsar.params"
#' @param lambda.log Default TRUE
#'
#' @return Returns a LIMON object with SpiecEasi networks and Covariance matrices by time added to it. SpiecEasi_Time_i: SpiecEasi object per time point. CovMatrix_Time_i: Covariance matrix per time point extracted from Spiec Easi using getOptCov(). CovMatrix_Diff_ji: Difference in covariance matrix between time j and time i.
#' @export
#'
LIMON_NetInf_Time <- function(Obj, method,
                              sel.criterion,
                              lambda.min.ratio,
                              nlambda ,
                              pulsar.select=TRUE,
                              pulsar.params=list(),
                              icov.select = pulsar.select,
                              icov.select.params = pulsar.params,
                              lambda.log = TRUE) {

  # Load required packages
  library(SpiecEasi)

  #empty list to store the Timepoint covariances matrices in
  CovMatrix_Time <- list()
  SpeicEasi_Time <- list()
  CovMatrix_Diff <- list()

  # Iterate through each Corrected_Counts_Time and run SpiecEasi
  for (i in seq_along(Obj[["Corrected_Counts_Time"]])) {
    count_table <- Obj[["Corrected_Counts_Time"]][[i]]

    # Fit a basic Spiec Easi Association
    net_data <- spiec.easi(data.matrix(count_table),
                           method = method,
                           sel.criterion = sel.criterion,
                           lambda.min.ratio = lambda.min.ratio,
                           nlambda = nlambda,
                           pulsar.select=pulsar.select,
                           pulsar.params=pulsar.params,
                           icov.select = icov.select,
                           icov.select.params = icov.select.params,
                           lambda.log = lambda.log)

    # Convert Covariance matrix to objects based on time point number
    network_name <- paste0("Net_", i)
    SpeicEasi_Time[[network_name]] <- net_data

    # Convert Covariance matrix to objects based on time point number
    cov_matrix_name <- paste0("Cov_Matrix_", i)
    CovMatrix_Time[[cov_matrix_name]] <- as.matrix(getOptCov(net_data))
    colnames(CovMatrix_Time[[cov_matrix_name]]) <- rownames(CovMatrix_Time[[cov_matrix_name]]) <-   colnames(count_table)

  }
  Obj[["SpiecEasi_Time"]] <- SpeicEasi_Time
  Obj[["CovMatrix_Time"]] <- CovMatrix_Time
  # Calculate difference in Covariance Matrices between Visits


  # number of differences
  N <- length(Obj[["CovMatrix_Time"]])

  # Loop over all pairs of count tables and calculate differences
  for (i in 1:(N-1)) {
    for (j in (i + 1):N) {
      diff_table <- Obj[["CovMatrix_Time"]][[j]] - Obj[["CovMatrix_Time"]][[i]]
      # Set the name for the difference table
      diff_table_name <- paste0("CovMatrix_Diff_", j, i)
      # Add the difference table to the new list
      CovMatrix_Diff[[diff_table_name]] <- diff_table
    }
  }


  # Store the list of networks and covariance matrices in the og object
  Obj[["CovMatrix_Diff"]] <- CovMatrix_Diff
  Obj
}
