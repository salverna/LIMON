#' Individual Network Inferences
#'
#' # Function to perform network inferences following method in lionessR - https://github.com/mararie/lionessR. Note: Consider using the same methods and parameters as in LIMON_NetInf_Time(). See the \pkg{SpiecEasi} documentation for more details.
#'
#' @param Obj LIMON Object from LIMON_NetInf_Time()
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
#' @return LIMON object with individual networks per sample stored in Individual_Networks by sample name
#' @export
#'
LIMON_IndNet <- function(Obj, method,
                         sel.criterion,
                         lambda.min.ratio,
                         nlambda ,
                         pulsar.select=TRUE,
                         pulsar.params=list(),
                         icov.select = pulsar.select,
                         icov.select.params = pulsar.params,
                         lambda.log = TRUE){

  #Make and individual networks list
  Individual_Networks <- list()

  for (i in seq_along(Obj[["SpiecEasi_Time"]])){

    #Call the Network per time point
    net_all <- as.matrix(getOptCov(Obj[["SpiecEasi_Time"]][[i]]))

    #Get count table for time point i
    count_table <- Obj[["Corrected_Counts_Time"]][[i]]

    #Compute network minus each individual for time point i
    for (j in 1:nrow(count_table)){

      #Get sample name
      samplename <- rownames(count_table[j,])

      #Get number of samples
      nsamples <- nrow(count_table)

      #Apply Spiec Easi to each all the samples minus one
      net_single <- spiec.easi(data.matrix(count_table[-j,]),
                               data.matrix(count_table),
                               method = method,
                               sel.criterion = sel.criterion,
                               lambda.min.ratio = lambda.min.ratio,
                               nlambda = nlambda,
                               pulsar.select=pulsar.select,
                               pulsar.params=pulsar.params,
                               icov.select = icov.select,
                               icov.select.params = icov.select.params,
                               lambda.log = lambda.log)

      net_single_cov <- as.matrix(getOptCov(net_single))

      net_name <- samplename

      Individual_Networks[[net_name]] <- nsamples*(net_all - net_single_cov) + net_single_cov
      colnames(Individual_Networks[[net_name]]) <- rownames(Individual_Networks[[net_name]]) <- colnames(count_table)
    }

  }
  Obj[["Individual_Networks"]] <- Individual_Networks
  Obj
}
