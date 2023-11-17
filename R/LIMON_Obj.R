#' Counts to LIMON
#'
#'Create Formatted Objects to run LIMON from raw count tables
#'
#' @param Counts A genomic count table with samples as rows and taxa/genes/etc as columns
#' @param SampleData A sample data frame with row names matching "Counts"
#' @returns Returns a LIMON Object with the Counts and SampleData stored
#' @export
#'

LIMON_Obj <- function(Counts, SampleData){
  list(Counts = Counts, SampleData = SampleData)}

