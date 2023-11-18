#' Individual network characteristics
#'
#' Calculate network centralities from each individual graph and add to Sample Data from export. Returns average Degree Centrality, average Closeness Centrality, average Betweenness Centrality, Eigenvector Centrality, and total Communities per sample.
#'
#' @param Obj LIMON Object from LIMON_NetInf_Time()
#'
#' @return LIMON Object with updated SampleData that contains network centralities per sample. SampleData is written to Global environment as LIMON_SampleData
#' @export
#'
LIMON_Centralities <- function(Obj){
  SampleData <- Obj$SampleData
  Individual_Networks <- Obj$Individual_Networks

  #vectors to store centralities
  degree_values <- numeric(length(rownames(SampleData)))
  closeness_values <- numeric(length(rownames(SampleData)))
  betweenness_values <- numeric(length(rownames(SampleData)))
  eigenvector_values <- numeric(length(rownames(SampleData)))
  communities <- numeric(length(rownames(SampleData)))

  for (i in 1:nrow(SampleData)) {
    # Get the rowname from SampleData
    sample_id <- rownames(SampleData)[i]

    # Get the corresponding igraph object from Individual_Networks
    graph <- graph.adjacency(Individual_Networks[[sample_id]], mode = "undirected", weighted = NULL)
    # Degree Centrality
    degree_values[i] <- mean(degree(graph), na.rm = TRUE)

    # Closeness Centrality
    closeness_values[i] <- mean(closeness(graph), na.rm = TRUE)

    # Betweenness Centrality
    betweenness_values[i] <- mean(betweenness(graph), na.rm = TRUE)

    # Eigenvector Centrality
    eigenvector_value <- eigen_centrality(graph)
    eigenvector_values[i] <- eigenvector_value[["value"]]

    # Return the total number of communties from a cluster random walk
    communities[i] <- max(membership(cluster_walktrap(graph)))

  }

  #Add centralities to SampleData
  SampleData$DegreeCentrality <- degree_values
  SampleData$ClosenessCentrality <- closeness_values
  SampleData$BetweennessCentrality <- betweenness_values
  SampleData$EigenCentrality <- eigenvector_values
  SampleData$Communitiies <- communities

  Obj$SampleData <- SampleData

  # Write SampleData to the global environment
  assign("LIMON_SampleData", SampleData, envir = .GlobalEnv)

  return(Obj)
}
