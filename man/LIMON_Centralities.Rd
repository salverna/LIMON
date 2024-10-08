% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LIMON_Centralities.R
\name{LIMON_Centralities}
\alias{LIMON_Centralities}
\title{Individual network characteristics}
\usage{
LIMON_Centralities(Obj, threshold = 0.2)
}
\arguments{
\item{Obj}{LIMON Object from LIMON_NetInf_Time()}

\item{threshold}{(Optional) Same threshold parameter used in the LIMON_IndEdges step. [Default] = 0.2}
}
\value{
LIMON Object with updated SampleData that contains network centralities per sample. SampleData is written to Global environment as LIMON_SampleData. Note the number of random walks for the Communities is 30
}
\description{
Calculate network centralities from each individual graph and add to Sample Data from export. Returns average Degree Centrality, average Closeness Centrality, average Betweenness Centrality, Eigenvector Centrality, and total Communities per sample.
}
