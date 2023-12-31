% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NEC_Participant_data}
\alias{NEC_Participant_data}
\title{NEC subset of Participant Data from "Necrotizing enterocolitis is preceded by increased gut bacterial replication, Klebsiella, and fimbriae-encoding bacteria." Olm et al. Sci Adv 2019;5(12):eaax5727}
\format{
A data frame with 150 rows and 23 variables:
\describe{
  \item{\code{Participant_ID}}{character COLUMN_DESCRIPTION}
  \item{\code{Age.at.study..days...EUPATH_0009212.}}{integer COLUMN_DESCRIPTION}
  \item{\code{Stunted.at.birth..EUPATH_0009216.}}{character COLUMN_DESCRIPTION}
  \item{\code{Weight..g...IAO_0000414.}}{integer COLUMN_DESCRIPTION}
  \item{\code{Sepsis.diagnosed..EUPATH_0009210.}}{character COLUMN_DESCRIPTION}
  \item{\code{Sex..PATO_0000047.}}{character COLUMN_DESCRIPTION}
  \item{\code{Contents.of.diet..EUPATH_0050531.}}{character COLUMN_DESCRIPTION}
  \item{\code{Country..ENVO_00000009.}}{character COLUMN_DESCRIPTION}
  \item{\code{Maternal.antepartum.antibiotics.administered..EUPATH_0009191.}}{character COLUMN_DESCRIPTION}
  \item{\code{Maternal.antepartum.antibiotics.administered.time..EUPATH_0009192.}}{character COLUMN_DESCRIPTION}
  \item{\code{Chorioamnionitis..EUPATH_0009175.}}{character COLUMN_DESCRIPTION}
  \item{\code{Group.B.Streptococcus.infection..EUPATH_0009172.}}{character COLUMN_DESCRIPTION}
  \item{\code{Necrotizing.enterocolitis..EUPATH_0009177.}}{character COLUMN_DESCRIPTION}
  \item{\code{Gestational.age.at.birth..days...EUPATH_0009095.}}{integer COLUMN_DESCRIPTION}
  \item{\code{Delivery.mode..EUPATH_0000581.}}{character COLUMN_DESCRIPTION}
  \item{\code{Amnionicity..EUPATH_0009187.}}{character COLUMN_DESCRIPTION}
  \item{\code{Birth.plurality..EUPATH_0009186.}}{character COLUMN_DESCRIPTION}
  \item{\code{Chorionicity..EUPATH_0009188.}}{character COLUMN_DESCRIPTION}
  \item{\code{Twin.type..EUPATH_0009189.}}{character COLUMN_DESCRIPTION}
  \item{\code{Season.of.birth..EUPATH_0009209.}}{character COLUMN_DESCRIPTION}
  \item{\code{Sibling.set.ID..EUPATH_0009215.}}{integer COLUMN_DESCRIPTION}
  \item{\code{Study.group..EUPATH_0000639.}}{character COLUMN_DESCRIPTION}
  \item{\code{Survived..EUPATH_0009211.}}{character COLUMN_DESCRIPTION}
}
}
\source{
\url{https://microbiomedb.org/mbio/app/workspace/analyses/DS_84fcb69f4e/new/details}
}
\usage{
NEC_Participant_data
}
\description{
This was a NICU study that conducted deep metagenomic profiling of samples collected longitudinally from 150 pre-term infants. 30 of these infants developed necrotizing enterocolitis. For this example data, we used the MicrobiomeDB to pull the relative abundance metagenomic tables at the species level.
}
\keyword{datasets}
