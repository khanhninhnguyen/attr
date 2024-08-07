#' Example of dataset
#'
#' @format
#'  A data frame with 53940 rows and 10 variables
#' A list of length n where n is the number of nearby stations.
#' Each element of the list is a data frame representing a main-nearby station
#' pair and is named after the nearby station. contains seven
#' columns:
#'
#' \describe{
#'   \item{date}{dates of measurement in Date format (\%Y-\%m-\%d)}
#'   \item{GE}{GPS - ERA}
#'   \item{GGp}{GPS - GPS'}
#'  `\item{GEp}{GPS - ERA'}
#'  `\item{EEp}{ERA - ERA'}
#'   \item{GpEp}{GPS' - ERA'}
#'   \item{GpE}{GPS' - ERA}
#'   ...
#' }
#'
#' @docType data
#' @keywords datasets
#' @name dataset
#' @usage data(dataset)
"dataset"

