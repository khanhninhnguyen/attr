#' Check if there is a cluster in segmentation result
#'
#' @description
#' This function determine if any clusters of change-points exist in a station
#' based on segmentation results. A cluster is defined as a group of change-points
#' within 80 days.
#'
#' @param Series_df A data frame containing data for a main-nearby pair. It contains
#' at least two columns: the first labeled \code{Date} with dates in "YYYY-MM-DD" format
#' and the second being one of the six difference series, including: GPS - ERA, GPS - GPS',
#' GPS - ERA', ERA - ERA', GPS' - ERA', and GPS' - ERA, with names \code{GE, GGp,
#' GEp, EEp, GpEp,} and \code{GpE}, respectively.
#'
#' @param Name_series The name of the column in \code{Series_df} to check.
#'
#' @param CP A vector of dates in Date format: "%Y-%m-%d" representing
#' change-points in the series to check.
#'
#' @return Vector of break points in a cluster in Date format: "\%Y-\%m-\%d"
#'
#' @export
#'
check_cluster_CP <- function(Series_df, Name_series, CP) {

  if (length(CP) > 0 ) {
    begin_end = get_min_max_date(Series_df, Name_series)
    List_breaks <- sort(c(begin_end, CP))

    get_cluster <- function(X, thres = 80){
      output <- integer(length(X))
      counter <- 1 # Initialize a counter for numbering elements < 80

      for(i in 1:length(X)) {
        if(X[i] > thres) {
          output[i] <- 0
          counter <- 1 # Reset counter for the next cluster of elements < 80
        } else {
          output[i] <- counter
          counter <- counter + 1 # Increment counter for the next element in the current cluster
        }
      }

      return(output)
    }

    # compute length of segments

    length_seg <- numeric(length(List_breaks) - 1)
    length_seg[1] <- sum(Series_df$Date >= List_breaks[1] &
                           Series_df$Date <= List_breaks[2])

    for (i in 2:(length(List_breaks) - 1)) {
      length_seg[i] <- sum(Series_df$Date > List_breaks[i] &
                             Series_df$Date <= List_breaks[i + 1])
    }

    check_cluster = get_cluster(length_seg)

    clust = rep(0, length(List_breaks))

    if (any(check_cluster != 0)) {

      list_seg = which(check_cluster != 0)

      for (i in list_seg) {
        l = rep(check_cluster[i], 2)
        l[2] = l[2] + 1
        clust[i:(i+1)] <- l
      }
    }

    Res <- data.frame(breakpoints = CP,
                      cluster = clust[c(-1, -length(clust))])
  } else {
    Res <- data.frame(cluster = 0)
  }

  return(Res)

}






