#' Check if there is cluster(s) in the segmentation results
#'
#' @description
#' This function determines if any clusters of change-points exist in the segmentation
#' results for a series. A cluster is defined as a group of change-points
#' separated at least by 80 days
#'
#' @param Series_df a data frame containing data for a main-nearby pair.
#' It contains at least two columns: the first column,labeled \code{date}, contains the dates in Date format (\%Y-\%m-\%d)
#' and the second contains one of the six difference series, including GPS - ERA, GPS - GPS',
#' GPS - ERA', ERA - ERA', GPS' - ERA', and GPS' - ERA, with names \code{GE, GGp,
#' GEp, EEp, GpEp,} and \code{GpE}, respectively
#' @param Name_series the name of the column in \code{Series_df} to check
#' @param CP a vector of dates in Date format (\%Y-\%m-\%d) given the change-points in the series to check
#' @param threshold a number indicating the maximum length of a cluster
#'
#' @return a vector of change-points which belongs to a cluster in Date format
#'
#' @export
#'
check_cluster_CP <- function(Series_df, Name_series, CP, threshold = 79) {

  if (length(CP) > 0 ) {
    begin_end = get_min_max_date(Series_df, Name_series)
    List_breaks <- sort(c(begin_end, CP))

    get_cluster <- function(X, thres = threshold){
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
    length_seg[1] <- sum(Series_df$date >= List_breaks[1] &
                           Series_df$date <= List_breaks[2])

    for (i in 2:(length(List_breaks) - 1)) {
      length_seg[i] <- sum(Series_df$date > List_breaks[i] &
                             Series_df$date <= List_breaks[i + 1])
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

    Res <- data.frame(changepoints = CP,
                      cluster = clust[c(-1, -length(clust))])
  } else {
    Res <- data.frame(cluster = 0)
  }

  return(Res)

}






