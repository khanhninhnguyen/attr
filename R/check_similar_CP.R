#' To check for similar changes in nearby stations
#'
#' @description
#'
#' This function identifies "similar" change-points in a nearby station that occur
#' within a specified number of days of the change-points in the main station.
#' Data between change-points in the main station and "similar" change-points
#' in the nearby station is removed from the four non-collocated series for the
#' significance test of the change in mean.
#'
#' @param Series_df A data frame containing data for a main-nearby pair. It contains seven
#' columns: the first labeled \code{Date} with dates in "YYYY-MM-DD" format and the
#' remaining six columns containing numeric values representing six different
#' series: GPS - ERA, GPS - GPS', GPS - ERA', ERA - ERA', GPS' - ERA', and GPS' - ERA,
#' with names \code{GE, GGp, GEp, EEp, GpEp,} and \code{GpE}, respectively.
#'
#' @param main_cp A vector of dates in Date format: "%Y-%m-%d" representing
#' change-points in the main station.
#'
#' @param nearby_cp_one A vector of dates in Date format: "%Y-%m-%d" representing
#' change-points in the nearby station.
#'
#' @param threshold An integer specifying the number of points between two change-points,
#' which defines if it is a "similar" change-point. Default is 10.
#'
#' @return A vector of dates in Date format: "%Y-%m-%d" representing "similar"
#' change-points in the nearby station to the change-points in the main station.
#'
#' @importFrom tidyr drop_na
#'
#' @importFrom dplyr arrange
#'
#' @export
#'

check_similar_CP <- function(Series_df, main_cp, nearby_cp_one, threshold){

  Date <- GE <- NULL

  sim_changes <- as.Date(character(0))

  main_df <- Series_df %>%
    select(Date, GE) %>%
    drop_na()

  if(length( nearby_cp_one) != 0){
    for (i in c(1:length(main_cp))) {
      main_brp = main_cp[i]

      after <- main_df$Date[which(main_df$Date > main_brp)][threshold]
      before <- tail(main_df$Date[which(main_df$Date < main_brp)],threshold)[1]

      ind_sim = which( nearby_cp_one >= before &  nearby_cp_one <= after)

      if(length(ind_sim) > 0){
        sim_changes <- c(sim_changes,  nearby_cp_one[ind_sim])

      }
    }
  }

  return(sim_changes)
}
