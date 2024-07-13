#' Check for "similar" change-points in nearby stations
#'
#' @description
#' This function identifies "similar" change-points in a nearby station that occur
#' within a specified number of days of the change-points in the main station.
#' Data between change-points in the main station and "similar" change-points
#' in the nearby station is removed from the four non-collocated series for the
#' significance test of the change in mean
#'
#' @param Series_df a data frame containing data for a main-nearby pair. It contains seven
#' columns: the first column,labeled \code{date}, contains the dates in Date format (%Y-%m-%d) and the
#' remaining six columns contain the six numeric signals for the six different
#' series, including GPS - ERA, GPS - GPS', GPS - ERA', ERA - ERA', GPS' - ERA', and GPS' - ERA,
#' with names \code{GE, GGp, GEp, EEp, GpEp,} and \code{GpE}, respectively
#' @param main_cp a vector of dates in Date format (%Y-%m-%d) representing the change-points in the main station
#' @param nearby_cp_one a vector of dates in Date format (%Y-%m-%d) representing the change-points in the nearby station
#' @param threshold an integer specifying the number of points between two change-points,
#' which defines if it is a "similar" change-point. Default is 10
#'
#' @return a vector of dates in Date format (%Y-%m-%d) representing the "similar"
#' change-points in the nearby station to the change-points in the main station
#'
#' @importFrom tidyr drop_na
#' @importFrom dplyr select arrange %>%
#'
#' @export
#'

check_similar_CP <- function(Series_df, main_cp, nearby_cp_one, threshold){

  date <- GE <- NULL

  sim_changes <- as.Date(character(0))

  main_df <- Series_df %>%
    dplyr:: select(date, GE) %>%
    tidyr:: drop_na()

  if(length( nearby_cp_one) != 0) {
    for (i in c(1:length(main_cp))) {
      main_brp = main_cp[i]

      after <- main_df$date[which(main_df$date > main_brp)][threshold]
      before <- tail(main_df$date[which(main_df$date < main_brp)],threshold)[1]

      ind_sim = which( nearby_cp_one >= before &  nearby_cp_one <= after)

      if(length(ind_sim) > 0){
        sim_changes <- c(sim_changes,  nearby_cp_one[ind_sim])

      }
    }
  }

  if(length(sim_changes) == 0) {
    sim_changes = "No similar change-points"
  }

  return(sim_changes)
}
