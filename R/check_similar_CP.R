#' To check is there similar changes in the nearby stations
#'
#' @param Six_Series A data frame contains 6 time series differences (so 7 columns
#'in total) between a main station and a nearby station. The first column should
#' be the date (in Date type on this format: "\%Y-\%m-\%d"), and the next six columns
#' must be in this order:
#' G-E, G-G', G-E', E-E', G'-E', G'-E.
#'
#' @param main_cp A vector of main station breakpoints in Date type on this format:
#' "\%Y-\%m-\%d".
#'
#' @param nearby_cp_one A vector of main station breakpoints in Date type on this format:
#' "\%Y-\%m-\%d".
#'
#' @param threshold Integer number is the threshold to define the number of
#' point between 2 breakpoints in a cluster
#'
#' @return list of main changes and their similar in the nearby in Date type
#'
#' @importFrom tidyr drop_na
#' @importFrom dplyr arrange

#' @export
#' @keywords internal

check_similar_CP <- function(Six_Series, main_cp, nearby_cp_one, threshold){

  Date <- GE <- NULL

  sim_changes <- as.Date(character(0))

  main_df <- Six_Series %>%
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
