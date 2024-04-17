#' To check is there similar changes in the nearby stations
#'
#' @param Six_Series A data frame contains 6 time series differences (so 7 columns
#'in total) between a main station and a nearby station. The first column should
#' be the date (in Date type on this format: "\%Y-\%m-\%d"), and the next six columns
#' must be in this order:
#' G-E, G-G', G-E', E-E', G'-E', G'-E.
#'
#' @param List_break_main A vector of main station breakpoints in Date type on this format:
#' "\%Y-\%m-\%d".
#'
#' @param List_break_nearby A vector of main station breakpoints in Date type on this format:
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

check_similar_CP <- function(Six_Series, List_break_main, List_break_nearby, threshold){

  Date <- GE <- NULL

  sim_changes <- as.Date(character(0))

  main_df <- Six_Series %>%
    select(Date, GE) %>%
    drop_na()

  if(length(List_break_nearby) != 0){
    for (i in c(1:length(List_break_main))) {
      main_brp = List_break_main[i]

      after <- main_df$Date[which(main_df$Date > main_brp)][threshold]
      before <- tail(main_df$Date[which(main_df$Date < main_brp)],threshold)[1]

      ind_sim = which(List_break_nearby >= before & List_break_nearby <= after)

      if(length(ind_sim) > 0){
        sim_changes <- c(sim_changes, List_break_nearby[ind_sim])

      }
    }
  }

  return(sim_changes)
}
