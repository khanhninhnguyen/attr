#' find the begin and end of data frame without NA
#'
#' @param data A data frame containing at least 2 columns named "Date" and
#' "Name_series". The "Date" column must be in Date format with the format "%Y-%m-%d",
#' and the "Name_series" column must be numeric (can contain NA values).
#'
#' @param column_name The name of the column in "data" that we want to check,
#' as a character string.
#'
#' @return Vector of Date which are the beginning and the end of time series
#' with available data in column_name
#'
#' @importFrom stats na.omit
#'
#' @keywords internal

get_min_max_date <- function(data, column_name) {

  Date <-  NULL

  df <- data[,c("Date",  column_name)]
  data[["Date"]] <- as.Date(data[["Date"]])
  df = na.omit(df)
  min_date = df$Date[1]
  max_date = df$Date[nrow(df)]
  y = c(min_date, max_date)
  return(y)
}
