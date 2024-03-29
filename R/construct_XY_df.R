#' Create a data frame containing X and Y (named "signal") for regression
#'
#' @param data A data frame containing at least 2 columns named "Date" and
#' "Name_series". The "Date" column must be in Date format with the format "%Y-%m-%d",
#' and the "Name_series" column must be numeric (can contain NA values).
#'
#' @param name_series The name of the column in "data" that we want to test,
#' as a character string.
#'
#' @param break_ind The index of the break if the offset is included in the model,
#' as an integer.
#'
#' @param Fourier  A logical value indicating whether to include Fourier elements
#'  in the model (TRUE or FALSE).
#'
#' @return A data frame used for testing that includes a column of the signal (Y)
#' and columns of X: 8 Fourier components (if selected), a global column, and a
#' jump column (if break_ind is not NULL).
#'
#' @importFrom dplyr select rename filter mutate %>%
#'
#' @importFrom tidyr complete
#'
#' @keywords internal

construct_XY_df <- function(data, name_series, break_ind = NULL, Fourier = TRUE){

  signal <- Date <- complete.time <- NULL

  data$Date = as.Date(data$Date, format = "%Y-%m-%d")

  Data_XY <- data %>%
    dplyr:: select(any_of(c(name_series,"Date"))) %>%
    dplyr:: rename(signal = name_series) %>%
    dplyr:: filter(!is.na(signal)) %>%
    tidyr:: complete(Date = seq(min(Date), max(Date), by = "day"))

  one.year = 365

  Data_XY <- Data_XY %>%
    mutate(complete.time = 1:nrow(Data_XY)) %>%
    select(-Date)

  if(isTRUE(Fourier)){
    for (i in 1:4){
      eval(parse(text=paste0("Data_XY <- Data_XY %>% mutate(cos",i,
                             "=cos(i*complete.time*(2*pi)/one.year),sin",i,
                             "=sin(i*complete.time*(2*pi)/one.year))")))
    }
  }

  Data_XY <- Data_XY %>%
    select(-complete.time)

  n0 = nrow(Data_XY)

  if(!is.null(break_ind)){
    Data_XY$right = c(rep(0, break_ind), rep(1, (n0-break_ind)))
  }

  Data_XY$left = rep(1, nrow(Data_XY))

  return(Data_XY)
}
