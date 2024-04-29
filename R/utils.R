#-----------------------------------------------------------------------------
# Utility Functions for Date Ranges in Data Frames
#-----------------------------------------------------------------------------

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
#'
#' @noRd


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

#-----------------------------------------------------------------------------
# Generalized Least Squares (GLS) Estimation Functions
#-----------------------------------------------------------------------------

#' Fit GLS (known variance-covariance matrix)
#'
#' @param phi  A numeric value indicating the parameter of the AR component.
#'
#' @param theta A numeric value indicating the parameter of the MA component.
#'
#' @param var.t A variance matrix computed by the moving window.
#'
#' @param df_XY A data frame resulted from construct_XY_df that includes a column
#' of the signal (Y) and columns of X: 8 Fourier components (if selected), a global
#' mean column, and a jump column (if break_ind is not NULL).
#'
#' @return A list containing the estimates and their standard errors, t-values,
#'  and p-values.
#'
#' @importFrom stats ARMAacf na.omit pnorm toeplitz
#'
#' @keywords internal
#'
#' @noRd
#'
Fit_GLS <- function(phi, theta, var.t, df_XY){

  signal <- NULL

  # cor_matrix fc to produce the variance corvariance matrix
  cor_matrix <- function(phi, theta, n.full, n.true){

    if(phi ==0 & theta ==0){
      y = rep(1, n.full)
    }else{
      y = ARMAacf(ar = phi, ma = theta, lag.max = n.true-1, pacf = FALSE)
    }

    m = toeplitz(y)

    return(m)
  }

  # compute the variance covariance matrix
  ind1 = which(is.na(df_XY$signal)==FALSE)
  var.t.na = var.t[ind1]
  var.matrix = diag(sqrt(var.t.na))
  cor.matrix = cor_matrix(phi, theta, n.full = nrow(df_XY), n.true = nrow(na.omit(df_XY)))

  if(phi==0 & theta==0){
    cov.var= diag(var.t.na)
  }else{
    cov.var0 = var.matrix  %*%  cor.matrix
    cov.var = cov.var0 %*%  var.matrix
  }

  # estimate
  X = as.matrix(df_XY%>% dplyr::select(-signal))
  X = X[ind1,]
  term2 = solve(cov.var)
  term1 = t(X) %*% term2 %*% X
  term3 = solve(term1)
  beta = term3 %*% t(X) %*% term2  %*% (as.matrix(df_XY$signal[ind1]))
  var.beta = term3

  # form the frame of result as in the ols
  residual = df_XY$signal[ind1] - X%*%beta
  t.val = beta/sqrt((diag(var.beta)))
  p.val = round(pnorm(-abs(t.val), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 4)
  fitted.gls = X%*%beta

  fit.gls = data.frame(Estimate = format(beta, scientific = TRUE, digits = 3),
                       `Std. Error` = format(sqrt(diag(var.beta)), scientific = TRUE, digits = 3),
                       `t value` = format(t.val, scientific = TRUE, digits = 3),
                       `Pr(>|t|)` = format(p.val, scientific = TRUE, digits = 3),
                       check.names = FALSE
  )

  Res <- list(Coefficients = beta,
              t.table = fit.gls,
              vcov = var.beta,
              residual = residual,
              fitted = fitted.gls)

  return(Res)

}

#-----------------------------------------------------------------------------
# Statistical Computation Functions
#-----------------------------------------------------------------------------

#' Compute the sliding window standard deviation by ScaleTau estimator
#'
#' @param Y A data frame contains at least 2 columns named Date and
#' name.var. The Date column must be in Date type on this format: "%Y-%m-%d"),
#' and the Name_series column is numeric (can be NA).
#'
#' @param name.var The name of the column in "Y" that we want to test,
#' as a character string.
#'
#' @param length.wind An integer value specifying the length of the window used
#' to compute the sliding standard deviation. Default is 60.
#'
#' @return A vector of values of the sliding standard deviation computed by the
#' scaleTau estimator.
#'
#' @importFrom robustbase scaleTau2
#'
#' @importFrom tidyr complete
#'
#' @importFrom stats na.omit approx
#'
#' @keywords internal
#'
#' @noRd

Sliding_std <- function(Y, name.var, length.wind = 60){

  # require date in the dataset, return std at time t
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = nrow(Y1)
  sigma.est1 = rep(NA, n)

  for (i in c(1:n)) {
    begin = max(i-(length.wind-1),1)
    end = min(n, i+length.wind)
    x.i = x[begin:end]
    x.idiff = (x.i)
    thre = 30
    if(i < 30 | i > (n-30)){
      thre = 16
    }
    if(length(na.omit(x.idiff)) <= thre){
      sd <- NA
    }else{
      sd <- robustbase::scaleTau2(na.omit(x.idiff), consistency = "finiteSample")
    }
    sigma.est1[i] <- sd
  }
  # linear regression of the variance for gaps  MAYBE REPLACE BY INTERPOLATION FUNCTION
  s = sigma.est1
  if (sum(is.na(s)) != 0 & sum(is.na(s)) != length(s)){
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    if(na.ind[1] == 1){
      ind.stop = which(is.na(s)==FALSE)[1]-1
      na.ind <- na.ind[-c(1:ind.stop)]
    }else if (is.na(s[n]) == 1){
      m = which(is.na(s) == FALSE)
      ind.start = m[length(m)]
      na.ind <- na.ind[-which(na.ind %in% c(ind.start:n))]
    }
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  sigma.est = s[which(Y1$date %in% Y$date)]

  return(sigma.est)

}

#-----------------------------------------------------------------------------
# Data Preparation Functions
#-----------------------------------------------------------------------------

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
#' @noRd

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
