#' Test Jump Significance in Difference Series
#'
#' @description
#' This function makes the test of the jump at change-points in a series using Feasible
#' Generalized Least Squares (FGLS)
#'
#' @param Series_df a data frame containing data for a main-nearby pair.
#' It contains at least two columns: the first column,labeled \code{date}, contains the dates in Date format (%Y-%m-%d)
#' and the second contains one of the six difference series, including GPS - ERA, GPS - GPS',
#' GPS - ERA', ERA - ERA', GPS' - ERA', and GPS' - ERA, with names \code{GE, GGp,
#' GEp, EEp, GpEp,} and \code{GpE}, respectively
#' @param Name_series the name of the column in \code{Series_df} to check
#' @param CP a vector of dates in Date format (%Y-%m-%d) representing the change-points in the series to check
#' @param limit a integer specifying the number of points in the segments
#' before and after the change-point which will be used for testing. If NULL (default), the segments
#' before and after each tested change-point are defined by the nearest change-points
#' @param tol a numeric value indicating the stopping tolerance on the parameter in the iterative FGLS algorithm. Default is 0.01
#' @param noise_model a string indicating one noise model. It must be one of the following: \code{'AR(1)', 'MA(1)', 'ARMA(1,1)',
#' or 'White'}
#' @param length_win an integer value specifying the length of the window used to compute the moving variance. Default is 60
#' @param name_case the name of the file to save the detailed results. By default, the detailed results are not saved and default is NULL
#'
#' @return A list of 2 elements:
#' \itemize{
#' \item \code{Used_dates} a vector of dates used for the test,
#' \item \code{Summary_tab} a summary data frame of the test results size l x 4, where l
#' is the number of variables in the model, and the 4 columns correspond to the
#' estimates, their standard errors, t-values, and p-values respectively. By default, 10 variables
#' are included in the output, such as the global mean, periodic bias (modeled by 4 coupled
#' sine and cosine Fourier series), and the offset.
#' }
#'
#' @importFrom dplyr filter slice_tail slice_head arrange all_of
#' @importFrom tidyr complete
#' @importFrom stringr str_c
#'
#' @export
#'
Test_CP <- function(Series_df, Name_series, CP, limit = NULL,
                             tol = 0.01, noise_model, length_win = 60, name_case = NULL){

  date <- .data <- wts <- NULL

  check_model = length(noise_model %in% c("AR(1)", "MA(1)", "ARMA(1,1)", "White"))
  stopifnot("noise_model must be one of the following: 'AR(1)', 'MA(1)', 'ARMA(1,1)', or 'White'}" = (check_model == 1))
  stopifnot("Series_df must be a dataframe" = is.data.frame(Series_df))

  noise_model = convert_model_order(noise_model)
  #####
  # count time
  start_time <- Sys.time()

  # creat the XY dataframe
  Series_df <- Series_df %>%
    dplyr:: select(all_of(c(Name_series,"date"))) %>%
    dplyr:: filter(!is.na(.data[[Name_series]])) %>%
    tidyr:: complete(date = seq(min(date), max(date), by = "day"))

  # check if the time series cover the Breakpoint
  cover <- (Series_df$date[1] < CP) & (Series_df$date[length(Series_df$date)] > CP)
  stopifnot("Changepoint must be in the Series_df" = isTRUE(cover))

  if (!is.null(limit)){

    before_data <- Series_df %>%
      dplyr:: filter(date <= CP, !is.na(.data[[Name_series]])) %>%
      dplyr:: arrange(date) %>%
      dplyr:: slice_tail(n = limit)

    # Extract data 1000 points after the given date, including the closest points not NA
    after_data <- Series_df %>%
      dplyr:: filter(date > CP, !is.na({{Name_series}})) %>%
      dplyr:: arrange(date) %>%
      dplyr:: slice_head(n = limit)

    Series_df <- rbind(before_data, after_data) %>%
      tidyr:: complete(date = seq(min(before_data$date),
                                  max(after_data$date), by = "day"))
  }

  df_XY = construct_XY_df(Series_df,
                          name_series = Name_series,
                          CP_ind = which(Series_df$date == CP))

  coef_arma0 = list(phi = 0, theta = 0)
  n_full = nrow(df_XY)
  residual0 = rep(NA, n_full)
  ind_avai = which(!is.na(df_XY$signal))

  # call expression
  name_var <- colnames(df_XY)[colnames(df_XY) != "signal"]
  mod_X <- name_var %>% stringr::str_c(collapse = "+")
  mod_expression <- c("signal","~",mod_X, -1) %>% stringr::str_c(collapse = "")

  # Fit the Ordinary Least Square
  fit_OLS = lm( mod_expression, data = df_XY)
  residual0[ind_avai] <- fit_OLS$residuals
  coef_OLS = fit_OLS$coefficients

  # estimate initial moving variance

  Y0 = data.frame(date = Series_df$date, residus = residual0)
  std0 <- Sliding_std(Y = Y0, name.var = "residus")
  d1 = 10
  epsilon = 10
  i = 0

  while (epsilon > tol) {

    fit_WLS <- lm(mod_expression, data = cbind(df_XY, wts = 1/((std0)^2)), weights = wts)

    residual1 <- rep(NA, n_full)
    residual1[ind_avai] <- fit_WLS$residuals
    norm_residual = residual1/std0

    # Fit ARMA model
    coef_arma = list(phi = 0, theta =0)

    if (!is.null(noise_model)){
      arima.fit = arima(x = norm_residual, order = noise_model, include.mean = FALSE)
      coef_arma$phi = ifelse("ar1" %in% names(arima.fit$coef), arima.fit$coef["ar1"], 0)
      coef_arma$theta = ifelse("ma1" %in% names(arima.fit$coef), arima.fit$coef["ma1"], 0)
    }

    fit.gls = Fit_GLS(phi = coef_arma$phi,
                  theta = coef_arma$theta,
                  var.t = (std0)^2,
                  df_XY = df_XY)
    # condition to stop
    coef_diff_norm = (fit.gls$Coefficients - coef_OLS)/(sqrt(diag(fit.gls$vcov)))
    epsilon = d1 - sum((coef_diff_norm)^2)
    d1 = sum((coef_diff_norm)^2)
    coef_OLS = fit.gls$Coefficients

    # Stop with more than 10 iteration and
    # report if the change in coefficients is very small

    i = 1+i
    if(i>10){
      print("Check again noise model")
      print("Number of iteration: ",i)
      break
    }

    if (!is.null(noise_model)){
      if(identical(noise_model,c(0,0,1))){
        d3 = abs(coef_arma$theta - coef_arma0$theta)
      }else{
        d3 = abs(coef_arma$phi - coef_arma0$phi)
      }
      epsilon = d3
    }
    # update the residual and moving variance
    residual = rep(NA, n_full)
    residual[ind_avai] = fit.gls$residual
    Y0 = data.frame(date = Series_df$date, residus = residual)
    std0 = Sliding_std(Y = Y0, name.var = "residus")
    coef_arma0 = coef_arma

    # print(paste0("Iteration:  ", i))
  }

  fit_val = rep(NA, n_full)
  fit_val[ind_avai] = fit.gls$fitted

  df_XY = df_XY %>%
    mutate(Day_list = Series_df$date)

  end_time <- Sys.time()

  Res_FGLS = (list(coefficients = fit.gls$Coefficients,
                   mw_std = std0,
                   residual = residual0,
                   fit = fit_val,
                   t.table = fit.gls$t.table,
                   coef_arma = coef_arma,
                   no_iteration = i,
                   coef_diff_norm = coef_diff_norm,
                   epsilon = epsilon,
                   all_last_GLS = fit.gls,
                   t = (end_time - start_time),
                   df_XY = df_XY))

  if (!is.null(name_case)) {
    save(Res_FGLS,
         file = paste0("Res_FGLS_",
                       name_case,
                       Name_series,
                       CP,
                       ".RData"))
  }

  summary_tab <- as.data.frame(apply(fit.gls$t.table, 2, function(x) as.numeric(x)))
  rownames(summary_tab) <- rownames(fit.gls$t.table)

  Res <- list(Summary_tab = summary_tab,
              Used_dates = Series_df$date)

  return(Res)

}

