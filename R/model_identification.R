#' Identify Time Series Noise Model
#'
#' @param dataset A list of data frames, where each data frame contains 6 time series
#' differences (so 7 columns in total) between a main station and a nearby station. The
#' first column should be the date, and the next six columns must be in this order:
#' G-E, G-G', G-E', E-E', G'-E', G'-E. The name of each data frame must be the name
#' of the nearby station.
#'
#' @param main_break A vector of main station breakpoints in Date type on this format:
#' "\%Y-\%m-\%d".
#'
#' @param nearby_break A list of break point vectors for each nearby station, named
#' after the respective nearby station, in Date type on this format: "\%Y-\%m-\%d".
#'
#' @return  A noise model data frame with dimensions (6 * n), where n is the
#' number of nearby stations. The 6 rows correspond to the 6 difference series.
#' The noise model for each series is given as a character string, such as "AR(1)",
#' "MA(1)", "ARMA(1,1)", or "White"
#'
#' @importFrom stringr str_c
#'
#' @importFrom lmtest coeftest
#'
#' @importFrom forecast auto.arima arimaorder
#'
#' @importFrom tidyr complete
#'
#' @importFrom nlme gls varFixed
#'
#' @importFrom stats arima lm setNames

#' @export
#'
model_identification <- function(dataset, main_break, nearby_break){

  #####
  # pre-check
  stopifnot("dataset must be a list of dataframe" = is.list(dataset))

  stopifnot("dataset must be a list of dataframe" = all(sapply(dataset, is.data.frame)))

  stopifnot("nearby_break must be a list" = is.list(nearby_break))

  stopifnot("names of dataset must be the same to nearby_break" = setequal(names(nearby_break), names(dataset)))

  stopifnot("dataset and nearby_break must have the same length" = (length(dataset) == length(nearby_break)))


  if (!inherits( main_break, "Date") | any(!sapply(nearby_break, inherits, "Date"))) {
    stop("Must specify both main_break and nearby_break in Date type")
  }

  if( any(!sapply(dataset, function(df) {
    ncol(df) == 7 & inherits(df[,1], "Date")
  }))) {
    stop("Each data frame in the dataset need to have 7 columns &\
         the first column must be in Date type")
  }

  names_col = c("Date", "GE", "GGp", "GEp", "EEp", "GpEp", "GpE")

  dataset <- lapply(dataset, function(df) {
    colnames(df) <-  names_col
    return(df)
  })

  #####
  # Identify_model function
  identify_model <- function(Series_df, Name_series, Breakpoints,
                             begin_day = NULL, end_day = NULL, tol0 = 0.01, name_case = NULL, ...) {

    Date <- Freq <- NULL

    stopifnot("Six_Series must be a dataframe" = is.data.frame(Series_df))

    estim <- data.frame(
      order = rep(NA, 2),
      coef = rep(NA, 2),
      p_value = rep(NA, 2)
    )
    # check also begin and end is Date type

    if (is.null(begin_day) & is.null(end_day)){

      begin_end = get_min_max_date(Series_df, Name_series)
      List_breaks <- sort(c(begin_end, Breakpoints))
      longest_seg_ind = which.max(diff(List_breaks, lag = 1))
      list_day_avai = Series_df$Date[which(!is.na(Series_df[[Name_series]]))]
      begin_day_0 = List_breaks[longest_seg_ind]
      begin_day = list_day_avai[which(list_day_avai > begin_day_0)[1]]
      end_day = List_breaks[longest_seg_ind + 1]

    } else{

      if (!inherits(begin_day, "Date") | !inherits(end_day, "Date")) {
        stop("Must specify both begin_day and end_day in Date type")
      }

    }

    df_data <- Series_df %>%
      dplyr::filter(Date >= begin_day & Date <= end_day) %>%
      tidyr::complete(Date = seq(begin_day, end_day, by = "day"))

    # Run the regression

    m = construct_XY_df(df_data, name_series = Name_series)
    nna_ind = which(!is.na(m$signal))
    norm_res <- rep(NA, nrow(m))

    fit_igls = IWLS(df_XY = m[nna_ind,], tol = tol0, day.list = df_data$Date[nna_ind])
    norm_res[nna_ind] = unlist(fit_igls$residual)/sqrt(unlist(fit_igls$var))
    arima_fit = fit_arima(norm_res)
    estim$order = arima_fit$pq
    estim$p_value = round(arima_fit$p, digits = 4)
    estim$coef = round(arima_fit$coef, digits = 4)
    rownames(estim) <- c("ar", "ma")

    Res <- list(fit = estim,
                residual = unlist(fit_igls$residual),
                variance = unlist(fit_igls$var),
                coef = round(fit_igls$coefficients, digits = 4))

    if(!is.null(name_case)){
      save(Res, file = paste0("Noise_model_", name_case, ".RData"))
    }

    return(estim)

  }

  # IWLS: iterative weighted least square fc
  IWLS <- function(df_XY, tol, day.list){
    resi0 = rep(NA, nrow(df_XY))

    # call expression
    list.para <- colnames(df_XY)[2:dim(df_XY)[2]]
    mod.X <-  list.para %>% stringr::str_c(collapse = "+")
    mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")

    # ols
    ols.fit = lm( mod.expression, data = df_XY)
    resi0[which(is.na(df_XY$signal)==FALSE)] <- ols.fit$residuals
    old.coef = ols.fit$coefficients

    # estimate initial moving variance
    Y0 = data.frame(date = day.list, residus = resi0)
    w0 = Sliding_std(Y = Y0, name.var = "residus", length.wind = 60)
    change1 = 10
    i=0

    while (change1 > tol) {
      df_XY$w = w0^2
      gls.fit = eval(parse(text=paste0("nlme::gls(",mod.expression,",data=df_XY, correlation = NULL, na.action = na.omit, weights=varFixed(value = ~w)",")")))
      change1 = sum((gls.fit$coefficients - old.coef)^2)
      # deg =  cbind(rep(1, nrow(df_XY)), as.matrix(df_XY[,c(2:9)]))
      fit.val = as.matrix(df_XY[,list.para]) %*% as.matrix(gls.fit$coefficients)
      resi0 = df_XY$signal - fit.val
      Y0 = data.frame(date = day.list, residus = resi0)
      w0 = Sliding_std(Y = Y0, name.var = "residus", length.wind = 60)
      old.coef = gls.fit$coefficients
      i=1+i

      if (i>10) {
        break
      }
    }

    Res <- list( coefficients = gls.fit$coefficients,
                 var = w0^2,
                 residual = as.numeric(resi0))
    return(Res)
  }

  fit_arima <- function(signal.test, significant.level = 0.05){

    Date <- .data <- signal <- NULL

    pre_fit <- forecast::auto.arima( signal.test,
                                     d = 0,
                                     ic = "bic",
                                     seasonal = FALSE,
                                     stationary = TRUE,
                                     allowmean =FALSE,
                                     lambda = NULL,
                                     max.p = 1,
                                     max.q = 1,
                                     start.p = 0,
                                     trace = FALSE,
                                     allowdrift = FALSE,
                                     approximation=FALSE
    )

    pq = forecast::arimaorder(pre_fit)
    coeffient = rep(NA, 3)
    coeffient[which(pq!= 0)] = pre_fit$coef

    sig_test <- lmtest::coeftest(pre_fit)
    p_vals = sig_test[,4]
    p_vals =  p_vals[which(rownames(sig_test) %in% c("ar1", "ma1"))]
    sig_order = rownames(sig_test)[which(p_vals < significant.level)]

    if(any(sig_order %in% c("ar1"))){
      p = 1
    }else{
      p = 0
    }

    if(any(sig_order %in% c("ma1"))){
      q = 1
    }else{
      q = 0
    }

    if(sum(pq) != (p+q)){
      coeffient = rep(NA, 3)
      pq = c( p, 0, q)
      fitARIMA = try(arima( signal, pq, method="ML"), TRUE)
      coeffient[which(pq!= 0)] = fitARIMA$coef
      sig_test <- lmtest::coeftest(fitARIMA)
      p_vals = sig_test[,4]
      p_vals =  p_vals[which(rownames(sig_test) %in% c("ar1", "ma1"))]
    }

    p_values = rep(NA, 3)
    p_values[which(pq!=0)] <- p_vals
    pq = pq[-2]
    coef = coeffient[-2]
    p_values = p_values[-2]

    return(list(pq = pq, coef = coef, p = p_values))
  }


  #####
  transform_model <- function(model_order){
    model_order = as.numeric(model_order)
    if (identical(model_order, c(1,0))) {
      y = "AR(1)"
    } else if (identical(model_order, c(0,1))) {
      y = "MA(1)"
    } else if(identical(model_order, c(1,1))) {
      y = "ARMA(1,1)"
    } else {
      y = "White"
    }
    return(y)
  }

  # for the main series
  GE_mod = identify_model(dataset[[1]], Name_series = "GE", Breakpoints = main_break)
  main_model = transform_model(GE_mod$order)
  # for the other series
  all_5_model <- sapply(names(dataset), function(x){
    List_joint = sort(c(main_break, nearby_break[[x]]))
    Break_points_six_series <- list(GGp = List_joint,
                                    GEp = List_joint,
                                    EEp = List_joint,
                                    GpEp = nearby_break[[x]],
                                    GpE = List_joint
    )
    sapply(names_col[3:7], function(y){
      all_mod = identify_model(dataset[[x]],
                               Name_series = y,
                               Breakpoints = Break_points_six_series[[y]])
      transform_model(all_mod$order)
    })
  })

  noise_model <- setNames(data.frame(t(rep(main_model, length(dataset)))),
                          names(dataset))
  noise_model = rbind( noise_model,
                       all_5_model)
  rownames(noise_model) <- names_col[-1]

  return(noise_model)

}
