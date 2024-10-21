#' Predict the Origin of Change-Points
#'
#' @description
#'
#' This function attributes change-points detected in the difference series to one
#' or both contributing series. Assuming the primary series is GPS and the secondary
#' is ERA, change-points are identified in the GPS - ERA series. The
#' function determines whether these change-points are due to GPS, ERA, or both.
#' The method utilizes data from a nearby station providing GPS' and ERA' series.
#' Users must compute six difference series (GPS-ERA, GPS-GPS', GPS-ERA', ERA-ERA',
#'  GPS'-ERA', GPS'-ERA) for each main-nearby station pair to serve as input.
#' This function tests the significance of the jump at a given change-point in each of the six
#' difference series for each main-nearby station pair, predicts configurations
#' from the test results, and aggregates them from all nearby stations.
#' Verify the format of the input in the "tests/testthat/test-attribution.R" file.
#'
#' @details
#' The function operates in four main steps:
#'
#' \strong{1. Model Identification} (optional): identify the noise model for each difference series
#' using the \code{\link{NoiseModel_Id}} function.
#'
#' \strong{2. Significance Test}: test the significance of the jump in each difference series
#' using the \code{\link{Test_CP}} function.
#'
#' \strong{3. Configuration Prediction}: predict the configuration of each triplet (main, change-point, nearby) based
#' on their test result (t-value) of the six difference series using the \code{\link{Prediction_CP}} function.
#'
#' \strong{4. Aggregation}: combine the prediction results for each change-point from all nearby stations.
#'
#' Additionally, the function performs two preliminary checks:
#'
#' \strong{a. Cluster Check}: determine if any clusters of change-points exist across all stations
#' using the \code{\link{check_cluster_CP}} function. A
#' cluster is defined as a group of change-points within 80 days. If a cluster
#' is found, the \code{\link{Attribution_CP}} function stops.
#'
#' \strong{b. Coincident Check}: identify "similar" change-points in nearby stations that occur within 10 days
#' of change-points in the main station using the \code{\link{check_similar_CP}} function.
#' Data between change-points in the main station and "similar" change-points in the nearby station is removed from the
#' four non-collocated series for the significance test of the change in mean.
#'
#' @param dataset a list of length n where n is the number of nearby stations.
#' Each element of the list is a data frame representing a main-nearby station
#' pair and is named after the nearby station. contains seven
#' columns: the first column,labeled \code{date}, contains the dates in Date format (\%Y-\%m-\%d) and the
#' remaining six columns contain the six numeric signals for the six different
#' series, including GPS - ERA, GPS - GPS', GPS - ERA', ERA - ERA', GPS' - ERA', and GPS' - ERA,
#' with names \code{GE, GGp, GEp, EEp, GpEp,} and \code{GpE}, respectively
#' @param main_cp a vector of dates in Date format (\%Y-\%m-\%d) representing the change-points in the main station
#' @param nearby_cp A list with n elements, where n is the number of nearby stations.
#' Each element is named after the corresponding nearby station and is a vector of dates in
#' Date format (\%Y-\%m-\%d) representing the change-points in the corresponding nearby station
#' @param noise_model_fix a parameter specifying the noise model for the series. It can be:
#'
#'   - a string indicating one noise model for all series.
#'
#'   - a data frame with 6 rows and n columns, where each row represents a difference
#'   series (e.g, \code{GE, GGp}, etc) and each column represents a nearby station,
#'   specifying the noise model for each series.
#'
#'   - NULL, to determine noise models using the \code{\link{NoiseModel_Id}} function.
#'
#' The noise model must be one of the following: \code{'AR(1)', 'MA(1)', 'ARMA(1,1)',
#' or 'White'}.
#' @param nearby_weight a numeric vector of length n, each element representing
#' the weight assigned to a nearby station for aggregating predicted configurations.
#' If NULL, equal weight is assigned to all stations.
#' @param limit_2side an integer specifying the number of points in the segments
#' before and after the change-point used for testing. If NULL (default), the segments
#' before and after each tested change-point are defined by the nearest change-points.
#' @param save_result any values different to NULL specifying test results (summary
#' table and date used for test) will be saved
#' @param name_main a string specifying name of the main station (optional)
#' @param lmin an integer specifying minimum number of point in both sides for test.
#' At default, lmin = 0.
#'
#'
#' @return A list with four components:
#' \describe{
#'   \item{Noise_model}{
#' a data frame with size 6xn where:
#'
#'  * the 6 rows represent the six difference series,
#'
#'  * the n columns represent the n nearby stations,
#'
#'  * each cell contains a string indicating the noise model (e.g., 'AR(1)',
#' MA(1)') of the according series. }
#'
#'   \item{Test_results}{
#'.    a list of m data frames, where m is the number of change-points. Each data
#'     frame contains t-values for the test jump at the change-point of six difference series,
#'     assessed for all main-nearby station pairs. The data frames are named after
#'     each change-point and structured with 6 rows for the six difference series
#'     and n columns corresponding to the nearby stations.
#'   }
#'
#'   \item{Prediction_results}{
#'    a list of m data frames, where m is the number of change-points. Each data
#'    frame is named after the corresponding change-point and contains the predicted
#'    configurations for each change-point. The data frames include three columns:
#'    \code{max_freq_config} is the candidate configuration with the highest
#'    frequency from running 24 models, \code{freq} is the associated frequency,
#'    and \code{final_config} is the final configuration (see the \code{\link{Prediction_CP}} function).
#'    Each row corresponds to a nearby station.
#'   }
#'
#'   \item{Aggregated_results}{
#'     a list of m data frames, where m is the number of change-points. Each data
#'     frame contains the aggregated prediction results for each change-point, based
#'     on the weights specified in \code{nearby_weight}. The data frames include three columns:
#'    \code{max_freq_config} is the candidate configuration with the highest
#'    weighted frequency from all nearby stations, \code{freq} is the associated frequency,
#'    and \code{final_config} is the final configuration , i.e. \code{max_freq_config}
#'    if the \code{max_freq_config} is unique. If \code{max_freq_config} is not unique,
#'    the \code{final_config_all} is the one with the highest prior probability of the table given in the article
#'   }
#' }
#'
#'
#' @importFrom dplyr %>% mutate select rename any_of
#' @importFrom utils tail
#'
#' @export
#'
#'
Attribution_CP <- function(dataset,
                           main_cp,
                           nearby_cp,
                           noise_model_fix = NULL,
                           nearby_weight = NULL,
                           name_main = NULL,
                           save_result = 0,
                           limit_2side = 100,
                           lmin = 0){
  options(warn = 2)
  #####
  # pre-check
  date <- candidate_config <- max_freq_config <-  NULL

  stopifnot("dataset must be a list of dataframe" = is.list(dataset))

  stopifnot("dataset must be a list of dataframe" = all(sapply(dataset, is.data.frame)))

  stopifnot("nearby_cp must be a list" = is.list(nearby_cp))

  stopifnot("names of dataset must be the same to nearby_cp" = setequal(names(nearby_cp), names(dataset)))

  stopifnot("dataset and nearby_cp must have the same length" = (length(dataset) == length(nearby_cp)))


  if (!inherits( main_cp, "Date") | any(!sapply(nearby_cp, inherits, "Date"))) {
    stop("Must specify both main_cp and nearby_cp in date format")
  }

  if( any(!sapply(dataset, function(df) {
    ncol(df) == 7 & inherits(df[,1], "Date")
  }))) {
    stop("Each data frame in the dataset need to have 7 columns &\
         the first column must be in date format")

  }

  # check dataset include 6 series for each nearby
  list_six_diff = c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE")

  required_columns = sapply(c(1:length(dataset)), function(i) all(list_six_diff %in% names(dataset[[i]])))
  stopifnot("Columns in test_result must include \
            : date, GE, GGp, GEp, EEp, GpEp, GpE" = all(required_columns))

  nearby_name <- names(dataset)
  #####
  # check cluster

  cluster_main = check_cluster_CP(Series_df = dataset[[1]],
                                  Name_series = "GE",
                                  CP = main_cp)

  cluster_other = lapply(nearby_name, function(x){
    cluster_nearby = check_cluster_CP(Series_df = dataset[[x]],
                                      Name_series = "GpEp",
                                      CP = nearby_cp[[x]])
  })

  check_clusters <- sapply(cluster_other, function(df) any(df$cluster != 0))

  warn_text_cluster = " Must not have any cluster of change points"

  if(any(cluster_main$cluster != 0)) {
    stop(paste0(warn_text_cluster, " in the main series"))
  }

  if(any(check_clusters)) {
    stop(paste0(warn_text_cluster,
                " in the nearby series: ",
                nearby_name[check_clusters]))
  }

  #####
  # Model identification
  if (is.null(noise_model_fix)) {
    Noise_model <- NoiseModel_Id(dataset, main_cp, nearby_cp)
  } else if (class(noise_model_fix[1]) == "character") {
    Noise_model <- data.frame(matrix(noise_model_fix,
                                     ncol = length(nearby_name),
                                     nrow = 6))
    rownames(Noise_model) <- list_six_diff
    colnames(Noise_model) <- nearby_name
  } else {
    Noise_model = noise_model_fix
    rownames(Noise_model) <- list_six_diff
    colnames(Noise_model) <- nearby_name
  }

  #####
  # remove the similar change points in the nearby for the test

  similar_cp <- lapply(nearby_name, function(x) {
    check_similar_CP(Series_df = dataset[[x]],
                     main_cp = main_cp,
                     nearby_cp_one = nearby_cp[[x]],
                     threshold = 10)
  })

  check_similar <- sapply(similar_cp, function(x) length(x)>0)

  nearby_cp_m = nearby_cp

  if(any(check_similar)) {
    ind_change = which(check_similar)
    for (i in ind_change) {
      list_brp_nb <- nearby_cp[[nearby_name[i]]]
      nearby_cp_m[[nearby_name[i]]] <- list_brp_nb [ which(!list_brp_nb %in% similar_cp [[i]])]
    }
  }

  #####
  test_sig <- function(List_CP, CP, df, Name_series, noise_model, limit_period, name_case, save_result) {

    begin <-  List_CP[tail(which(List_CP < CP),1)]
    end <- List_CP[which( List_CP > CP)[1]]

    if (length(begin) == 0 | length(end) == 0) {

      t_val_jump <- NA

    } else if (length(begin) > 0 & length(end) > 0) {

      data_test <- df %>%
        dplyr:: filter(date > begin & date <= end)

      if (is.null(save_result)){
        name_case <- NULL
      }

      # Regression
      fit_fgls <- Test_CP(Series_df = data_test,
                          Name_series = Name_series,
                          CP = CP,
                          noise_model = noise_model,
                          limit = limit_period,
                          lmin = lmin,
                          save_result = save_result,
                          name_case = name_case)

      t_val <- fit_fgls$Summary_tab$`t value`
      t_val_jump <- t_val[which(rownames(fit_fgls$Summary_tab) == "jump")]
    }

    return(t_val_jump)
  }
  # Significance test for the jump in the main

  List_CP_main <- sort( c( main_cp,
                           get_min_max_date(data = dataset[[1]],
                                            column_name = "GE")))

  # Significance test for the other 5 test
  t_val_all = lapply(nearby_name, function(x) {
    t <- NA

    nearby_cp_z <- sort(c(nearby_cp_m[[x]],
                          get_min_max_date(data = dataset[[x]],
                                           column_name = "GpEp")))

    List_joint = sort(c(List_CP_main, nearby_cp_z))

    CP_six_series <- list( GE = List_joint,
                           GGp = List_joint,
                           GEp = List_joint,
                           EEp = List_joint,
                           GpEp = List_joint,
                           GpE = List_joint
    )

    t_val_all_brk <- sapply(main_cp, function(y){

      five_t_val <- sapply(list_six_diff, function(z){

        noise_model_z = unlist(Noise_model[which(rownames(Noise_model) == z), x])

        tryCatch({
          # Attempt to run the problematic function
          t <- test_sig(List_CP = CP_six_series[[z]],
                        CP = y,
                        df = dataset[[x]],
                        Name_series = z,
                        noise_model = noise_model_z,
                        limit_period = limit_2side,
                        name_case = paste0(x,"-", name_main),
                        save_result = save_result)

        }, error = function(e) {
          # If there's an error, return NA
          t <- NA
        })
        setNames(t, z)
        return(t)
      })
      return(five_t_val)

    })

    colnames(t_val_all_brk) <- as.character(main_cp)

    return(t_val_all_brk)

  })
  names(t_val_all) <- nearby_name

  #####
  # Prediction
  t_val_all <- lapply(main_cp, function(x) {
    z = sapply(t_val_all, function(y) {
      y[,which(as.Date(colnames(y),format = "%Y-%m-%d") == x)]
    })
  })
  names(t_val_all) <- main_cp

  # help functions
  most_frequent_with_ties <- function(x) {
    uniq_vals <- unique(x)
    counts <- table(x)
    max_counts <- max(counts)
    most_freq_vals <- names(counts)[counts == max_counts]
    return(paste(most_freq_vals, collapse=", "))
  }

  select_final_c <- function(x) {

    if (isTRUE(grepl(",",x))) {
      list_of_elements <- as.numeric(unlist((strsplit(x, ","))))
      prob <- c(0.18225,0.010125,0.010125,0.0005625,0.010125,0.0005625,0.0005625,
                0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,
                0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.0005625,
                0.010125,0.18225,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
                0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.0405,0.00225,
                0.000125,0.000125)
      y = list_of_elements[which.max(prob[list_of_elements])]
    } else {
      y = x
    }
    return(y)
  }
  max_freq <- function(x) {
    uniq_vals <- unique(x)
    counts <- table(x)
    max_counts <- round(max(counts)/length(x), digits = 2)
    return(max_counts)
  }

  Predict_ind_result <- lapply(as.character(main_cp), function(x) {

    new_df = na.omit(as.data.frame(t(t_val_all[[x]])))

    if (nrow(new_df) >0) {
      prediction = Prediction_CP(new_df)
      rownames(prediction) <- rownames(new_df)
    } else{
      prediction = NA
    }
    return(prediction)
  })
  names(Predict_ind_result) <- main_cp

  #####
  # Aggregation
  if (is.null(nearby_weight)) {
    final_config = lapply (Predict_ind_result, function(x) {

      if (all(is.na(x))) {
        z = NA
      } else {
        y = as.numeric(x$final_config)
        z = data.frame(max_freq_config = most_frequent_with_ties(y),
                       freq = max_freq(y)) %>%
          dplyr:: mutate(final_config = select_final_c(max_freq_config))

      }
    })
  }

  Res <- list(Noise_model = Noise_model,
              Test_result = t_val_all,
              Prediction_all = Predict_ind_result,
              Prediction_agg = final_config)

  return(Res)
}
