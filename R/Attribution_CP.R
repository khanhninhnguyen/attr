#' Predict the Origin of Change-Points
#'
#' @description
#' This function predicts the origin of change-points in a station, which detected
#' on the difference series GPS - reference (ERA). The ideas of this method is to use
#' nearby station, which give GPS' and ERA' series. User need to compute the
#' six difference series  (GPS-ERA, GPS-GPS', GPS-ERA', ERA-ERA', GPS'-ERA', GPS'-ERA)
#' for each main-nearby station pair. This function tests the significance of changes
#' in mean in each of the six difference series for each main-nearby station pair,
#' predicts configurations from the test results, and aggregates them from all nearby stations
#' for each change-points.
#'
#' @details
#' The function operates in four main steps:
#'
#' \strong{1. Model Identification} (optional): Uses the `NoiseModel_Id.R` function to
#' to identify the noise model for each difference series.
#'
#' \strong{2. Significance Test}: Uses the `Test_CP.R` function to test the significance
#' of changes in mean in each difference series.
#'
#' \strong{3. Configuration Prediction}:
#'
#'    - Uses the `Prediction_CP.R` function to predict the configuration of each
#' triplet (main, changepoint, nearby) based on their test result of the six difference series.
#'
#'    - Uses 24 pre-trained predictive models from `Predictiverule.RDS`.
#'
#' \strong{4. Aggregation}: Aggregates the prediction results for each change-point
#' from all nearby stations.
#'
#' Additionally, the function performs two preliminary checks:
#'
#' \strong{a. Cluster Check}:
#'
#'    - Looks for clusters of change-points (within 80 days) in main and nearby stations.
#'
#'    - If a cluster exists, the function stops.
#'
#' \strong{b. Coincident Check}:
#' \itemize{
#'   \item Identifies "coincident" change-points (within 10 days) between main and nearby stations.
#'   \item Coincident points in nearby stations are not used for selection of segment for significance test.
#'   \item Data between main change-points and coincident points is removed from four non-collocated series.
#' }
#'
#' @param dataset A list of length n where n is the number of nearby station.
#' Each element of the list is a data frame representing a main-nearby
#' station pair and is named after the nearby station. The data frames have 7 columns.
#' The first column is labeled "Date" and contains dates in "YYYY-MM-DD" format.
#' The other 6 columns contain numeric values. These columns represent six difference
#' series: GPS - ERA, GPS - GPS', GPS - ERA', ERA - ERA', GPS' - ERA', and GPS' - ERA.
#'
#' @param main_cp A vector of day representing change-points in the main station
#' in Date format: "\%Y-\%m-\%d".
#'
#' @param nearby_cp A list of length n where n is the number of nearby station.
#' Each element is a vector contains of day representing change-points for the
#' corresponding nearby station in Date format: "\%Y-\%m-\%d".
#'
#'   - Each vector named after the corresponding station.
#' @param noise_model_fix Can be:
#'
#'   - A string specifying one noise model for all series.
#'
#'   - A data frame (6 rows x n columns) specifying models for each series (row) and nearby station (column).
#'
#'   - NULL to determine noise models using NoiseModel_Id.R.
#'   Valid models: 'AR(1)', 'MA(1)', 'ARMA(1,1)', 'White'.
#'
#' @param nearby_weight A numeric vector where:
#'
#'   - Length equals the number of nearby stations
#'
#'   - Each element is the weight for a nearby station
#'
#'   - If NULL, all stations are assigned equal weight
#'
#' @param limit_2side An integer that:
#'
#'   - Specifies the number of points used before and after each changepoint for testing
#'
#'   - If NULL (default), all points are used
#'
#' Please verify the format of input in the tests/testthat/test-attribution.R
#' CHECK HOW TO MAKE IT LINE BY LINE
#' @return A list with four components:
#'
#' \describe{
#'   \item{noise_model}{
#'A data frame (6 rows x n columns) where:
#'
#'  * Rows represent the six difference series
#'
#'  * Columns represent nearby stations
#'
#'  * Each cell contains a string indicating the noise model (e.g., 'AR(1)',
#' MA(1)') of the according series }
#'
#'   \item{test_results}{A list of data frames, where each data frame contains the test
#' results for each changepoint and each nearby station. The columns are the names of the
#' nearby stations, and the rows are the six difference series.}
#'   \item{prediction_results}{A list of data frames containing the predicted
#' configuration for each changepoint and nearby station. Each data frame is named
#' after the corresponding nearby station and contains three columns: the
#' candidate configuration with the highest frequency after running 24 models,
#' the probability of the candidate, and the final configuration. Each row
#' corresponds to a nearby station.}
#'   \item{aggregated_results}{ A data frame of the aggregated prediction results
#' for each changepoint, based on the weights specified in nearby_weight. This data
#' frame is presented in a similar format to the prediction_results data frame for
#' a single case}
#' }
#'
#'
#' @importFrom dplyr %>% mutate select rename any_of
#'
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
                           limit_2side = 100){
  #####
  # pre-check
  Date <- candidate_config <- max_freq_config <-  NULL

  stopifnot("dataset must be a list of dataframe" = is.list(dataset))

  stopifnot("dataset must be a list of dataframe" = all(sapply(dataset, is.data.frame)))

  stopifnot("nearby_cp must be a list" = is.list(nearby_cp))

  stopifnot("names of dataset must be the same to nearby_cp" = setequal(names(nearby_cp), names(dataset)))

  stopifnot("dataset and nearby_cp must have the same length" = (length(dataset) == length(nearby_cp)))


  if (!inherits( main_cp, "Date") | any(!sapply(nearby_cp, inherits, "Date"))) {
    stop("Must specify both main_cp and nearby_cp in Date type")
  }

  if( any(!sapply(dataset, function(df) {
    ncol(df) == 7 & inherits(df[,1], "Date")
  }))) {
    stop("Each data frame in the dataset need to have 7 columns &\
         the first column must be in Date type")

  }

  # check dataset include 6 series for each nearby
  list_six_diff = c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE")

  required_columns = sapply(c(1:length(dataset)), function(i) all(list_six_diff %in% names(dataset[[i]])))
  stopifnot("Columns in test_result must include \
            : Date, GE, GGp, GEp, EEp, GpEp, GpE" = all(required_columns))

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

  # convert to model order
  convert_model_order <- function(model_name) {
    if(model_name == "AR(1)") {
      y = c(1,0,0)
    } else if(model_name == "MA(1)") {
      y = c(0,0,1)
    } else if(model_name == "ARMA(1,1)") {
      y = c(1,0,1)
    } else if(model_name == "White") {
      y = c(0,0,0)
    }
    return(y)
  }
  Noise_model_order = as.data.frame(sapply(nearby_name, function(x) {
    lapply(Noise_model[[x]], convert_model_order)
  }))
  rownames(Noise_model_order) <- list_six_diff

  #####
  # remove the similar change points in the nearby for the test

  similar_cp <- lapply(nearby_name, function(x) {
    check_similar_CP(Six_Series = dataset[[x]],
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
  # Significance test for the jump in the main
  test_sig <- function(List_CP, CP, df, Name_series, noise_model, limit_period) {

    begin <-  List_CP[tail(which(List_CP < CP),1)]
    end <- List_CP[which( List_CP > CP)[1]]

    if (length(begin) == 0 | length(end) ==0) {

      t_val_jump <- NA

    } else if (length(begin) > 0 & length(end) > 0) {

      data_test <- df %>%
        dplyr:: filter(Date > begin & Date <= end)
      fit_fgls <- Test_CP(Series_df = data_test,
                                   Name_series = Name_series,
                                   CP = CP,
                                   noise_model = noise_model,
                                   limit = limit_period,
                                   name_case = NULL)

      t_val <- fit_fgls$summary$`t value`
      t_val_jump <- t_val[which(rownames(fit_fgls$summary) == "right")]

    }
    return(t_val_jump)
  }

  List_CP_main <- sort( c( main_cp,
                              get_min_max_date(data = dataset[[1]],
                                               column_name = "GE")))

  t_values_main = sapply(main_cp, function(x) {
    t_val <- NA
    tryCatch({
      t_val <- test_sig(List_CP = List_CP_main,
                        CP = x,
                        df = dataset[[1]],
                        Name_series = "GE",
                        noise_model = unlist(Noise_model_order[1,1]),
                        limit_period = limit_2side)

    }, error = function(e) {

    })
    return(t_val)
  })
  # Significance test for the other 5 test
  t_values_other = lapply(nearby_name, function(x) {
    t <- NA
    nearby_cp_z <- sort(c(nearby_cp[[x]],
                             get_min_max_date(data = dataset[[x]],
                                              column_name = "GpEp")))

    List_joint = sort(c(List_CP_main, nearby_cp_z))

    CP_five_series <- list(GGp = List_joint,
                           GEp = List_joint,
                           EEp = List_joint,
                           GpEp = nearby_cp_z,
                           GpE = List_joint
    )

    t_val_all_brk <- sapply(main_cp, function(y){

      five_t_val <- sapply(list_six_diff[-1], function(z){

        noise_model_z = unlist(Noise_model_order[which(rownames(Noise_model_order) == z), x])

        tryCatch({
          # Attempt to run the problematic function
          t <- test_sig(List_CP = CP_five_series[[z]],
                        CP = y,
                        df = dataset[[x]],
                        Name_series = z,
                        noise_model = noise_model_z,
                        limit_period = limit_2side)
        }, error = function(e) {
          # If there's an error, return NA
          t <- NA
        })
        setNames(t, z)
        return(t)
      })

      return(five_t_val)

    })
    return(t_val_all_brk)

  })
  names(t_values_other) <- nearby_name

  t_val_all <- lapply(t_values_other, function(x) {
    new_row <- matrix(t_values_main, nrow = 1) # Create a new row of zeros (change this as needed)
    y <- rbind(new_row, x) # Bind the new row at the top
    rownames(y)[1] <- "GE"
    colnames(y) <- as.character(main_cp)
    y
  })

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
      prediction = Prediction_CP(new_df)[,c(25:27)]
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
        z = data.frame(candidate_config = most_frequent_with_ties(y),
                       proba = max_freq(y)) %>%
          mutate(final_config = select_final_c(candidate_config))

      }
    })
  }

  Res <- list(Noise_model = Noise_model,
              Test_result = t_val_all,
              Prediction_all = Predict_ind_result,
              Prediction_agg = final_config)

  return(Res)
}
