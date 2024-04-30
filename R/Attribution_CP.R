#' Predict the origin of breakpoints in the main station
#'
#' @description
#' This function predicts the origin of each breakpoint in the main station by testing
#' the significance of jumps across six difference series for each main-nearby station
#' pair. It then predicts the configuration using a predictive rule and aggregates the
#' results for each breakpoint from all nearby stations.
#'
#' @details
#' The function has four main steps:
#'
#' 1. Model identification: Identify the noise model for each of the six difference
#' series between the main station and a nearby station.
#'
#' 2. Significance testing: Test the significance of jumps in each of the six difference
#' series using the identified noise model.
#'
#' 3. Configuration prediction: Predict the configuration of the main station and the
#' nearby station based on the test results.
#'
#' 4. Aggregation: Aggregate the prediction results for each breakpoint from all nearby
#' stations.
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
#' @param noise_model_fix A list of noise models, either fixed or to be identified. If
#' fixed, each model should be  of the following four: 'AR(1)', 'MA(1)', 'ARMA(1,1)',
#' or 'White'. If not fixed, the noise model will be identified using the longest
#' segment.
#'
#' @param nearby_weight A vector of numeric weights for each nearby station, used in
#' aggregating the final prediction. If not fixed, all nearby stations will be assigned
#' equal weight. The length of the vector should be the same as the
#' number of nearby stations.
#' @param limit_2side An integer specifying the number of points to be used in the
#' significance testing before and after each breakpoint. This can be used to speed up
#' the testing process. Default is NULL, which means all points are used.
#'
#' @return A list with four compnts:
#'
#' \describe{
#'   \item{noise_model}{A data frame of the identified noise models for each of the six
#' difference series between the main station and each nearby station. Each element is
#' a string indicating the noise model (e.g, 'AR(1)', 'MA(1)', etc)}
#'   \item{test_results}{A list of data frames, where each data frame contains the test
#' results for each breakpoint and each nearby station. The columns are the names of the
#' nearby stations, and the rows are the six difference series.}
#'   \item{prediction_results}{A list of data frames containing the predicted
#' configuration for each breakpoint and nearby station. Each data frame is named
#' after the corresponding nearby station and contains three columns: the
#' candidate configuration with the highest frequency after running 24 models,
#' the probability of the candidate, and the final configuration. Each row
#' corresponds to a nearby station.}
#'   \item{aggregated_results}{ A data frame of the aggregated prediction results
#' for each breakpoint, based on the weights specified in nearby_weight. This data
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
                           main_break,
                           nearby_break,
                           noise_model_fix = NULL,
                           nearby_weight = NULL,
                           limit_2side = 100){
  #####
  # pre-check
  Date <- candidate_config <- max_freq_config <-  NULL

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

  # check dataset include 6 series for each nearby
  list_six_diff = c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE")

  required_columns = lapply(c(1:length(dataset)), function(i) all(list_six_diff %in% names(dataset[[i]])))
  stopifnot("Columns in test_result must include \
            : Date, GE, GGp, GEp, EEp, GpEp, GpE" = all(required_columns))

  nearby_name <- names(dataset)
  #####
  # check cluster

  cluster_main = check_cluster_CP(Series_df = dataset[[1]],
                               Name_series = "GE",
                               Break_points = main_break)

  cluster_other = lapply(nearby_name, function(x){
    cluster_nearby = check_cluster_CP(Series_df = dataset[[x]],
                                   Name_series = "GpEp",
                                   Break_points = nearby_break[[x]])
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
    Noise_model <- NoiseModel_Id(dataset, main_break, nearby_break)
  } else {
    Noise_model <- data.frame(matrix(noise_model_fix,
                                     ncol = length(nearby_name),
                                     nrow = 6))
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
                     List_break_main = main_break,
                     List_break_nearby = nearby_break[[x]],
                     threshold = 10)
  })

  check_similar <- sapply(similar_cp, function(x) length(x)>0)

  nearby_break_m = nearby_break

  if(any(check_similar)) {
    ind_change = which(check_similar)
    for (i in ind_change) {
      list_brp_nb <- nearby_break[[nearby_name[i]]]
      nearby_break_m[[nearby_name[i]]] <- list_brp_nb [ which(!list_brp_nb %in% similar_cp [[i]])]
    }
  }

  #####
  # Significance test for the jump in the main
  test_sig <- function(List_break, Breakpoint, df, Name_series, noise_model, limit_period) {

    begin <-  List_break[tail(which(List_break < Breakpoint),1)]
    end <- List_break[which( List_break > Breakpoint)[1]]

    if (length(begin) == 0 | length(end) ==0) {

      t_val_jump <- NA

    } else if (length(begin) > 0 & length(end) > 0) {

      data_test <- df %>%
        dplyr:: filter(Date > begin & Date <= end)
      fit_fgls <- Test_CP(Series_df = data_test,
                                   Name_series = Name_series,
                                   Break_point = Breakpoint,
                                   noise_model = noise_model,
                                   limit = limit_period,
                                   name_case = NULL)

      t_val <- fit_fgls$summary$`t value`
      t_val_jump <- t_val[which(rownames(fit_fgls$summary) == "right")]

    }
    return(t_val_jump)
  }

  List_break_main <- sort( c( main_break,
                              get_min_max_date(data = dataset[[1]],
                                               column_name = "GE")))

  t_values_main = sapply(main_break, function(x) {
    t_val <- NA
    tryCatch({
      t_val <- test_sig(List_break = List_break_main,
                        Breakpoint = x,
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
    nearby_break_z <- sort(c(nearby_break[[x]],
                             get_min_max_date(data = dataset[[x]],
                                              column_name = "GpEp")))

    List_joint = sort(c(List_break_main, nearby_break_z))

    Break_points_five_series <- list(GGp = List_joint,
                                     GEp = List_joint,
                                     EEp = List_joint,
                                     GpEp = nearby_break_z,
                                     GpE = List_joint
    )

    t_val_all_brk <- sapply(main_break, function(y){

      five_t_val <- sapply(list_six_diff[-1], function(z){

        noise_model_z = unlist(Noise_model_order[which(rownames(Noise_model_order) == z), x])

        tryCatch({
          # Attempt to run the problematic function
          t <- test_sig(List_break = Break_points_five_series[[z]],
                        Breakpoint = y,
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
    colnames(y) <- as.character(main_break)
    y
  })

  #####
  # Prediction
  t_val_all <- lapply(main_break, function(x) {
    z = sapply(t_val_all, function(y) {
      y[,which(as.Date(colnames(y),format = "%Y-%m-%d") == x)]
    })
  })
  names(t_val_all) <- main_break

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
      prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
                0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
                0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
                0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
                0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
                0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)
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

  Predict_ind_result <- lapply(as.character(main_break), function(x) {

    new_df = na.omit(as.data.frame(t(t_val_all[[x]])))

    if (nrow(new_df) >0) {
      prediction = Prediction_CP(new_df)[,c(25:27)]
      rownames(prediction) <- rownames(new_df)
    } else{
      prediction = NA
    }
    return(prediction)
  })
  names(Predict_ind_result) <- main_break

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
