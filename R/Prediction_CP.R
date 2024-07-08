#' Apply the predictive rule to the test result
#'
#' @description
#' This function applies the trained predictive rules using Random Forest
#' on the test results, which include 6 t-values for each change-point and 1 nearby
#' station.  The 24 predictive rules stored in the "extdata/PredictiveRule.RDS"
#' file are trained with the IWV converted from the Nevada Geodetic Laboratory solution
#' and ERA5 data.
#'
#' @param test_result A data frame contains t-values for changes in the mean of
#' six difference series, assessed for all main-nearby station pairs. The data frames
#' are structured with 6 rows for the six difference series and n columns corresponding
#' to the nearby stations.
#'
#' @return A data frame contains the predicted configurations for each change-point.
#' The data frames include three columns: \code{candidate_config} is the candidate
#' configuration with the highest frequency from running 24 models, \code{proba}
#' is the probability of the candidate, and \code{final_config} is the final configuration.
#' Each row corresponds to a nearby station.
#'
#' @import randomForest
#'
#' @importFrom dplyr mutate
#'
#' @importFrom stats predict

#' @export
#' @keywords internal
#'
Prediction_CP <- function(test_result){

  # load predictive models
  predictive_models = readRDS(file = system.file("extdata", "PredictiveRule.RDS", package = "attr"))

  Freq <- max_freq_config <- NULL

  # Validate the input form and column's names
  stopifnot("test_result must be a data frame" = is.data.frame(test_result))
  stopifnot("test_result must have 6 columns" = (ncol(test_result)==6))
  required_columns = all(c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE") %in% names(test_result))
  stopifnot("Columns in test_result must be : GE, GGp, GEp, EEp, GpEp, GpE" = isTRUE(required_columns))

  # Apply RF models
  candidate_config <- data.frame(matrix(NA,
                                        ncol = length(predictive_models),
                                        nrow = nrow(test_result)))

  for (i in c(1:length(predictive_models))) {
    set.seed(1)
    Model_i = predictive_models[[paste0("model",i)]]
    candidate_config[,i] = predict(Model_i, newdata = test_result)
  }

  # Helper function to find the most frequent values with ties
  most_frequent_with_ties <- function(x) {
    uniq_vals <- unique(x)
    counts <- table(x)
    max_counts <- max(counts)
    most_freq_vals <- names(counts)[counts == max_counts]
    return(paste(most_freq_vals, collapse=", "))
  }

  # Helper function to select the final configuration based on probabilities
  select_final_c <- function(x) {
    prob <- c(0.18225,0.010125,0.010125,0.0005625,0.010125,0.0005625,0.0005625,
              0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,
              0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.0005625,
              0.010125,0.18225,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
              0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.0405,0.00225,
              0.000125,0.000125)
    prob <- prob / sum(prob)

    if (isTRUE(grepl(",",x))) {
      list_of_elements <- as.numeric(unlist((strsplit(x, ","))))
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

  Res = candidate_config %>%
    mutate(max_freq_config = apply(candidate_config, 1, most_frequent_with_ties),
           proba = apply(candidate_config, 1, max_freq),
           final_config = sapply(max_freq_config, select_final_c))

  return(Res)

}
