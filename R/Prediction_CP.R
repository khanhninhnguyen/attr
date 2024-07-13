#' Apply the predictive rule to the test result
#'
#' @description
#' This function applies trained predictive rules, built from the Random Forest algorithm,
#' on test results for a given change-point of a given station and its available nearby stations.
#' 24 predictive rules are used and stored in the "extdata/PredictiveRule.RDS" file.
#' They are trained with the IWV converted from the Nevada Geodetic Laboratory solution and ERA5,
#' and all lead to same and smallest classification error
#'
#' @param test_result a data frame with size n x 6 containing the t-values of the existence test
#' of the jump at a given change-point for all main-nearby station pairs with a given main station and its n nearby stations.
#'
#' @return a data frame with size n x 3 containing the predicted configurations for the given change-point
#' and the n nearby stations.The three columns are: \code{max_freq_config} is the candidate
#' configuration with the highest frequency from running the 24 predictive rules, \code{freq}
#' is the associated frequency of this candidate, and \code{final_config_one} is the final configuration, i.e.
#' \code{max_freq_config} if the \code{max_freq_config} is unique. If \code{max_freq_config} is not unique,
#' the \code{final_config} is the one with the highest prior probability of the table given in the article
#'
#'
#' @import randomForest
#' @importFrom dplyr mutate
#' @importFrom stats predict
#'
#' @export
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

  Res = data.frame(max_freq_config = apply(candidate_config, 1, most_frequent_with_ties),
                   freq = apply(candidate_config, 1, max_freq),
                   final_config = sapply(max_freq_config, select_final_c))

  return(Res)

}
