#' Apply the predictive rule to the test result
#'
#' @param test_result A data frame with dimension (n*6) where  where n is the
#' number of nearby stations. The 6 columns correspond to the 6 difference series.
#'
#' @param predictive_models A list of 24 Random Forest pre-trained and given in
#' "inst/extdata.Rf.RDS"
#'
#' @return A dataframe of 3 columns the  candidate configuration with the highest
#' frequency after running 24 models, the probability of the candidate, and the
#' final configuration. Each row corresponds to a nearby station.
#'
#' @import randomForest
#'
#' @importFrom dplyr mutate
#'
#' @importFrom stats predict

#' @export
#' @keywords internal
#'
predictive_rule <- function(test_result, predictive_models){

  Freq <- max_freq_config <- NULL

  stopifnot("test_result must be a data frame" = is.data.frame(test_result))

  stopifnot("test_result must have 6 columns" = (ncol(test_result)==6))

  condition_name = identical(names(test_result),
                             c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE"))

  stopifnot("Columns in test_result must be \
            : GE, GGp, GEp, EEp, GpEp, GpE" = isTRUE(condition_name))

  candidate_config <- data.frame(matrix(NA,
                                        ncol = length(predictive_models),
                                        nrow = nrow(test_result)))

  for (i in c(1:length(predictive_models))) {
    set.seed(1)
    Model_i = predictive_models[[paste0("model",i)]]
    candidate_config[,i] = predict(Model_i, newdata = test_result)

  }

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

  Res = candidate_config %>%
    mutate(max_freq_config = apply(candidate_config, 1, most_frequent_with_ties),
           proba = apply(candidate_config, 1, max_freq),
           final_config = sapply(max_freq_config, select_final_c))

  return(Res)
}
