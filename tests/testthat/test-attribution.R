rm(list=ls())
# devtools::load_all()
library(testthat)
# load data and make them in a good form

dataset <- readRDS(file = system.file("extdata", "dataset.RDS",
                                       package = "attr"))

# add change-points to test
convert_dates <- function(dates_list) {
  lapply(dates_list, function(dates) as.Date(dates, format = "%Y-%m-%d"))
}

main_cp <- as.Date("2016-04-04", format = "%Y-%m-%d")

nearby_cp <- convert_dates(list(
  bces = c("2010-02-25", "2012-08-13", "2017-06-19"),
  bcns = c("2010-02-24", "2011-12-17", "2012-08-18", "2017-06-29"),
  pgc5 = c("2009-07-24", "2011-12-17", "2017-11-20"),
  ptaa = c(),
  p435 = c("2011-10-12")
))

# Premilinary test --------------------------------------------------------
## Test if is there any cluster (group of change points within a period e.g, 80 days)
## in the main and nearby station

### when main station has no cluster
main_cp0a = as.Date( c("2016-04-04"),
                   format = "%Y-%m-%d")


### when main station has 1 cluster including the first 2 change points
main_cp0b = as.Date( c("2016-04-04","2016-03-04", "2017-04-04"),
                     format = "%Y-%m-%d")

### run test
test_that("Check cluster", {
  expect_identical(sum(check_cluster_CP(Series_df = dataset[[1]],
                                        Name_series = "GE",
                                        CP = main_cp0a)$cluster), 0)

  expect_identical(sum(check_cluster_CP(Series_df = dataset[[1]],
                                           Name_series = "GE",
                                           CP = main_cp0b)$cluster), 3)
})

## Test if any change point in the nearby is close to the change points in the main
## station, i.e.
#' \code{max_freq_config_one} if the \code{max_freq_config_one} is unique. If \code{max_freq_config_one} is not unique,
#' the \code{final_config_one} is the one with the highest prior probability of the table given in the article

nearby_cp0a = list(
  bcns = as.Date( c("2010-02-25",
                    "2011-12-17",
                    "2012-08-18",
                    "2017-06-29"),
                  format = "%Y-%m-%d")
  )

nearby_cp0b = list(
  bcns = as.Date( c("2010-02-25",
                    "2011-12-17",
                    "2012-08-18",
                    "2016-04-10"),
                  format = "%Y-%m-%d")
  )

test_that("Check similar change points", {
  expect_identical(check_similar_CP(Series_df = dataset[[1]],
                                 main_cp = main_cp,
                                 nearby_cp_one = nearby_cp0a[[1]],
                                 threshold = 10),
                   "No similar change-points")

  expect_identical(unlist(check_similar_CP(Series_df = dataset[[1]],
                                           main_cp = main_cp,
                                           nearby_cp_one = nearby_cp0b[[1]],
                                           threshold = 10)),
                   as.Date("2016-04-10", format = "%Y-%m-%d"))
})


# Test Model Identification -----------------------------------------------

exp_NoiseModel = data.frame(
  bces = c("ARMA(1,1)", "AR(1)", "ARMA(1,1)", "AR(1)", "ARMA(1,1)", "ARMA(1,1)"),
  ptaa = c("ARMA(1,1)", rep("AR(1)", 3), "ARMA(1,1)", "AR(1)"),
  row.names = c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE")
)

test_that("Check model identification", {
  expect_identical(NoiseModel_Id(dataset = dataset[c("bces", "ptaa")],
                                 main_cp = main_cp,
                                 nearby_cp = nearby_cp[c("bces", "ptaa")]),
                   exp_NoiseModel)
})


# Test the test of significance of change-point (Test_CP.R) ------------------

test_res = Test_CP(dataset$bces, Name_series = "GE", CP = main_cp[1],
                 noise_model = c("ARMA(1,1)"), limit = 100)

test_that("Check test of one change-point", {
  expect_identical(test_res$Summary_tab$`t value`[which(rownames(test_res$Summary_tab) == "jump")],
                   -2.58)
})


# Test the Attribution_CP -------------------------------------------------

## Case 1: fix 1 noise model for all series

test_Attribution_CP_c1 = Attribution_CP(dataset,
                        main_cp = main_cp,
                        nearby_cp = nearby_cp,
                        noise_model_fix = "AR(1)")

test_that("Attribution", {
  expect_identical(unique(as.vector(as.matrix(test_Attribution_CP_c1$Noise_model))),
                   "AR(1)")
  expect_identical(unlist(test_Attribution_CP_c1$Test_result[[as.character(main_cp)]])[4,4],
                   -0.386)
  expect_identical(as.character(test_Attribution_CP_c1$Prediction_agg[[as.character(main_cp)]])[1],
                   "23")
})

## Case 2: unknown noise model for all series
test_Attribution_CP_c2 = Attribution_CP(dataset,
                                        main_cp = main_cp,
                                        nearby_cp = nearby_cp,
                                        noise_model_fix = NULL)

test_that("Attribution", {
  expect_identical(test_Attribution_CP_c2$Noise_model[,c("bces", "ptaa")],
                   exp_NoiseModel)
  expect_identical(unlist(test_Attribution_CP_c2$Test_result[[as.character(main_cp)]])[6,4],
                   -1.7)
  expect_identical(as.character(test_Attribution_CP_c2$Prediction_agg[[as.character(main_cp)]])[1],
                   "23")
})

## Case 3: known noise model for each series
test_Attribution_CP_c3 = Attribution_CP(dataset[c("bces", "ptaa")],
                                        main_cp = main_cp,
                                        nearby_cp = nearby_cp[c("bces", "ptaa")],
                                        noise_model_fix = exp_NoiseModel)

test_that("Attribution", {
  expect_identical(test_Attribution_CP_c2$Noise_model[,c("bces", "ptaa")],
                   exp_NoiseModel)
  expect_identical(unlist(test_Attribution_CP_c2$Test_result[[as.character(main_cp)]])[6,4],
                   -1.7)
  expect_identical(as.character(test_Attribution_CP_c2$Prediction_agg[[as.character(main_cp)]])[1],
                   "23")
})
# Test prediction ---------------------------------------------------------

test_that("Prediction", {
  expect_error(Prediction_CP(c(2,2,2,1,1,1)),
               "test_result must be a data frame", fixed = TRUE)
  expect_error(Prediction_CP(data.frame(2,2,2,1,1,1)),
               "Columns in test_result must be : GE, GGp, GEp, EEp, GpEp, GpE",
               fixed = TRUE)
  expect_identical(
    as.character(
      Prediction_CP(data.frame(GE = 2, GGp = 2, GEp = 2, EEp = 1, GpEp = 1, GpE = 1))$final_config),
    "1")
})



