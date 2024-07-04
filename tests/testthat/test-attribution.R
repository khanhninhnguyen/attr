rm(list=ls())
# devtools::load_all()

# load data and make them in a good form

predictive_models = readRDS(file = system.file("extdata", "Rf.RDS",
                                               package = "attr"))

dataset <- get(load(file = system.file("extdata", "test.RData",
                                       package = "attr")))

dataset <- lapply(dataset, function(x) {
  names(x) <- c("Date", "GE", "GGp", "GEp", "EEp", "GpEp", "GpE")
  return(x)
})

# add change-points to test

main_cp = as.Date( c("2016-04-04"),
                      format = "%Y-%m-%d")
nearby_cp = list(
  bces = as.Date( c("2010-02-25",
                    "2012-08-13",
                    "2017-06-19"),
                  format = "%Y-%m-%d"),
  bcns = as.Date( c("2010-02-24",
                    "2011-12-17",
                    "2012-08-18",
                    "2017-06-29"),
                  format = "%Y-%m-%d"),
  pgc5 =  as.Date( c( "2009-07-24",
                      "2011-12-17",
                      "2017-11-20"),
                   format = "%Y-%m-%d"),
  ptaa = as.Date( c(),
                  format = "%Y-%m-%d"),
  p435 = as.Date( c("2011-10-12"),
                  format = "%Y-%m-%d")
)


# Premilinary test --------------------------------------------------------
## Test if is there any cluster (group of change points within a period e.g, 80 days)
## in the main and nearby station

### when main station has no cluster
main_cp0a = as.Date( c("2016-04-04", "2017-04-04"),
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
## station

nearby_cp0a = list(
  bcns = as.Date( c("2010-02-24",
                    "2011-12-17",
                    "2012-08-18",
                    "2017-06-29"),
                  format = "%Y-%m-%d")
  )

nearby_cp0b = list(
  bcns = as.Date( c("2010-02-24",
                    "2011-12-17",
                    "2012-08-18",
                    "2016-04-10"),
                  format = "%Y-%m-%d")
  )

test_that("Check similar change points", {
  expect_length(check_similar_CP(Six_Series = dataset[[1]],
                                 main_cp = main_cp,
                                 nearby_cp_one = nearby_cp0a[[1]],
                                 threshold = 10), 0)

  expect_identical(unlist(check_similar_CP(Six_Series = dataset[[1]],
                                           main_cp = main_cp,
                                           nearby_cp_one = nearby_cp0b[[1]],
                                           threshold = 10)),
                   as.Date("2016-04-10", format = "%Y-%m-%d"))
})


# Test Model Identification -----------------------------------------------


exp_NoiseModel = data.frame(
  bces = c(rep("ARMA(1,1)", 3), "AR(1)", rep("ARMA(1,1)", 2)),
  ptaa = c("ARMA(1,1)", rep("AR(1)", 5)),
  row.names = c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE")
)

test_that("Check model identification", {
  expect_identical(NoiseModel_Id(dataset = dataset[c("bces", "ptaa")],
                                 main_cp = main_cp,
                                 nearby_cp = nearby_cp[c("bces", "ptaa")]),
                   exp_NoiseModel)
})

# # First test on the clear case
# ##with given noise model
# # test0a = Attribution_CP(dataset,
# #                        main_break,
# #                        nearby_break,
# #                        noise_model_fix = "AR(1)")
# #
# # ## identify model by itself
# # test0b = Attribution_CP(dataset,
# #                              main_break,
# #                              nearby_break,
# #                              noise_model_fix = NULL)
# #
# # # test 1 : there is cluster in main, get error
# # main_break1 = c(main_break, main_break +10)
# # test1 = Attribution_CP(dataset,
# #                              main_break =  main_break1,
# #                              nearby_break)
# #
# # # test 2 : there is cluster in nearby, get error
# #
# # nearby_break1 <- nearby_break
# # nearby_break1$bces <- c(nearby_break1$bces, nearby_break1$bces[1] +10)
# # test2 = Attribution_CP(dataset,
# #                              main_break = main_break,
# #                              nearby_break = nearby_break1)
# #
# #
# #
