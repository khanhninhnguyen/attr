# rm(list=ls())
# devtools::load_all()
# predictive_models = readRDS(file = system.file("extdata", "Rf.RDS",
#                                                package = "attr"))
#
# dataset <- get(load(file = system.file("extdata", "test.RData",
#                                        package = "attr")))
# a <- NULL
#
# main_break = as.Date( c("2016-04-04"),
#                       format = "%Y-%m-%d")
# nearby_break = list(
#   bces = as.Date( c("2010-02-25",
#                     "2012-08-13",
#                     "2017-06-19"),
#                   format = "%Y-%m-%d"),
#   bcns = as.Date( c("2010-02-24",
#                     "2011-12-17",
#                     "2012-08-18",
#                     "2017-06-29"),
#                   format = "%Y-%m-%d"),
#   pgc5 =  as.Date( c( "2009-07-24",
#                       "2011-12-17",
#                       "2017-11-20"),
#                    format = "%Y-%m-%d"),
#   ptaa = as.Date( c(),
#                   format = "%Y-%m-%d"),
#   p435 = as.Date( c("2011-10-12"),
#                   format = "%Y-%m-%d")
# )
#
# # First test on the clear case
# ##with given noise model
# test0a = attribute_breakpoint(dataset,
#                              main_break,
#                              nearby_break,
#                              noise_model_fix = "AR(1)")
#
# ## identify model by itself
# test0b = attribute_breakpoint(dataset,
#                              main_break,
#                              nearby_break,
#                              noise_model_fix = NULL)
#
# # test 1 : there is cluster in main, get error
# main_break1 = c(main_break, main_break +10)
# test1 = attribute_breakpoint(dataset,
#                              main_break =  main_break1,
#                              nearby_break)
#
# # test 2 : there is cluster in nearby, get error
#
# nearby_break1 <- nearby_break
# nearby_break1$bces <- c(nearby_break1$bces, nearby_break1$bces[1] +10)
# test2 = attribute_breakpoint(dataset,
#                              main_break = main_break,
#                              nearby_break = nearby_break1)
#
#
#
