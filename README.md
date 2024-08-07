
<!-- README.md is generated from README.Rmd. Please edit that file -->

# attr

This function attributes change-points detected in the difference series
to one or both contributing series. Assuming the primary series is GPS
and the secondary is ERA, change-points are identified in the GPS - ERA
series. The function determines whether these change-points are due to
GPS, ERA, or both.

## Table of Contents

-   [Installation](#installation)
-   [Usage](#usage)
-   [Examples](#examples)

## Installation

You can install the development version of `attr` from GitHub with:

``` r
# Install devtools if you haven't already
# install.packages("devtools")

# Install the package from GitHub
# devtools::install_github("khanhninhnguyen/attr")
```

## Usage

The main function of this package is `Attribution_CP`. The method
utilizes data from a nearby station providing GPS’ and ERA’ series.
Users must compute six difference series (GPS-ERA, GPS-GPS’, GPS-ERA’,
ERA-ERA’, GPS’-ERA’, GPS’-ERA) for each main-nearby station pair to use
as input in the `dataset` variable. An example dataset is given in
“attr/inst/extdata/dataset.RData”:

``` r
library(attr)

data("dataset")

str(dataset) 
#> List of 5
#>  $ bces:'data.frame':    10309 obs. of  7 variables:
#>   ..$ date: Date[1:10309], format: "1994-01-02" "1994-01-03" ...
#>   ..$ GE  : num [1:10309] 0.33 1.07 -0.69 0.14 0.16 -0.86 -0.57 -0.31 -0.15 0.93 ...
#>   ..$ GGp : num [1:10309] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GEp : num [1:10309] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ EEp : num [1:10309] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpEp: num [1:10309] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpE : num [1:10309] NA NA NA NA NA NA NA NA NA NA ...
#>  $ bcns:'data.frame':    10320 obs. of  7 variables:
#>   ..$ date: Date[1:10320], format: "1994-01-02" "1994-01-03" ...
#>   ..$ GE  : num [1:10320] 0.33 1.07 -0.69 0.14 0.16 -0.86 -0.57 -0.31 -0.15 0.93 ...
#>   ..$ GGp : num [1:10320] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GEp : num [1:10320] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ EEp : num [1:10320] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpEp: num [1:10320] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpE : num [1:10320] NA NA NA NA NA NA NA NA NA NA ...
#>  $ pgc5:'data.frame':    10353 obs. of  7 variables:
#>   ..$ date: Date[1:10353], format: "1994-01-02" "1994-01-03" ...
#>   ..$ GE  : num [1:10353] 0.33 1.07 -0.69 0.14 0.16 -0.86 -0.57 -0.31 -0.15 0.93 ...
#>   ..$ GGp : num [1:10353] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GEp : num [1:10353] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ EEp : num [1:10353] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpEp: num [1:10353] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpE : num [1:10353] NA NA NA NA NA NA NA NA NA NA ...
#>  $ ptaa:'data.frame':    10313 obs. of  7 variables:
#>   ..$ date: Date[1:10313], format: "1994-01-02" "1994-01-03" ...
#>   ..$ GE  : num [1:10313] 0.33 1.07 -0.69 0.14 0.16 -0.86 -0.57 -0.31 -0.15 0.93 ...
#>   ..$ GGp : num [1:10313] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GEp : num [1:10313] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ EEp : num [1:10313] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpEp: num [1:10313] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpE : num [1:10313] NA NA NA NA NA NA NA NA NA NA ...
#>  $ p435:'data.frame':    10339 obs. of  7 variables:
#>   ..$ date: Date[1:10339], format: "1994-01-02" "1994-01-03" ...
#>   ..$ GE  : num [1:10339] 0.33 1.07 -0.69 0.14 0.16 -0.86 -0.57 -0.31 -0.15 0.93 ...
#>   ..$ GGp : num [1:10339] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GEp : num [1:10339] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ EEp : num [1:10339] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpEp: num [1:10339] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ GpE : num [1:10339] NA NA NA NA NA NA NA NA NA NA ...
```

For example, we want to attribute the change-points in the main station:

``` r
str(main_cp)
#>  Date[1:1], format: "2016-04-04"
```

The nearby stations have some change-points too:

``` r
str(nearby_cp)
#> List of 5
#>  $ bces: Date[1:3], format: "2010-02-25" "2012-08-13" ...
#>  $ bcns: Date[1:4], format: "2010-02-24" "2011-12-17" ...
#>  $ pgc5: Date[1:3], format: "2009-07-24" "2011-12-17" ...
#>  $ ptaa: 'Date' num(0) 
#>  $ p435: Date[1:1], format: "2011-10-12"
```

The attribution method involves three main tasks, each implemented in a
function that can run independently. Understanding these steps will help
you understand the `Attribution_CP` function :

### 1. Model Identification for All Time Series

The `NoiseModel_Id` function fits the noise model for six difference
series in `dataset` according to three ARMA models (AR(1), MA(1) and
ARMA(1,1)) and the white noise.

``` r
NoiseModel_Id(dataset = dataset,
              main_cp = main_cp,
              nearby_cp = nearby_cp)
#>           bces      bcns      pgc5      ptaa      p435
#> GE   ARMA(1,1) ARMA(1,1) ARMA(1,1) ARMA(1,1) ARMA(1,1)
#> GGp      AR(1)     AR(1)     AR(1)     AR(1)     AR(1)
#> GEp  ARMA(1,1) ARMA(1,1)     AR(1)     AR(1)     AR(1)
#> EEp      AR(1)     AR(1)     MA(1)     AR(1)     AR(1)
#> GpEp ARMA(1,1) ARMA(1,1)     AR(1) ARMA(1,1) ARMA(1,1)
#> GpE  ARMA(1,1)     AR(1)     AR(1)     AR(1)     AR(1)
```

### 2. Test of the Jump at Change-Points in Every Series

The `Test_CP` function tests the jump at change-points in a series by
calling its name in the elementary dataframe in `dataset` :

``` r
TestRes = Test_CP(dataset$bces, Name_series = "GE", CP = main_cp[1],
                 noise_model = c("ARMA(1,1)"), limit = 100)
str(TestRes)
#> List of 2
#>  $ Summary_tab:'data.frame': 10 obs. of  4 variables:
#>   ..$ Estimate  : num [1:10] 0.23 -1.46 -0.758 -0.792 -0.363 0.222 -0.134 0.109 -0.823 0.854
#>   ..$ Std. Error: num [1:10] 1.24 7.97 5.14 1.63 1.11 2.32 0.617 0.401 0.319 4.57
#>   ..$ t value   : num [1:10] 0.185 -0.183 -0.147 -0.487 -0.326 0.0955 -0.216 0.271 -2.58 0.187
#>   ..$ Pr(>|t|)  : num [1:10] 0.853 0.855 0.883 0.626 0.745 0.924 0.829 0.786 0.0099 0.852
#>  $ Used_dates : Date[1:200], format: "2015-12-27" "2015-12-28" ...
```

### 3. Configuration Prediction

`Prediction_CP` function runs predictions by applying predictors trained
from our large dataset. The input is a dataframe of t-values for six
series.

``` r
# arbitrary example: assume that the t values of 6 tests are : 2,2,2,1,1,1
ArbPred = Prediction_CP(data.frame(GE = 2, GGp = 2, GEp = 2, EEp = 1, GpEp = 1, GpE = 1))
ArbPred
#>   max_freq_config freq final_config
#> 1               1    1            1
```

### Examples

Finally, the main function includes all these steps. The noise model can
be either fixed for all time series :

    #> List of 4
    #>  $ Noise_model   :'data.frame':  6 obs. of  5 variables:
    #>   ..$ bces: chr [1:6] "AR(1)" "AR(1)" "AR(1)" "AR(1)" ...
    #>   ..$ bcns: chr [1:6] "AR(1)" "AR(1)" "AR(1)" "AR(1)" ...
    #>   ..$ pgc5: chr [1:6] "AR(1)" "AR(1)" "AR(1)" "AR(1)" ...
    #>   ..$ ptaa: chr [1:6] "AR(1)" "AR(1)" "AR(1)" "AR(1)" ...
    #>   ..$ p435: chr [1:6] "AR(1)" "AR(1)" "AR(1)" "AR(1)" ...
    #>  $ Test_result   :List of 1
    #>   ..$ 2016-04-04: num [1:6, 1:5] -2.5 0.233 -2.14 2.13 -2.16 -2.52 -2.5 -0.206 -1.62 0.508 ...
    #>   .. ..- attr(*, "dimnames")=List of 2
    #>   .. .. ..$ : chr [1:6] "GE" "GGp" "GEp" "EEp" ...
    #>   .. .. ..$ : chr [1:5] "bces" "bcns" "pgc5" "ptaa" ...
    #>  $ Prediction_all:List of 1
    #>   ..$ 2016-04-04:'data.frame':   5 obs. of  3 variables:
    #>   .. ..$ max_freq_config: chr [1:5] "38" "23" "15" "23" ...
    #>   .. ..$ freq           : num [1:5] 1 0.62 0.58 1 1
    #>   .. ..$ final_config   : Named chr [1:5] "38" "23" "15" "23" ...
    #>   .. .. ..- attr(*, "names")= chr [1:5] "38" "23" "15" "23" ...
    #>  $ Prediction_agg:List of 1
    #>   ..$ 2016-04-04:'data.frame':   1 obs. of  3 variables:
    #>   .. ..$ max_freq_config: chr "23"
    #>   .. ..$ freq           : num 0.6
    #>   .. ..$ final_config   : chr "23"

Or specified for all time series:

``` r
NoiseModel = data.frame(
  bces = c("ARMA(1,1)", "AR(1)", "ARMA(1,1)", "AR(1)", "ARMA(1,1)", "ARMA(1,1)"),
  ptaa = c("ARMA(1,1)", rep("AR(1)", 3), "ARMA(1,1)", "AR(1)"),
  row.names = c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE")
)
Attribution_CP( dataset = dataset[c("bces", "ptaa")],
                main_cp = main_cp,
                nearby_cp = nearby_cp[c("bces", "ptaa")],
                noise_model_fix = NoiseModel)
#> $Noise_model
#>           bces      ptaa
#> GE   ARMA(1,1) ARMA(1,1)
#> GGp      AR(1)     AR(1)
#> GEp  ARMA(1,1)     AR(1)
#> EEp      AR(1)     AR(1)
#> GpEp ARMA(1,1) ARMA(1,1)
#> GpE  ARMA(1,1)     AR(1)
#> 
#> $Test_result
#> $Test_result$`2016-04-04`
#>        bces   ptaa
#> GE   -2.580 -2.580
#> GGp   0.233 -0.484
#> GEp  -2.170 -2.150
#> EEp   2.130 -0.386
#> GpEp -2.190 -2.480
#> GpE  -2.570 -1.700
#> 
#> 
#> $Prediction_all
#> $Prediction_all$`2016-04-04`
#>      max_freq_config freq final_config
#> bces              38    1           38
#> ptaa              23    1           23
#> 
#> 
#> $Prediction_agg
#> $Prediction_agg$`2016-04-04`
#>   max_freq_config freq final_config
#> 1          23, 38  0.5           23
```

Or set to `NULL`, and the noise models of all time series will be
estimated using the `NoiseModel_Id` function :

``` r
Attribution_CP( dataset,
                main_cp = main_cp,
                nearby_cp = nearby_cp,
                noise_model_fix = NULL)
#> $Noise_model
#>           bces      bcns      pgc5      ptaa      p435
#> GE   ARMA(1,1) ARMA(1,1) ARMA(1,1) ARMA(1,1) ARMA(1,1)
#> GGp      AR(1)     AR(1)     AR(1)     AR(1)     AR(1)
#> GEp  ARMA(1,1) ARMA(1,1)     AR(1)     AR(1)     AR(1)
#> EEp      AR(1)     AR(1)     MA(1)     AR(1)     AR(1)
#> GpEp ARMA(1,1) ARMA(1,1)     AR(1) ARMA(1,1) ARMA(1,1)
#> GpE  ARMA(1,1)     AR(1)     AR(1)     AR(1)     AR(1)
#> 
#> $Test_result
#> $Test_result$`2016-04-04`
#>        bces   bcns   pgc5   ptaa   p435
#> GE   -2.580 -2.580 -2.580 -2.580 -2.580
#> GGp   0.233 -0.206 -0.439 -0.484  0.145
#> GEp  -2.170 -1.690 -1.610 -2.150 -1.620
#> EEp   2.130  0.508  0.362 -0.386 -0.123
#> GpEp -2.190 -1.990 -1.720 -2.480 -2.370
#> GpE  -2.570 -1.970 -1.610 -1.700 -1.640
#> 
#> 
#> $Prediction_all
#> $Prediction_all$`2016-04-04`
#>      max_freq_config freq final_config
#> bces              38 1.00           38
#> bcns              23 1.00           23
#> pgc5              15 0.58           15
#> ptaa              23 1.00           23
#> p435              23 1.00           23
#> 
#> 
#> $Prediction_agg
#> $Prediction_agg$`2016-04-04`
#>   max_freq_config freq final_config
#> 1              23  0.6           23
```
