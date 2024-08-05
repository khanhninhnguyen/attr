---
title: "README.md"
author: "Ninh"
date: "02/08/2024"
output: html_document
---

# attr

This function attributes change-points detected in the difference series to one or both contributing series. Assuming the primary series is GPS and the secondary is ERA, change-points are identified in the GPS - ERA series. 
The function determines whether these change-points are due to GPS, ERA, or both.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)

## Installation

You can install the development version of `attr` from GitHub with:

```{r Install, echo=FALSE}

# Install devtools if you haven't already
# install.packages("devtools")

# Install the package from GitHub
# devtools::install_github("khanhninhnguyen/attr")
library(attr)
```

## Usage 

The main function of this package is `Attribution_CP`. The method utilizes data from a nearby station providing GPS' and ERA' series. Users must compute six difference series (GPS-ERA, GPS-GPS', GPS-ERA', ERA-ERA', GPS'-ERA', GPS'-ERA) for each main-nearby station pair to use as input in the `dataset` variable. An example dataset is given in "attr/inst/extdata/dataset.RData":

```{r dataset, echo=FALSE}
dataset <- readRDS(file = system.file("extdata", "dataset.RDS",
                                       package = "attr"))
```

```{r}
str(dataset) 
```

For example, we want to attribute the change-points in the main station:

```{r, echo=FALSE}
convert_dates <- function(dates_list) {
  lapply(dates_list, function(dates) as.Date(dates, format = "%Y-%m-%d"))
}

main_cp <- as.Date("2016-04-04", format = "%Y-%m-%d")
```

```{r}
str(main_cp)
```
The nearby stations have some change-points too:

```{r, echo=FALSE}
nearby_cp <- convert_dates(list(
  bces = c("2010-02-25", "2012-08-13", "2017-06-19"),
  bcns = c("2010-02-24", "2011-12-17", "2012-08-18", "2017-06-29"),
  pgc5 = c("2009-07-24", "2011-12-17", "2017-11-20"),
  ptaa = c(),
  p435 = c("2011-10-12")
))

```
```{r}
str(nearby_cp)
```

The attribution method involves three main tasks, each implemented in a function that can run independently. Understanding these steps will help you understand the `Attribution_CP` function :

### 1. Model Identification for All Time Series

The `NoiseModel_Id` function fits the noise model for six difference series in `dataset` according to three ARMA models (AR(1), MA(1) and ARMA(1,1)) and the white noise.

```{r Model identification}
NoiseModel_Id(dataset = dataset,
              main_cp = main_cp,
              nearby_cp = nearby_cp)
```

### 2. Test of the Jump at Change-Points in Every Series

The `Test_CP` function tests the jump at change-points in a series by calling its name in the elementary dataframe in `dataset` : 

```{r Significance test}
TestRes = Test_CP(dataset$bces, Name_series = "GE", CP = main_cp[1],
                 noise_model = c("ARMA(1,1)"), limit = 100)
str(TestRes)

```

### 3. Configuration Prediction

`Prediction_CP` function runs predictions by applying predictors trained from our large dataset.  The input is a dataframe of t-values for six series.

```{r Prediction}
# arbitrary example: assume that the t values of 6 tests are : 2,2,2,1,1,1
ArbPred = Prediction_CP(data.frame(GE = 2, GGp = 2, GEp = 2, EEp = 1, GpEp = 1, GpE = 1))
ArbPred

```
### Examples
Finally, the main function includes all these steps. The noise model can be either fixed for all time series :

```{r Attribution1, echo=FALSE}

test_Attribution_CP_c1 = Attribution_CP(dataset,
                        main_cp = main_cp,
                        nearby_cp = nearby_cp,
                        noise_model_fix = "AR(1)")
str(test_Attribution_CP_c1)
```

Or specified for all time series:

```{r Attribution2}
NoiseModel = data.frame(
  bces = c("ARMA(1,1)", "AR(1)", "ARMA(1,1)", "AR(1)", "ARMA(1,1)", "ARMA(1,1)"),
  ptaa = c("ARMA(1,1)", rep("AR(1)", 3), "ARMA(1,1)", "AR(1)"),
  row.names = c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE")
)
Attribution_CP( dataset = dataset[c("bces", "ptaa")],
                main_cp = main_cp,
                nearby_cp = nearby_cp[c("bces", "ptaa")],
                noise_model_fix = NoiseModel)
```

Or set to `NULL`, and the noise models of all time series will be estimated using the `NoiseModel_Id` function :

```{r Attribution3}
Attribution_CP( dataset,
                main_cp = main_cp,
                nearby_cp = nearby_cp,
                noise_model_fix = NULL)
```

