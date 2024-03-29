#' Compute the sliding window standard deviation by ScaleTau estimator
#'
#' @param Y A data frame contains at least 2 columns named Date and
#' name.var. The Date column must be in Date type on this format: "%Y-%m-%d"),
#' and the Name_series column is numeric (can be NA).
#'
#' @param name.var The name of the column in "Y" that we want to test,
#' as a character string.
#'
#' @param length.wind An integer value specifying the length of the window used
#' to compute the sliding standard deviation. Default is 60.
#'
#' @return A vector of values of the sliding standard deviation computed by the
#' scaleTau estimator.
#'
#' @export
#'
#' @importFrom robustbase scaleTau2
#'
#' @importFrom tidyr complete
#'
#' @importFrom stats na.omit approx
#'
#' @keywords internal

Sliding_std <- function(Y, name.var, length.wind = 60){

  # require date in the dataset, return std at time t
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = nrow(Y1)
  sigma.est1 = rep(NA, n)

  for (i in c(1:n)) {
    begin = max(i-(length.wind-1),1)
    end = min(n, i+length.wind)
    x.i = x[begin:end]
    x.idiff = (x.i)
    thre = 30
    if(i < 30 | i > (n-30)){
      thre = 16
    }
    if(length(na.omit(x.idiff)) <= thre){
      sd <- NA
    }else{
      sd <- robustbase::scaleTau2(na.omit(x.idiff), consistency = "finiteSample")
    }
    sigma.est1[i] <- sd
  }
  # linear regression of the variance for gaps  MAYBE REPLACE BY INTERPOLATION FUNCTION
  s = sigma.est1
  if (sum(is.na(s)) != 0 & sum(is.na(s)) != length(s)){
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    if(na.ind[1] == 1){
      ind.stop = which(is.na(s)==FALSE)[1]-1
      na.ind <- na.ind[-c(1:ind.stop)]
    }else if (is.na(s[n]) == 1){
      m = which(is.na(s) == FALSE)
      ind.start = m[length(m)]
      na.ind <- na.ind[-which(na.ind %in% c(ind.start:n))]
    }
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  sigma.est = s[which(Y1$date %in% Y$date)]

  return(sigma.est)

}
