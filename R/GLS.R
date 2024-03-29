#' Fit GLS (known variance-covariance matrix)
#'
#' @param phi  A numeric value indicating the parameter of the AR component.
#'
#' @param theta A numeric value indicating the parameter of the MA component.
#'
#' @param var.t A variance matrix computed by the moving window.
#'
#' @param df_XY A data frame resulted from construct_XY_df that includes a column
#' of the signal (Y) and columns of X: 8 Fourier components (if selected), a global
#' mean column, and a jump column (if break_ind is not NULL).
#'
#' @return A list containing the estimates and their standard errors, t-values,
#'  and p-values.
#'
#' @importFrom stats ARMAacf na.omit pnorm toeplitz
#'
#' @export
#'
#' @keywords internal
GLS <- function(phi, theta, var.t, df_XY){

  signal <- NULL

  # cor_matrix fc to produce the variance corvariance matrix
  cor_matrix <- function(phi, theta, n.full, n.true){

    if(phi ==0 & theta ==0){
      y = rep(1, n.full)
    }else{
      y = ARMAacf(ar = phi, ma = theta, lag.max = n.true-1, pacf = FALSE)
    }

    m = toeplitz(y)

    return(m)
  }

  # compute the variance covariance matrix
  ind1 = which(is.na(df_XY$signal)==FALSE)
  var.t.na = var.t[ind1]
  var.matrix = diag(sqrt(var.t.na))
  cor.matrix = cor_matrix(phi, theta, n.full = nrow(df_XY), n.true = nrow(na.omit(df_XY)))

  if(phi==0 & theta==0){
    cov.var= diag(var.t.na)
  }else{
    cov.var0 = var.matrix  %*%  cor.matrix
    cov.var = cov.var0 %*%  var.matrix
  }

  # estimate
  X = as.matrix(df_XY%>% dplyr::select(-signal))
  X = X[ind1,]
  term2 = solve(cov.var)
  term1 = t(X) %*% term2 %*% X
  term3 = solve(term1)
  beta = term3 %*% t(X) %*% term2  %*% (as.matrix(df_XY$signal[ind1]))
  var.beta = term3

  # form the frame of result as in the ols
  residual = df_XY$signal[ind1] - X%*%beta
  t.val = beta/sqrt((diag(var.beta)))
  p.val = round(pnorm(-abs(t.val), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 4)
  fitted.gls = X%*%beta

  fit.gls = data.frame(Estimate = format(beta, scientific = TRUE, digits = 3),
                       `Std. Error` = format(sqrt(diag(var.beta)), scientific = TRUE, digits = 3),
                       `t value` = format(t.val, scientific = TRUE, digits = 3),
                       `Pr(>|t|)` = format(p.val, scientific = TRUE, digits = 3),
                       check.names = FALSE
  )

  Res <- list(Coefficients = beta,
              t.table = fit.gls,
              vcov = var.beta,
              residual = residual,
              fitted = fitted.gls)

  return(Res)

}
