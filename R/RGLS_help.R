

### creating sample data

library(lavaan)

mymodel<-'
f1 =~x1 + x2 + x3 + x4 + x5
f2 =~x6 + x7 + x8 + x9 + x10
f3 =~x11 + x12 + x13 + x14 + x15
f3 ~~ f1 + f2
f1 ~~ f2'



pop.model<-'

f1 =~ 0.7*x1 + 0.7*x2 + 0.75*x3 + 0.8*x4 + 0.8*x5
f2 =~ 0.7*x6 + 0.7*x7 + 0.75*x8 + 0.8*x9 + 0.8*x10
f3 =~ 0.7*x11 + 0.7*x12 + 0.75*x13 + 0.8*x14 + 0.8*x15
f3 ~~  0.4*f1 + 0.5*f2
f1 ~~ 0.3*f2

x1 ~~ 0.51*x1
x2 ~~ 0.51*x2
x3 ~~ 0.4375*x3
x4 ~~ 0.36*x4
x5 ~~ 0.36*x5
x6 ~~ 0.51*x6
x7 ~~ 0.51*x7
x8 ~~ 0.4375*x8
x9 ~~ 0.36*x9
x10 ~~ 0.36*x10
x11 ~~ 0.51*x11
x12 ~~ 0.51*x12
x13 ~~ 0.4375*x13
x14 ~~ 0.36*x14
x15 ~~ 0.36*x15
'

rgls_data<-simulateData(pop.model, model.type = "sem", meanstructure = FALSE, sample.nobs=1000L)


#' Regularized Generalized Least Squares (RGLS)
#'
#' This function conducts a regularized generalized least squares test for a structural equation model (SEM).
#'
#' @param model An object specifying the SEM model.
#' @param data A data frame containing the observed variables.
#' @param meanstructure Logical. Indicate whether to estimate the means.
#' @param int.ov.free Logical. Indicate whether to free intercepts for latent variables.
#' @param std.lv Logical. Indicate whether to standardize latent variables.
#' @return A list with test statistics, including the regularized likelihood ratio (T_RGLS) and p-value.
#' @examples
#' \dontrun{
#'   result <- rgls(model = mymodel, data = mydata)
#'   print(result)
#' }
#'
#' @keywords structural equation model SEM GLS regularized test
#' @import lavaan
#' @export
rgls <- function(model, data, meanstructure = FALSE, int.ov.free = FALSE, std.lv = FALSE) {
  fit_gls <- cfa(model = model, data = data, estimator = "GLS")

  S_N <- var(data)
  P <- ncol(data)
  N <- nrow(data)
  S <- cov(data)
  ev <- eigen(S)
  di <- ev$values
  Q <- ev$vectors
  tQ <- t(Q)
  alpha <- (1 + (median(ev$values))^2)^(-1)
  lam <- 10
  Ei <- (-N + sqrt(N^2 + 4 * lam * alpha * (N * di + lam * (1 - alpha)))) / (2 * lam * alpha)
  new_Ei <- diag(x = Ei, nrow = P, ncol = P)
  Sigma_hat <- Q %*% new_Ei %*% tQ

  Sigma_gls0 <- fitted(fit_gls)$cov
  Sigma_gls <- matrix(unlist(Sigma_gls0), ncol = P, byrow = TRUE)
  R <- (S_N - Sigma_gls) %*% solve(Sigma_hat)
  T_RGLS <- (N - 1) / 2 * sum(diag(R %*% R))
  df_RGLS <- lavInspect(fit_gls, "test")[[1]]$df

  # Calculate p-value
  p_value <- 1 - pchisq(T_RGLS, df_RGLS)

  # Return a list with results
  result <- list(T_RGLS = T_RGLS, p_value = p_value)
  return(result)
}



