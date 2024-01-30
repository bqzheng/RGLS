


### RGLS Function ###

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


