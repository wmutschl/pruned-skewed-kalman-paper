#' Compute pruned parameters of CSN distribution
#'
#' This function computes the pruned parameters of a CSN distributed random variable
#' by removing correlations below a specified threshold.
#'
#' @param sigma Scale parameter of CSN distribution, [p x p] matrix
#' @param gamma First skewness parameter of CSN distribution, [q x p] matrix
#' @param nu Second skewness parameter of CSN distribution, [q x 1] vector
#' @param delta Third skewness parameter of CSN distribution, [q x 1] matrix
#' @param tol Threshold at which correlations are pruned, scalar between 0 and 1. Default is 0.1
#'
#' @return A list containing the pruned parameters:
#' \itemize{
#'   \item Sigma - Pruned scale parameter
#'   \item Gamma - Pruned first skewness parameter
#'   \item nu - Pruned second skewness parameter
#'   \item Delta - Pruned third skewness parameter
#' }
#'
#' @copyright 2022-2025 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
#'
#' @note This is free software: you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation, either version 3 of the License, or
#' (at your option) any later version.
#'
#' This is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#' GNU General Public License for more details.
#'
#' This file is part of the replication files for the paper "Pruned Skewed
#' Kalman Filter and Smoother: With Application to DSGE models" by
#' Gaygysyz Guljanov, Willi Mutschler, Mark Trede.
pruned_csn_params <- function(sigma, gamma, nu, delta, tol) {
  p <- dim(sigma)[1]
  q <- dim(delta)[1]
  p1 <- sigma
  p2 <- sigma %*% t(gamma)
  p3 <- gamma %*% sigma
  p4 <- delta + gamma %*% sigma %*% t(gamma)
  cov_p <- rbind(cbind(p1, p2), cbind(p3, p4))
  maxval <- apply(cov2cor(cov_p)[(p + 1):(p + q), 1:p, drop = FALSE], 1, max)
  cutvector <- (maxval < tol)
  if (all(cutvector)) {
    gamma <- matrix(0, 1, 1)
    nu <- 0
    delta <- matrix(1, 1, 1)
  } else if (any(cutvector)) {
    p3 <- p3[!cutvector, , drop = FALSE]
    p4 <- p4[!cutvector, !cutvector, drop = FALSE]
    gamma <- p3 %*% solve(sigma)
    nu <- nu[!cutvector]
    delta <- p4 - p3 %*% t(gamma)
  }
  return(list(sigma = sigma,
              gamma = gamma,
              nu = nu,
              delta = delta))
}