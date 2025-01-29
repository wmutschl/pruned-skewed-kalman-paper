#' Density functions of bivariate CSN distributions with different skewness parameters
#'
#' This script replicates Figure 2 of the paper, showing density plots for bivariate CSN
#' distributions with varying skewness parameters.
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
#' Gaygysyz Guljanov, Willi Mutschler, and Mark Trede.

# housekeeping
setwd(this.path::here()) # set working directory
rm(list = ls()) # clear all variables and unload libraries

# load libraries
library(csn)
library(ggplot2)
library(latex2exp)
library(patchwork)

# define plotting function
pltobj <- function(xx, mu, sigma, gamma, nu, delta, ttl = "") {
  text_size <- 12
  title_size <- 13
  z <- dcsn(xx, mu, sigma, gamma, nu, delta)
  p <- ggplot(data.frame(x1 = xx[, 1], x2 = xx [, 2], f = z), aes(x1, x2, fill = z)) +
    geom_raster(interpolate = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = TeX("$x_1$"), y = TeX("$x_2$")) +
    guides(fill = "none") +
    scale_fill_gradient(low = "white", high = "black") +
    ggtitle(ttl) +
    geom_vline(xintercept = 0, colour = "black") +
    geom_hline(yintercept = 0, colour = "black") +
    theme(
      text = element_text(size = text_size),
      axis.title.y = element_text(size = text_size),
      axis.title.x = element_text(size = text_size),
      plot.title = element_text(size = title_size, hjust = 0.5)
    )
  return(p)
}

# define grid
x1 <- seq(-1.2, 2, length = 51)
x2 <- seq(-1.2, 2, length = 53)
xx <- cbind(rep(x1, times = length(x2)), rep(x2, each = length(x1)))

# define parameters
mu_a <- c(0, 0)
mu_b <- c(0, 0)
mu_c <- c(0, 0)
mu_d <- c(0, 0)

sigma_a <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
sigma_b <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
sigma_c <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
sigma_d <- matrix(c(1, 0.7, 0.7, 1), 2, 2)

gamma_a <- matrix(c(6, 0, 0, 6), 2, 2)
gamma_b <- matrix(c(6, 0, 0, -6), 2, 2)
gamma_c <- matrix(c(6, 6, 6, 6), 2, 2)
gamma_d <- matrix(c(6, 0, 0, 6), 2, 2)

nu_a <- c(0, 0)
nu_b <- c(0, 0)
nu_c <- c(0, 0)
nu_d <- c(-6, -6)

delta_a <- diag(2)
delta_b <- diag(2)
delta_c <- diag(2)
delta_d <- diag(2)

# create plots
p_a <- pltobj(xx, mu = mu_a, sigma = sigma_a, gamma = gamma_a, nu = nu_a, delta = delta_a, "(a)")
p_b <- pltobj(xx, mu = mu_b, sigma = sigma_b, gamma = gamma_b, nu = nu_b, delta = delta_b, "(b)")
p_c <- pltobj(xx, mu = mu_c, sigma = sigma_c, gamma = gamma_c, nu = nu_c, delta = delta_c, "(c)")
p_d <- pltobj(xx, mu = mu_d, sigma = sigma_d, gamma = gamma_d, nu = nu_d, delta = delta_d, "(d)")

fig_2_top <- (p_a + p_b) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_2_bottom <- (p_c + p_d) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_2 <- (fig_2_top / fig_2_bottom)
print(fig_2)

# save plot
ggsave(
  filename = "figure_2.eps",
  plot = fig_2,
  dpi = 300,
  units = "px",
  width = 3000,
  height = 2000
)