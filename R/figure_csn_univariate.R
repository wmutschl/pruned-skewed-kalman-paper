#' Density functions of univariate CSN distributions with different skewness parameters
#'
#' This script replicates Figure 1 of the paper, showing density plots for CSN distributions
#' with varying skewness parameters.
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

# function to compute theoretical moments of univariate CSN
compute_moments <- function(mu, sigma, gamma, nu, delta) {
  S <- delta + gamma^2 * sigma
  kappa <- nu / sqrt(S)
  beta <- (gamma * sigma) / sqrt(S)
  r0 <- dnorm(-kappa) / pnorm(-kappa)
  kappa1 <- mu + beta * r0 # first cumulant
  kappa2 <- sigma + beta^2 * (kappa * r0 - r0^2) # second cumulant
  kappa3 <- beta^3 * (r0 * (kappa^2 - 1) - 3 * kappa * r0^2 + 2 * r0^3) # third cumulant
  return(list(mean = kappa1, var = kappa2, skewness = kappa3 / (kappa2^(3/2))))
}

# function to create CSN density plot
plot_csn_density <- function(mu, sigma, gamma, nu, delta,
                             subplot_label = "",
                             x_limits = c(-3.6, 3.6),
                             x_label = "x",
                             y_label = "",
                             text_size = 14,
                             title_size = 15,
                             axis_text_size = 12) {
  # Compute moments
  csn_moments <- compute_moments(mu, sigma, gamma, nu, delta)
  # Create x sequence and compute density
  x <- seq(x_limits[1], x_limits[2], length = 200)
  f <- dcsn(matrix(x, 200, 1), mu, sigma, gamma, nu, delta)
  # Create plot
  p <- ggplot(data.frame(x = x, f = f), aes(x, f)) +
    geom_line(linewidth = 1.2) +
    labs(
      title = TeX(sprintf("%s $\\mu = %.2f$, $\\Sigma = %.2f$, $\\Gamma = %.2f$, $\\nu = %.2f$, $\\Delta = %.2f$",
                          subplot_label, mu, sigma, gamma, nu, delta)),
      subtitle = TeX(sprintf("$E[X] = %.2f$, $Var[X] = %.2f^{2}$, $Skew[X] = %.2f$",
                             csn_moments$mean, sqrt(csn_moments$var), csn_moments$skewness)),
      x = x_label,
      y = y_label
    ) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    theme_minimal(base_size = 14) +
    theme(
      # Clean white background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      # Add subtle panel border for better definition
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.5),
      # Grid styling for better readability
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      # Typography settings
      text = element_text(size = text_size),
      axis.title.y = element_text(size = text_size, margin = margin(r = 10)),
      axis.title.x = element_text(size = text_size, margin = margin(t = 10)),
      axis.text = element_text(size = axis_text_size, color = "grey30"),
      plot.title = element_text(size = title_size, hjust = 0.5, face = "plain", margin = margin(b = 5)),
      plot.subtitle = element_text(size = text_size - 1, hjust = 0.5, margin = margin(b = 10)),
      # Plot margins
      plot.margin = margin(10, 10, 10, 10)
    )
  return(p)
}

################################
# Figure 1 (a): CSN(0,1,0,0,1) #
################################
p_a <- plot_csn_density(
  mu = 0, sigma = 1, gamma = 0, nu = 0, delta = 1,
  subplot_label = "(a)", x_label = "", y_label = "Density"
)
print(p_a)

################################
# Figure 1 (b): CSN(0,1,5,0,1) #
################################
p_b <- plot_csn_density(
  mu = 0, sigma = 1, gamma = 5, nu = 0, delta = 1,
  subplot_label = "(b)", x_label = "", y_label = ""
)
print(p_b)

####################################
# Figure 1 (c): CSN(mu0,2,10,20,1) #
####################################
# mu0 such that E[X] = 0
csn_moments <- compute_moments(mu = 0, sigma = 2, gamma = 10, nu = 20, delta = 1)
mu0 <- -csn_moments$mean
p_c <- plot_csn_density(
  mu = mu0, sigma = 2, gamma = 10, nu = 20, delta = 1,
  subplot_label = "(c)", x_label = "", y_label = ""
)
print(p_c)

####################################
# Figure 1 (d): CSN(0,1,5,-8,0.01) #
####################################
p_d <- plot_csn_density(
  mu = 0, sigma = 1, gamma = 5, nu = -8, delta = 0.01,
  subplot_label = "(d)", x_label = "x", y_label = "Density"
)
print(p_d)

#################################
# Figure 1 (e): CSN(0,1,-5,0,1) #
#################################
p_e <- plot_csn_density(
  mu = 0, sigma = 1, gamma = -5, nu = 0, delta = 1,
  subplot_label = "(e)", x_label = "x", y_label = ""
)
print(p_e)

########################################
# Figure 1 (f): CSN(mu0,5,-10,20,0.01) #
########################################
# mu0 such that E[X] = 0
csn_moments <- compute_moments(mu = 0, sigma = 5, gamma = -10, nu = 20, delta = 0.01)
mu0 <- -csn_moments$mean
p_f <- plot_csn_density(
  mu = mu0, sigma = 5, gamma = -10, nu = 20, delta = 0.01,
  subplot_label = "(f)", x_label = "x", y_label = ""
)
print(p_f)

############
# Figure 1 #
############
fig_1_top <- (p_a + p_b + p_c)
fig_1_bottom <- (p_d + p_e + p_f)
fig_1 <- (fig_1_top / fig_1_bottom) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.spacing = unit(0.5, "cm"),
    legend.margin = margin(t = 10, b = 5),
    plot.margin = margin(10, 10, 10, 10)
  )
print(fig_1)

# save plot
ggsave(
  filename = "figure_csn_univariate.pdf",
  plot = fig_1,
  dpi = 300,
  units = "px",
  width = 4000,
  height = 2000
)
