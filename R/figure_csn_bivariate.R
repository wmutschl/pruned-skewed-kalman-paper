#' Density functions of bivariate CSN distributions with different skewness parameters
#'
#' This script replicates Figure 2 of the paper, showing density plots for bivariate CSN
#' distributions with varying skewness parameters.
#'
#' @copyright 2022-2026 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
pltobj <- function(xx, mu, sigma, gamma, nu, delta, ttl = "", x_label = "", y_label = "") {
  text_size <- 14
  title_size <- 15
  axis_text_size <- 12
  z <- dcsn(xx, mu, sigma, gamma, nu, delta)

  p <- ggplot(data.frame(x1 = xx[, 1], x2 = xx [, 2], f = z), aes(x1, x2, fill = z)) +
    geom_raster(interpolate = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = x_label,
      y = y_label,
      title = ttl
    ) +
    guides(fill = "none") +
    scale_fill_gradient(low = "white", high = "black") +
    geom_vline(xintercept = 0, colour = "black") +
    geom_hline(yintercept = 0, colour = "black") +
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

# define grid
x1 <- seq(-1.2, 2, length = 51)
x2 <- seq(-1.2, 2, length = 53)
xx <- cbind(rep(x1, times = length(x2)), rep(x2, each = length(x1)))

# define parameters (mu and sigma are the same for all)
mu <- c(0, 0)
sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)

# (a) Γ = diag(6), ν = 0, Δ = I
gamma_a <- matrix(c(6, 0, 0, 6), 2, 2)
nu_a <- c(0, 0)
delta_a <- diag(2)

# (b) Γ = diag(6), ν = 0, Δ = [[1, 0.99], [0.99, 1]]
gamma_b <- matrix(c(6, 0, 0, 6), 2, 2)
nu_b <- c(0, 0)
delta_b <- matrix(c(1, 0.99, 0.99, 1), 2, 2)

# (c) Γ = diag(6), ν = [-3, -3], Δ = I
gamma_c <- matrix(c(6, 0, 0, 6), 2, 2)
nu_c <- c(-3, -3)
delta_c <- diag(2)

# (d) Γ = diag(6, -6), ν = 0, Δ = I
gamma_d <- matrix(c(6, 0, 0, -6), 2, 2)
nu_d <- c(0, 0)
delta_d <- diag(2)

# (e) Γ = [[6, 6], [6, 6]], ν = 0, Δ = I
gamma_e <- matrix(c(6, 6, 6, 6), 2, 2)
nu_e <- c(0, 0)
delta_e <- diag(2)

# (f) Γ = diag(12, 18), ν = [-6, -9], Δ = diag(4, 9)
gamma_f <- matrix(c(12, 0, 0, 18), 2, 2)
nu_f <- c(-6, -9)
delta_f <- matrix(c(4, 0, 0, 9), 2, 2)

# create plots

p_a <- pltobj(xx, mu = mu, sigma = sigma, gamma = gamma_a, nu = nu_a, delta = delta_a,
              ttl = TeX("(a) $\\Gamma = diag(6,6)$, $\\nu = (0,0)'$, $\\Delta = I_{2}$"),
              x_label = "", y_label = TeX("$x_2$"))

p_b <- pltobj(xx, mu = mu, sigma = sigma, gamma = gamma_b, nu = nu_b, delta = delta_b,
              ttl = TeX("(b) $\\Gamma = diag(6,6)$, $\\nu = (0,0)'$, $vech(\\Delta) = (1,0.99,1)'$"),
              x_label = "", y_label = "")

p_c <- pltobj(xx, mu = mu, sigma = sigma, gamma = gamma_c, nu = nu_c, delta = delta_c,
              ttl = TeX("(c) $\\Gamma = diag(6,6)$, $\\nu = -(3,3)'$, $\\Delta = I_{2}$"),
              x_label = "", y_label = "")

p_d <- pltobj(xx, mu = mu, sigma = sigma, gamma = gamma_d, nu = nu_d, delta = delta_d,
              ttl = TeX("(d) $\\Gamma = diag(6,-6)$, $\\nu = (0,0)'$, $\\Delta = I_{2}$"),
              x_label = TeX("$x_1$"), y_label = TeX("$x_2$"))

p_e <- pltobj(xx, mu = mu, sigma = sigma, gamma = gamma_e, nu = nu_e, delta = delta_e,
              ttl = TeX("(e) $\\Gamma = diag(6,6)$, $\\nu = (0,0)'$, $\\Delta = I_{2}$"),
              x_label = TeX("$x_1$"), y_label = "")

p_f <- pltobj(xx, mu = mu, sigma = sigma, gamma = gamma_f, nu = nu_f, delta = delta_f,
              ttl = TeX("(f) $\\Gamma = diag(12,18)$, $\\nu = -(6,9)'$, $\\Delta = diag(4,9)$"),
              x_label = TeX("$x_1$"), y_label = "")

fig_2_top <- (p_a + p_b + p_c)
fig_2_bottom <- (p_d + p_e + p_f)
fig_2 <- (fig_2_top / fig_2_bottom) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.spacing = unit(0.5, "cm"),
    legend.margin = margin(t = 10, b = 5),
    plot.margin = margin(10, 10, 10, 10)
  )
print(fig_2)

# save plot
ggsave(
  filename = "figure_csn_bivariate.pdf",
  plot = fig_2,
  dpi = 300,
  units = "px",
  width = 4500,
  height = 3000
)
