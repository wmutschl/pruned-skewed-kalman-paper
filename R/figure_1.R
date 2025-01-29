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

# set parameters
mu <- 0
nu <- 0
sigma <- matrix(1, 1, 1)
delta <- matrix(1, 1, 1)

# define grid
x <- seq(-4, 4, length = 200)

# compute densities
f1 <- dcsn(matrix(x, 200, 1), mu, sigma, gamma = matrix(-5, 1, 1), nu, delta)
f2 <- dcsn(matrix(x, 200, 1), mu, sigma, gamma = matrix(0, 1, 1), nu, delta)
f3 <- dcsn(matrix(x, 200, 1), mu, sigma, gamma = matrix(5, 1, 1), nu, delta)
f4 <- dcsn(matrix(x, 200, 1), mu, sigma, gamma = matrix(5, 1, 1), -8, delta)

# create plots
text_size <- 12
title_size <- 13

p_a <- ggplot(data.frame(x = x, f = f1), aes(x, f)) +
  geom_line(linewidth = 1.2) +
  ggtitle(TeX("(a) \\Gamma = -5, \\nu = 0")) +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  ) +
  ylab("density")
print(p_a)

p_b <- ggplot(data.frame(x = x, f = f2), aes(x, f)) +
  geom_line(linewidth = 1.2) +
  ggtitle(TeX("(b) \\Gamma = 0, \\nu = 0")) +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  ) +
  ylab("density")
print(p_b)

p_c <- ggplot(data.frame(x = x, f = f3), aes(x, f)) +
  geom_line(linewidth = 1.2) +
  ggtitle(TeX("(c) \\Gamma = 5, \\nu = 0")) +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  ) +
  ylab("density")
print(p_c)

p_d <- ggplot(data.frame(x = x, f = f4), aes(x, f)) +
  geom_line(linewidth = 1.2) +
  ggtitle(TeX("(d) \\Gamma = 5, \\nu = -8")) +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  ) +
  ylab("density")
print(p_d)

fig_1_top <- (p_a + p_b) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_1_bottom <- (p_c + p_d) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_1 <- (fig_1_top / fig_1_bottom)
print(fig_1)

# save plot
ggsave(
  filename = "figure_1.eps",
  plot = fig_1,
  dpi = 300,
  units = "px",
  width = 3000,
  height = 2000
)
