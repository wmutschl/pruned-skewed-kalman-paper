#' Illustration of pruned distribution
#'
#' This script replicates Figure 3 of the paper, comparing the pdf and cdf
#' of an univariate CSN distributed variable with its pruned approximation.
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
library(ggplot2)
library(csn)
library(patchwork)

# load pruned_csn_params function
source("csn_prune_params.R")

# define parameters
sigma <- matrix(1, 1, 1)
gamma <- matrix(c(6, 0.1), 2, 1)
nu <- matrix(c(0, 0), 2, 1)
delta <- matrix(c(1, -0.1, -0.1, 1), 2, 2)

# compute unpruned parameters
cov_p <- rbind(cbind(sigma, sigma %*% t(gamma)), cbind(gamma %*% sigma, delta + gamma %*% sigma %*% t(gamma)))
corr_r <- cov2cor(cov_p)
print(corr_r)

# compute pruned parameters
pruned_params <- pruned_csn_params(sigma, gamma, nu, delta, 0.1)
sigma_pruned <- pruned_params$sigma
gamma_pruned <- pruned_params$gamma
nu_pruned <- pruned_params$nu
delta_pruned <- pruned_params$delta

# compute densities
x <- seq(-1, 4, length = 200)
densities <- data.frame(
  x = x,
  pdf_unpruned = dcsn(x, 0, sigma,        gamma,        as.vector(nu),        delta),
  pdf_pruned   = dcsn(x, 0, sigma_pruned, gamma_pruned, as.vector(nu_pruned), delta_pruned),
  cdf_unpruned = pcsn(x, 0, sigma,        gamma,        as.vector(nu),        delta),
  cdf_pruned   = pcsn(x, 0, sigma_pruned, gamma_pruned, as.vector(nu_pruned), delta_pruned)
)

# create plots
text_size <- 12
title_size <- 13

p_pdfs <- ggplot(densities, aes(x, pdf_unpruned)) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = pdf_pruned), linetype = "dashed", linewidth = 1.2) +
  ggtitle("Probability density function") +
  ylab("") +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  )
print(p_pdfs)

p_cdfs <- ggplot(densities, aes(x, cdf_unpruned)) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = cdf_pruned), linetype = "dashed", linewidth = 1.2) +
  ggtitle("Cumulative distribution function") +
  ylab("") +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  )
print(p_cdfs)

fig_3 <- (p_pdfs + p_cdfs) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
print(fig_3)

# save plot
ggsave(
  filename = "figure_3.eps",
  plot = fig_3,
  dpi = 300,
  units = "px",
  width = 3000,
  height = 1000
)
