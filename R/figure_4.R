#' Estimated probability density functions of shocks
#'
#' This script replicates Figure 4 of the paper, comparing the
#' estimated probability density functions of shocks based on the
#' maximum likelihood estimates for the Gaussian and CSN models.
#'
#' @copyright 2025 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
library(latex2exp)
library(R.matlab)

# load Gaussian estimates
params_gauss <- readMat("../ireland2004_ml_1_gaussian/Output/ireland2004_ml_1_gaussian_shock_params.mat")
sigma_gauss <- params_gauss$csn[2][[1]]

# load CSN estimates
params_csn <- readMat("../ireland2004_ml_3_csn/Output/ireland2004_ml_3_csn_shock_params.mat")
mu_csn <- params_csn$csn[1][[1]]
sigma_csn <- params_csn$csn[2][[1]]
gamma_csn <- params_csn$csn[3][[1]]
nu_csn <- params_csn$csn[4][[1]]
delta_csn <- params_csn$csn[5][[1]]

# define grid for plotting
x_pref <- seq(-0.1 * 100, 0.1 * 100, length = 200)
x_cost <- seq(-0.001 * 100, 0.001 * 100, length = 200)
x_prod <- seq(-0.03 * 100, 0.03 * 100, length = 200)
x_monp <- seq(-0.01 * 100, 0.01 * 100, length = 200)

# compute densities
densities <- data.frame(
  f_pref_gauss = dnorm(x_pref, 0, sqrt(sigma_gauss[1, 1])),
  f_cost_gauss = dnorm(x_cost, 0, sqrt(sigma_gauss[2, 2])),
  f_prod_gauss = dnorm(x_prod, 0, sqrt(sigma_gauss[3, 3])),
  f_monp_gauss = dnorm(x_monp, 0, sqrt(sigma_gauss[4, 4])),
  f_pref_csn = dcsn(x_pref, mu_csn[1, 1], sigma_csn[1, 1], gamma_csn[1, 1], nu_csn[1, 1], delta_csn[1, 1]),
  f_cost_csn = dcsn(x_cost, mu_csn[2, 1], sigma_csn[2, 2], gamma_csn[2, 2], nu_csn[2, 1], delta_csn[2, 2]),
  f_prod_csn = dcsn(x_prod, mu_csn[3, 1], sigma_csn[3, 3], gamma_csn[3, 3], nu_csn[3, 1], delta_csn[3, 3]),
  f_monp_csn = dcsn(x_monp, mu_csn[4, 1], sigma_csn[4, 4], gamma_csn[4, 4], nu_csn[4, 1], delta_csn[4, 4])
)

# plot settings
text_size <- 12
title_size <- 13

# preference shock plot
p_pref <- ggplot(densities, aes(x_pref, f_pref_csn)) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = f_pref_gauss), linetype = "dashed", linewidth = 1.2) +
  ylab("") +
  xlab(TeX(r"($\eta_{a}$)")) +
  ggtitle("Preference shock") +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  )
print(p_pref)

# cost-push shock plot
p_cost <- ggplot(densities, aes(x_cost, f_cost_csn)) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = f_cost_gauss), linetype = "dashed", linewidth = 1.2) +
  ylab("") +
  xlab(TeX(r"($\eta_{e}$)")) +
  ggtitle("Cost-push shock") +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  )
print(p_cost)

# productivity shock plot
p_prod <- ggplot(densities, aes(x_prod, f_prod_csn)) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = f_prod_gauss), linetype = "dashed", linewidth = 1.2) +
  ylab("") +
  xlab(TeX(r"($\eta_{z}$)")) +
  ggtitle("Productivity shock") +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  )
print(p_prod)

# monetary policy shock plot
p_monp <- ggplot(densities, aes(x_monp, f_monp_csn)) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(y = f_monp_gauss), linetype = "dashed", linewidth = 1.2) +
  ylab("") +
  xlab(TeX(r"($\eta_{r}$)")) +
  ggtitle("Monetary policy shock") +
  theme(
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    plot.title = element_text(size = title_size, hjust = 0.5)
  )
print(p_monp)

# create figure
fig_4_top <- (p_pref + p_cost) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_4_bottom <- (p_prod + p_monp) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_4 <- (fig_4_top / fig_4_bottom)
print(fig_4)

# save figure
ggsave(
  filename = "figure_4.eps",
  plot = fig_4,
  dpi = 300,
  units = "px",
  width = 3000,
  height = 2000
)
