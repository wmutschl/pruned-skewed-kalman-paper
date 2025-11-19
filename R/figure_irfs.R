#' Impulse response functions to a monetary policy shock
#'
#' This script replicates Figure 5 of the paper, comparing the
#' impulse response functions to a monetary policy shock estimated
#' using the Gaussian and CSN model variants.
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
library(dplyr)
library(tidyr)

irfs <- readMat("../ireland2004_irfs/Output/ireland2004_irfs.mat")
gaussian_neg <- irfs$irfs.gaussian.neg[, , 1]
gaussian_pos <- irfs$irfs.gaussian.pos[, , 1]
csn_neg <- irfs$irfs.csn.neg[, , 1]
csn_pos <- irfs$irfs.csn.pos[, , 1]

# get impulse response functions from simulations done in Dynare
irf_length <- 15

shock_suffix <- c(
  "Preference" = "a",
  "Cost-Push" = "e",
  "Productivity" = "z",
  "Monetary Policy" = "r"
)

variable_prefix <- c(
  "Annualized Interest Rate" = "rAnnualized.eta",
  "Output Gap" = "xhat.eta",
  "Output Growth" = "ghat.eta",
  "Annualized Inflation" = "piAnnualized.eta"
)

models <- c("Gaussian" = "gaussian", "CSN" = "csn")
signs  <- c("16th" = "neg", "84th" = "pos")

# Create the "skeleton" tibble:
data_irfs <- expand.grid(
  distribution = names(models),
  sign         = names(signs),
  shock        = names(shock_suffix),
  variable     = names(variable_prefix),
  time         = seq_len(irf_length),
  stringsAsFactors = TRUE
)

# Fill in the `value` column using rowwise + get().
# This approach constructs the object name (e.g. "gaussian_neg")
# and the column name (e.g. "rhat.eta.a") for each row.
data_irfs <- data_irfs %>%
  rowwise() %>%
  mutate(
    dataset_name = paste0(models[distribution], "_", signs[sign]),
    col_name     = paste0(variable_prefix[variable], ".", shock_suffix[shock]),
    value        = get(dataset_name)[[col_name]][time]
  ) %>%
  ungroup() %>%
  select(distribution, sign, shock, variable, time, value)

# check if the value is correct (also check in MATLAB irfs_csn_pos.xhat_eta_r)
data_irfs %>%
  filter(
    distribution == "CSN",
    sign         == "84th",
    shock        == "Monetary Policy",
    variable     == "Output Gap"
  ) %>%
  pull(value) == csn_pos$xhat.eta.r

# plot options
text_size <- 14
axis_text_size <- 12
strip_size <- 13
legend_size <- 14

# IRFs with respect to monetary policy shock
p_irfs <- ggplot(
  data_irfs %>%
    filter(shock == "Monetary Policy",
           variable %in% c("Output Gap", "Annualized Inflation", "Annualized Interest Rate")) %>%
    mutate(value = value * 100, time = time - 1) %>%
    mutate(variable = factor(variable,
                             levels = c("Annualized Interest Rate", "Annualized Inflation", "Output Gap"))),
  aes(x = time, y = value, color = sign, linetype = distribution)
) +
  geom_line(size = 1.2) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  labs(
    x = "Time horizon (in quarters)",
    y = "Percentage deviation from steady-state (in %)",
    color = "Shock size (percentile)",
    linetype = "Shock distribution"
  ) +
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
    strip.text = element_text(size = strip_size),
    strip.background = element_blank(),
    # Legend settings
    legend.title = element_text(size = legend_size),
    legend.text = element_text(size = text_size),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.spacing = unit(0.5, "cm"),
    legend.margin = margin(t = 10, b = 5),
    # Plot margins
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_color_manual(values = c("16th" = "grey50", "84th" = "grey20")) +
  scale_linetype_manual(values = c("Gaussian" = "dotdash", "CSN" = "solid")) +
  scale_x_continuous(limits = c(0, irf_length - 1), breaks = seq(0, irf_length - 1, by = 2)) +
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, by = 0.1))
print(p_irfs)

ggsave(
  filename = "figure_irfs.pdf",
  plot = p_irfs,
  dpi = 300,
  units = "px",
  width = 3000,
  height = 1500
)
