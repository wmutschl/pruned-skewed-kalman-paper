#' Data visualization for Ireland (2004) dataset
#'
#' This script creates a figure showing the demeaned time series data
#' for output growth, inflation, and interest rate from 1980 to 2003,
#' along with summary statistics (standard deviation and skewness).
#'
#' @copyright 2026 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
library(tidyr)
library(e1071) # for skewness
library(gridExtra)
library(grid)

# load data
ireland_data <- read.table("../data/ireland2004_gpr.dat", header = FALSE)
colnames(ireland_data) <- c("Output growth", "Inflation", "Interest rate")

# create timeline (quarterly data from 1948 Q2 to 2003 Q4)
timeline <- seq(from = 1948.25, to = 2003, by = 0.25)

# add timeline to data
ireland_data$Year <- timeline

# filter for POST_1980 sample (as in Ireland 2004)
post_1980 <- ireland_data[ireland_data$Year >= 1980, ]

# demean the data
demean <- function(x) x - mean(x)
post_1980$`Output growth` <- demean(post_1980$`Output growth`)
post_1980$Inflation <- demean(post_1980$Inflation)
post_1980$`Interest rate` <- demean(post_1980$`Interest rate`)

# convert to percentage
post_1980$`Output growth` <- post_1980$`Output growth` * 100
post_1980$Inflation <- post_1980$Inflation * 100
post_1980$`Interest rate` <- post_1980$`Interest rate` * 100

# compute summary statistics
# use type = 2 for skewness (SAS/SPSS style, with bias correction)
stats <- data.frame(
  Series = c("Output growth", "Inflation", "Interest rate"),
  `Std-dev` = c(
    round(sd(post_1980$`Output growth`), 3),
    round(sd(post_1980$Inflation), 3),
    round(sd(post_1980$`Interest rate`), 3)
  ),
  Skew = c(
    round(skewness(post_1980$`Output growth`, type = 2), 3),
    round(skewness(post_1980$Inflation, type = 2), 3),
    round(skewness(post_1980$`Interest rate`, type = 2), 3)
  ),
  check.names = FALSE
)
print(stats)

# reshape data for ggplot
plot_data <- pivot_longer(
  post_1980,
  cols = c("Output growth", "Inflation", "Interest rate"),
  names_to = "Variable",
  values_to = "value"
)

# set factor levels for consistent ordering and legend
plot_data$Variable <- factor(
  plot_data$Variable,
  levels = c("Output growth", "Inflation", "Interest rate")
)

# plot settings
text_size <- 12
title_size <- 12
axis_text_size <- 9

# create statistics table as a grob with border (left aligned)
stats_table <- tableGrob(
  stats,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 7,
    core = list(
      fg_params = list(fontsize = 7, hjust = 0, x = 0.1),
      bg_params = list(fill = "white", col = NA),
      padding = unit(c(2, 2), "mm")
    ),
    colhead = list(
      fg_params = list(fontsize = 7, fontface = "bold", hjust = 0, x = 0.1),
      bg_params = list(fill = "grey95", col = NA),
      padding = unit(c(2, 2), "mm")
    )
  )
)

# add border around the table
stats_table <- gtable::gtable_add_grob(
  stats_table,
  grobs = grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey70", lwd = 1)),
  t = 1, b = nrow(stats_table), l = 1, r = ncol(stats_table)
)

# create main time series plot with annotation box
fig_data <- ggplot(plot_data, aes(x = Year, y = value, color = Variable, shape = Variable)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_color_manual(
    values = c("Output growth" = "darkgoldenrod", "Inflation" = "steelblue", "Interest rate" = "forestgreen")
  ) +
  scale_shape_manual(
    values = c("Output growth" = 16, "Inflation" = 17, "Interest rate" = 15)  # circle, triangle, square
  ) +
  scale_x_continuous(
    breaks = seq(1980, 2000, by = 5),
    limits = c(1980, 2003)
  ) +
  labs(
    x = "Year",
    y = "Demeaned values (in %)"
  ) +
  # add statistics table as annotation in upper right corner
  annotation_custom(
    grob = stats_table,
    xmin = 1995, xmax = 2003,
    ymin = 1.3, ymax = 2.3
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # Clean white background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    # Panel border
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.5),
    # Grid styling
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    # Typography
    text = element_text(size = text_size),
    axis.title.y = element_text(size = axis_text_size, margin = margin(r = 10)),
    axis.title.x = element_text(size = axis_text_size, margin = margin(t = 10)),
    axis.text = element_text(size = axis_text_size, color = "grey30"),
    # Legend inside panel at bottom
    legend.position = c(0.5, 0.05),
    legend.justification = c(0.5, 0),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    legend.title = element_blank(),
    legend.text = element_text(size = text_size - 4),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    legend.margin = margin(3, 5, 3, 5),
    # Plot margins
    plot.margin = margin(10, 15, 10, 10)
  )
print(fig_data)

# save figure
ggsave(
  filename = "figure_data.pdf",
  plot = fig_data,
  dpi = 300,
  units = "px",
  width = 2400,
  height = 1100
)
