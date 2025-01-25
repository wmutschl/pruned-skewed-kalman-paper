# Replicates Figure 1: "Density functions of univariate CSN distributions with
# different skewness parameters"
# =========================================================================
# Copyright © 2022-2023 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# -------------------------------------------------------------------------
# This file is part of the replication files for the paper "Pruned Skewed
# Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
# With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
# by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
# =========================================================================
library(csn)
library(ggplot2)
library(gridExtra)
library(viridis)
library(latex2exp)
library(patchwork)
# set working directory
setwd(this.path::here())

source("R/csn_prune_params.R")

mu <- 0
nu <- 0
Sigma <- matrix(1,1,1)
Delta <- matrix(1,1,1)

x <- seq(-4,4,length=200)
f1 <- dcsn(matrix(x,200,1),mu,Sigma,gamma=matrix(-5,1,1),nu,Delta)
f2 <- dcsn(matrix(x,200,1),mu,Sigma,gamma=matrix(0,1,1),nu,Delta)
f3 <- dcsn(matrix(x,200,1),mu,Sigma,gamma=matrix(5,1,1),nu,Delta)
f4 <- dcsn(matrix(x,200,1),mu,Sigma,gamma=matrix(5,1,1),-8,Delta)

text_size = 12
title_size = 13

p1 <- ggplot(data.frame(x=x,f=f1),aes(x,f)) +
    geom_line(linewidth=1.2) +
    ggtitle(TeX("(a) \\Gamma = -5, \\nu = 0")) +
    theme(
      text = element_text(size=text_size),
      axis.title.y = element_text(size=text_size),
      axis.title.x = element_text(size=text_size),
      plot.title = element_text(size=title_size, hjust=0.5)
    )+
    ylab("density")
print(p1)
p2 <- ggplot(data.frame(x=x,f=f2),aes(x,f))+
    geom_line(linewidth=1.2)+
    ggtitle(TeX("(b) \\Gamma = 0, \\nu = 0"))+
    theme(
      text = element_text(size=text_size),
      axis.title.y = element_text(size=text_size),
      axis.title.x = element_text(size=text_size),
      plot.title = element_text(size=title_size, hjust=0.5)
    )+
    ylab("density")
p2
p3 <- ggplot(data.frame(x=x,f=f3),aes(x,f))+
    geom_line(linewidth=1.2)+
    ggtitle(TeX("(c) \\Gamma = 5, \\nu = 0"))+
    theme(
      text = element_text(size=text_size),
      axis.title.y = element_text(size=text_size),
      axis.title.x = element_text(size=text_size),
      plot.title = element_text(size=title_size, hjust=0.5)
    )+
    ylab("density")
p4 <- ggplot(data.frame(x=x,f=f4),aes(x,f))+
    geom_line(linewidth=1.2)+
    ggtitle(TeX("(d) \\Gamma = 5, \\nu = -8"))+
    theme(
      text = element_text(size=text_size),
      axis.title.y = element_text(size=text_size),
      axis.title.x = element_text(size=text_size),
      plot.title = element_text(size=title_size, hjust=0.5)
    )+
    ylab("density")

grid.arrange(p1, p2, p3, p4, ncol=2)

fig_1_top    <- (p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_1_bottom <- (p3 + p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_1 <- (fig_1_top / fig_1_bottom)
print(fig_1)
ggsave(filename = "plots/R/fig_1.pdf", plot = fig_1, dpi = 300, units = "px", width = 3000, height = 2000)
