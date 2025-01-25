# Replicates Figure 3: "Probability density functions and cumulative distribution
# functions of a CSN distributed random variable with two skewness dimensions
# and the approximating pruned CSN distribution with one skewness dimenstions"
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
library(ggplot2)
library(gridExtra)
library(csn)
library(latex2exp)
library(R.matlab)
# set working directory
setwd(this.path::here())

source("R/csn_prune_params.R")

paramsGauss <- readMat("ireland2004_ml_1_gaussian/Output/ireland2004_ml_1_gaussian_shock_params.mat")
SigmaGauss <- paramsGauss$csn[2][[1]]

paramsCSN   <- readMat("ireland2004_ml_3_csn/Output/ireland2004_ml_3_csn_shock_params.mat")
muCSN <- paramsCSN$csn[1][[1]]
SigmaCSN <- paramsCSN$csn[2][[1]]
GammaCSN <- paramsCSN$csn[3][[1]]
nuCSN <- paramsCSN$csn[4][[1]]
DeltaCSN <- paramsCSN$csn[5][[1]]

x1 <- seq(-0.1*100,0.1*100,length=200)
x2 <- seq(-0.001*100,0.001*100,length=200)
x3 <- seq(-0.03*100,0.03*100,length=200)
x4 <- seq(-0.01*100,0.01*100,length=200)

H <- data.frame(f1_gauss=dnorm(x1,0,sqrt(SigmaGauss[1,1])),
                f2_gauss=dnorm(x2,0,sqrt(SigmaGauss[2,2])),
                f3_gauss=dnorm(x3,0,sqrt(SigmaGauss[3,3])),
                f4_gauss=dnorm(x4,0,sqrt(SigmaGauss[4,4])),
                f1=dcsn(x1,muCSN[1,1],SigmaCSN[1,1],GammaCSN[1,1],nuCSN[1,1],DeltaCSN[1,1]),
                f2=dcsn(x2,muCSN[2,1],SigmaCSN[2,2],GammaCSN[2,2],nuCSN[2,1],DeltaCSN[2,2]),
                f3=dcsn(x3,muCSN[3,1],SigmaCSN[3,3],GammaCSN[3,3],nuCSN[3,1],DeltaCSN[3,3]),
                f4=dcsn(x4,muCSN[4,1],SigmaCSN[4,4],GammaCSN[4,4],nuCSN[4,1],DeltaCSN[4,4])
)

text_size = 12
title_size = 13

p1 <- ggplot(H, aes(x1,f1)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f1_gauss),linetype="dashed",linewidth=1.2) +
  ylab("") +
  xlab(TeX(r"($\eta_{a}$)")) +
  ggtitle("Preference shock") +
  theme(
    text = element_text(size=text_size),
    axis.title.y = element_text(size=text_size),
    axis.title.x = element_text(size=text_size),
    plot.title = element_text(size=title_size, hjust=0.5)
  )
print(p1)

p2 <- ggplot(H, aes(x2,f2)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f2_gauss),linetype="dashed",linewidth=1.2) +
  ylab("") +
  xlab(TeX(r"($\eta_{e}$)")) +
  ggtitle("Cost-push shock") +
  theme(
    text = element_text(size=text_size),
    axis.title.y = element_text(size=text_size),
    axis.title.x = element_text(size=text_size),
    plot.title = element_text(size=title_size, hjust=0.5)
  )
print(p2)

p3 <- ggplot(H, aes(x3,f3)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f3_gauss),linetype="dashed",linewidth=1.2) +
  ylab("") +
  xlab(TeX(r"($\eta_{z}$)")) +
  ggtitle("Productivity shock") +
  theme(
    text = element_text(size=text_size),
    axis.title.y = element_text(size=text_size),
    axis.title.x = element_text(size=text_size),
    plot.title = element_text(size=title_size, hjust=0.5)
  )
print(p3)

p4 <- ggplot(H, aes(x4,f4)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f4_gauss),linetype="dashed",linewidth=1.2) +
  ylab("") +
  xlab(TeX(r"($\eta_{r}$)")) +
  ggtitle("Monetary policy shock") +
  theme(
    text = element_text(size=text_size),
    axis.title.y = element_text(size=text_size),
    axis.title.x = element_text(size=text_size),
    plot.title = element_text(size=title_size, hjust=0.5)
  )
print(p4)

grid.arrange(p1,p2,p3,p4,ncol=2)

fig_4_top    <- (p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_4_bottom <- (p3 + p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_4 <- (fig_4_top / fig_4_bottom)
print(fig_4)
ggsave(filename = "plots/R/fig_4.pdf", plot = fig_4, dpi = 300, units = "px", width = 3000, height = 2000)
