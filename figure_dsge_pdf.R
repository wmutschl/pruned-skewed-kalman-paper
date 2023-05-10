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
source("R/csn_prune_params.R")
savedev <- function(OS,file){
  if(OS=="linux" | OS=="mac"){
    dev.copy2pdf(file=file)
  }
  if(OS=="windows") {
    savePlot(filename=file,"pdf")
  }
}
params <- readMat("results/results_ireland2004_params_for_pdfs.mat")

x1 <- seq(-0.1,0.1,length=200)
x2 <- seq(-0.001,0.001,length=200)
x3 <- seq(-0.03,0.03,length=200)
x4 <- seq(-0.01,0.01,length=200)
H <- data.frame(f1_gauss=dnorm(x1,0,sqrt(params$Sigma.gauss[1,1])),
                f2_gauss=dnorm(x2,0,sqrt(params$Sigma.gauss[2,2])),
                f3_gauss=dnorm(x3,0,sqrt(params$Sigma.gauss[3,3])),
                f4_gauss=dnorm(x4,0,sqrt(params$Sigma.gauss[4,4])),
                f1=dcsn(x1,params$mu.eta[1,1],params$Sigma.eta[1,1],params$Gamma.eta[1,1],params$nu.eta[1,1],params$Delta.eta[1,1]),
                f2=dcsn(x2,params$mu.eta[2,1],params$Sigma.eta[2,2],params$Gamma.eta[2,2],params$nu.eta[2,1],params$Delta.eta[2,2]),
                f3=dcsn(x3,params$mu.eta[3,1],params$Sigma.eta[3,3],params$Gamma.eta[3,3],params$nu.eta[3,1],params$Delta.eta[3,3]),
                f4=dcsn(x4,params$mu.eta[4,1],params$Sigma.eta[4,4],params$Gamma.eta[4,4],params$nu.eta[4,1],params$Delta.eta[4,4])
)

p1 <- ggplot(H, aes(x1,f1)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f1_gauss),linetype="dashed",linewidth=1.2) +
  ylab("density") +
  xlab(TeX(r"($\eta_{a}$)")) +
  theme(text = element_text(size = 20))
p1
savedev("linux",file="../Paper/plots/dsge_pdf_eta_a.pdf")

p2 <- ggplot(H, aes(x2,f2)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f2_gauss),linetype="dashed",linewidth=1.2) +
  ylab("density") +
  xlab(TeX(r"($\eta_{e}$)")) +
  theme(text = element_text(size = 20))
p2
savedev("linux",file="../Paper/plots/dsge_pdf_eta_e.pdf")

p3 <- ggplot(H, aes(x3,f3)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f3_gauss),linetype="dashed",linewidth=1.2) +
  ylab("density") +
  xlab(TeX(r"($\eta_{z}$)")) +
  theme(text = element_text(size = 20))
p3
savedev("linux",file="../Paper/plots/dsge_pdf_eta_z.pdf")

p4 <- ggplot(H, aes(x4,f4)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f4_gauss),linetype="dashed",linewidth=1.2) +
  ylab("density") +
  xlab(TeX(r"($\eta_{r}$)")) +
  theme(text = element_text(size = 20))
p4
savedev("linux",file="../Paper/plots/dsge_pdf_eta_r.pdf")

grid.arrange(p1,p2,p3,p4,ncol=2)
