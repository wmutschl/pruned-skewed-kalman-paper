# Replicates Figure 3: "Probability density functions and cumulative distribution
# functions of a CSN distributed random variable with two skewness dimensions
# and the approximating pruned CSN distribution with one skewness dimenstions"
# of the paper
# "Pruned Skewed Kalman Filter and Smoother: With Application to the Yield Curve"
# by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
# =========================================================================
library(ggplot2)
library(gridExtra)
library(csn)
library(latex2exp)
source("R/csn_prune_params.R")

Sigma_gauss <- diag(c(0.030168028450613,0.000247656247964,0.008864707956171,0.002790305261243))

mu_eta <- matrix(c(0.020056896328599,0.000163664051014,0.010320181911290,-0.003493511614978),4,1)
Sigma_eta <- diag(c(0.001042540680430,0.000000068274142,0.000169265857962,0.000020243790347))
Gamma_eta <- diag(c(-38.418941358883224,-4850.019235169948,-7089.714197084982,939.5377937249717))
nu_eta <- matrix(c(0,0,0,0),4,1)
Delta_eta <- diag(c(1,1,1,1))

x1 <- seq(-0.1,0.1,length=200)
x2 <- seq(-0.001,0.001,length=200)
x3 <- seq(-0.03,0.03,length=200)
x4 <- seq(-0.01,0.01,length=200)
H <- data.frame(x=x,
                f1_gauss=dnorm(x1,0,Sigma_gauss[1,1]),
                f2_gauss=dnorm(x2,0,Sigma_gauss[2,2]),
                f3_gauss=dnorm(x3,0,Sigma_gauss[3,3]),
                f4_gauss=dnorm(x4,0,Sigma_gauss[4,4]),
                f1=dcsn(x1,mu_eta[1,1],Sigma_eta[1,1],Gamma_eta[1,1],nu_eta[1,1],Delta_eta[1,1]),
                f2=dcsn(x2,mu_eta[2,1],Sigma_eta[2,2],Gamma_eta[2,2],nu_eta[2,1],Delta_eta[2,2]),
                f3=dcsn(x3,mu_eta[3,1],Sigma_eta[3,3],Gamma_eta[3,3],nu_eta[3,1],Delta_eta[3,3]),
                f4=dcsn(x4,mu_eta[4,1],Sigma_eta[4,4],Gamma_eta[4,4],nu_eta[4,1],Delta_eta[4,4])
)

p1 <- ggplot(H, aes(x1,f1)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f1_gauss),linetype="dashed",linewidth=1.2) +
  ylab("density") +
  xlab(TeX(r"($\eta_{a}$)")) +
  theme(text = element_text(size = 20))
p1
pdf(".../Paper/plots/dsge_pdf_eta_a.pdf")

p2 <- ggplot(H, aes(x2,f2)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f2_gauss),linetype="dashed",linewidth=1.2) +
  ylab("density") +
  xlab(TeX(r"($\eta_{e}$)")) +
  theme(text = element_text(size = 20))
p2
pdf(".../Paper/plots/dsge_pdf_eta_e.pdf")

p3 <- ggplot(H, aes(x3,f3)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f3_gauss),linetype="dashed",linewidth=1.2) +
  ylab("density") +
  xlab(TeX(r"($\eta_{z}$)")) +
  theme(text = element_text(size = 20))
p3
pdf(".../Paper/plots/dsge_pdf_eta_z.pdf")

p4 <- ggplot(H, aes(x4,f4)) +
  geom_line(linewidth=1.2) +
  geom_line(aes(y=f4_gauss),linetype="dashed",linewidth=1.2) +
  ylab("density") +
  xlab(TeX(r"($\eta_{r}$)")) +
  theme(text = element_text(size = 20))
p4
pdf(".../Paper/plots/dsge_pdf_eta_r.pdf")

grid.arrange(p1,p2,p3,p4,ncol=2)

