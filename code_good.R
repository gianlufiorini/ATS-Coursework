#loading libraries ####
library(ggplot2)
library(tseries)
library(car)
library(MASS)

# Loading Dataset ####
setwd("C:/Users/comet/Desktop/ATS/ATS-coursework")
load("C:/Users/comet/Desktop/ATS/ATS-coursework/fMRI-ROI-time-series.RData")
l <- 404

#gaussian
y1 <- Y[59,,16,1]
#non-gaussian
y2 <- Y[59,,13,1]

is.na(y1) |> any()
is.na(y2) |> any()

## Plotting ####

ts.plot(y2, col = "red")
lines(y1)

## Loading Functions ####
source("ats_lab_2_functions.R")

# Estimating Kalman Filter ####
#Setting initial values
theta0 <- c(.5,2,2)

theta1_mle <- estimator(y1,theta0)
theta2_mle <- estimator(y2, theta0)

KF1 <- KF(y1, 
          theta1_mle$theta_list$hat_phi,
          theta1_mle$theta_list$hat_sigma_e,
          theta1_mle$theta_list$hat_sigma_eta)
KF1

KF2 <- KF(y2, 
          theta2_mle$theta_list$hat_phi,
          theta2_mle$theta_list$hat_sigma_e,
          theta2_mle$theta_list$hat_sigma_eta)
KF2

## Computing innovations ####
v1 <- y1 - KF1$mu_pred
v2 <- y2 - KF2$mu_pred

## Plotting ####
par(mfrow = c(1,2))
ts.plot(y1, ylim = c(min(min(y1,y2)),
                     max(max(y1,y2))))
lines(KF1$mu_pred, col = "green", lwd = 2)

ts.plot(y2, ylim = c(min(min(y1,y2)),
                     max(max(y1,y2))))
lines(KF2$mu_pred, col = "green", lwd = 2)

par(mfrow = c(1,1))
ts.plot(v1, ylim = c(min(min(v1,v2)),
                     max(max(v1,v2))))
lines(v2, ylim = c(min(min(v1,v2)),
                   max(max(v1,v2))), col = "red")

KF1$Pt1_t
KF2$Pt1_t

alpha <- .05
F_alpha_2 <- qchisq(alpha/2, df = 1)
F_1_alpha_2 <- qchisq(1 - alpha/2, df = 1)

NIS1 <- v1^2/KF1$Pt1_t
mean(NIS1)

NIS2 <- v2^2/KF2$Pt1_t
mean(NIS2)

plot(NIS1, ylim = c(0, max(max(NIS1, NIS2))))
lines(NIS2, col = "red")

KF1$llk
KF2$llk

acf(v1)
acf(v2)

shapiro.test(v1)
shapiro.test(v2)

jarque.bera.test(v1)
jarque.bera.test(v2)

ks.test(v1, y = "pnorm")
ks.test(v2, y = "pnorm")

qqPlot(v1)
qqPlot(v2)

