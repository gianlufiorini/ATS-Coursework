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
source("ats_lab_2_functions_with_mean.R")

# Estimating Kalman Filter ####
#Setting initial values
theta0 <- c(0,.5,2,2)

theta1_mle <- estimator(y1,theta0)
theta2_mle <- estimator(y2, theta0)

KF1 <- KF(y1, 
          theta1_mle$theta_list$hat_omega,
          theta1_mle$theta_list$hat_phi,
          theta1_mle$theta_list$hat_sigma_e,
          theta1_mle$theta_list$hat_sigma_eta)
KF1

KF2 <- KF(y2, 
          theta2_mle$theta_list$hat_omega,
          theta2_mle$theta_list$hat_phi,
          theta2_mle$theta_list$hat_sigma_e,
          theta2_mle$theta_list$hat_sigma_eta)
KF2

## Computing innovations ####
v1 <- y1 - KF1$mu_pred
v2 <- y2 - KF2$mu_pred

#Standardized innovations. If assumptions hold they should be IID N(0,1)
st_v1 <- v1 / KF1$Ft[-l]
st_v2 <- v2 / KF1$Ft[-l]

## Plotting ####
par(mfrow = c(1,2))
ts.plot(y1, ylim = c(min(min(y1,y2)),
                     max(max(y1,y2))))
lines(KF1$mu_pred, col = "green", lwd = 2)

ts.plot(y2, ylim = c(min(min(y1,y2)),
                     max(max(y1,y2))))
lines(KF2$mu_pred, col = "green", lwd = 2)

par(mfrow = c(1,1))
ts.plot(st_v1, ylim = c(min(min(st_v1,st_v2)),
                     max(max(st_v1,st_v2))))
lines(st_v2, ylim = c(min(min(st_v1,st_v2)),
                   max(max(st_v1,st_v2))), col = "red")
abline(h = 2, lty = 2, lwd = 1.5)
abline(h = -2, lty = 2, lwd = 1.5)

KF1$Pt1_t
KF2$Pt1_t

KF1$llk
KF2$llk

acf(st_v1, lag = 400)
acf(st_v2, lag = 400)

Box.test(st_v1, lag = 20, fitdf = 3, type = "Ljung-Box")
Box.test(st_v2, lag = 20, fitdf = 3, type = "Ljung-Box")

shapiro.test(st_v1)
shapiro.test(st_v2)

jarque.bera.test(st_v1)
jarque.bera.test(st_v2)

ks.test(st_v1, y = "pnorm")
ks.test(st_v2, y = "pnorm")

qqPlot(st_v1)
qqPlot(st_v2)

