#loading libraries ####
library(ggplot2)
library(car)
library(MASS)

# Loading Dataset ####
setwd("C:/Users/comet/Desktop/ATS/ATS-coursework")
load("C:/Users/comet/Desktop/ATS/ATS-coursework/fMRI-ROI-time-series.RData")

#Y[j,,i,d]
#scan d from subject j, on region i. 404 measurements at lags of 14000ms
#d=1,2
#j=1,...,70
j <- sample(1:70, 2)
#i=1,...,24
i <- sample(1:24, 2)

l <- 404

y1 <- Y[j[1],,i[1],1]
y2 <- Y[j[2],,i[2],1]

var(y1)
var(y2)

ts.plot(y1)
ts.plot(y2, col = "red")

hist(y1, breaks = 15)
hist(y2, breaks = 15)

theta_0 <- c(.5,2,2)
theta_mle_1 <- estimator(y1, theta_0)
theta_mle_2 <- estimator(y2, theta_0)

theta_mle_1
theta_mle_2

KF1 <- KF(y1, 
          theta_mle_1$theta_list$hat_phi,
          theta_mle_1$theta_list$hat_sigma_e,
          theta_mle_1$theta_list$hat_sigma_eta)
KF1

KF2 <- KF(y2, 
          theta_mle_2$theta_list$hat_phi,
          theta_mle_2$theta_list$hat_sigma_e,
          theta_mle_2$theta_list$hat_sigma_eta)
KF2

v1 <- y1 - KF1$mu_pred
v2 <- y2 - KF2$mu_pred

ts.plot(v1)
lines(v2, col = "red")

par(mfrow = c(1,2))
hist(v1, breaks = 15, freq = F)
hist(v2, breaks = 15, add = T, density = 30, col = "red", freq =F)

qqPlot(v1)
qqPlot(v2)

acf(v1)
acf(v2)

kurtosis(v1)
kurtosis(v2)

shapiro.test(v1)
shapiro.test(v2)

