rm(list=ls())
library(ggplot2)

#Preliminary Analysis on the presumed gaussian####
load("fMRI-ROI-time-series.RData")
y1 <- Y[59,,16,1]
df1 <- data.frame(time = seq_along(y1), value = as.numeric(y1))

fit1 <- lowess(y1, f = 3/4)$y
fit1 <- loess( y1 ~ time, data = data.frame(y1, time =  1:length(y1)), family = "symmetric",
               degree = 1)$fitted

df_sd_y1<-data.frame(time = seq_along(y1), fit = fit1, value=sqrt(lowess((y1 - fit1)^2, f = 3/4)$y)) %>% 
  mutate(ub1 = fit1 + 2*value, lb1 = fit1 - 2*value)
df_sd_y1<-data.frame(time = seq_along(y1), fit = fit1, value=sqrt(loess(sqrt(devs) ~ time, 
                                                                        data = data.frame(devs = (y1 - fit1)^2, time = 1:length(y1)),,
                                                                        degree = 1, family = "symmetric")$fitted)) %>% 
  mutate(ub1 = fit1 + 2*value, lb1 = fit1 - 2*value)

g1<-ggplot(df1, aes(x = time, y = value)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_line(data=df_sd_y1, aes(x = time, y = fit1), color = "green", linewidth = 0.75)+
  geom_line(data=df_sd_y1, aes(x = time, y = ub1), color = "purple", linewidth = 0.75)+
  geom_line(data=df_sd_y1, aes(x = time, y = lb1), color = "purple", linewidth = 0.75)+ 
  theme_minimal() + labs(x = "Time Index", y = "y2")
g1

# from the plot the time series appears stationary for the mean and for the standard deviation, consequently
# for the variance too
# the green line represents the local means and we can grasp from the image that is fluctuating
# around zero
# the purple lines are the upper and lower for a 95% CI for the mean.
# since they are nearly flat, the standard deviation doesn't change too much over time

#table of summary statistics
xbar <- mean(y1) # very close to 0, suggesting that the series is stationary wrt mean
var <- var(y1) # 1.69,
kur <- kurtosis(y1) # close to 3, ok!
skew <- skewness(y1) # again close to 0
#The value of the sample mean suggests that the underlying process is in fact zero-mean.
#from the values of kurtosis close to 3 and skewness close to zero we are hints that
#the time series might follow a normal distribution

acf(y1)
pacf(y1)
# The acf graph of our process shows an exponential decay of the spikes, the correlations
# start high at lag 1 (around 0.6) and trail off toward zero as the lags increase.
# The pacf has a significant spike at lag 1 but then at lag 2 the decay is fast to 0.
# These graphs are compatible both with an AR(1) or ARMA(1,1).

qqnorm(y1)
qqline(y1)
#there is no an evident departure of the sample quantiles from the values of the theoretical ones

shapiro_result <- shapiro.test(y1)
print(shapiro_result)
#H0: population is normally distributed
#test confirms that we do not reject the null hypothesis corresponding to the normality assumption

#############################################################################################
#############################################################################################

y2 <- Y[59,,13,1]

df2 <- data.frame(time = seq_along(y2), value = as.numeric(y2))

fit2 <- lowess(y2, f = 3/4)$y
fit2 <- loess( y2 ~ time, data = data.frame(y2, time =  1:length(y2)),
               family = "symmetric", degree = 1)$fitted

df_sd_y2<-data.frame(time = seq_along(y2), fit = fit2, value=sqrt(lowess((y2 - fit2)^2, f = 3/4)$y)) %>% 
  mutate(ub2 = fit2 + 2*value, lb2 = fit2 - 2*value)
df_sd_y2<-data.frame(time = seq_along(y2), fit = fit2, value=sqrt(loess(sqrt(devs) ~ time, 
                      data = data.frame(devs = (y2 - fit2)^2, time = 1:length(y2)), family = "symmetric",
                      degree = 1)$fitted)) %>% 
                      mutate(ub2 = fit2 + 2*value, lb2 = fit2 - 2*value)

g2<-ggplot(df2, aes(x = time, y = value)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_line(data=df_sd_y2, aes(x = time, y = fit2), color = "green", linewidth = 0.75)+
  geom_line(data=df_sd_y2, aes(x = time, y = ub2), color = "purple", linewidth = 0.75)+
  geom_line(data=df_sd_y2, aes(x = time, y = lb2), color = "purple", linewidth = 0.75)+ 
  theme_minimal() + labs(x = "Time Index", y = "y2")
g2

#again the time series appears stationary for the mean while the sd appears to very slight increase over 
#time but this is mostly due to the presence of an outlier observation

xbar <- mean(y2) # very close to 0, suggesting that the series is stationary wrt mean
var <- var(y2) # 11.82
kur <- kurtosis(y2) # close to 3, ok!
skew <- skewness(y2) # again close to 0
#The value of the sample mean suggests that the underlying process is in fact zero-mean.
#the value of kurtosis is a bit higher than the one of a normal while skewness is still close to zero 


acf(y2)
pacf(y2)
#The acf graph of our process shows an exponential decay of the spikes, the correlations 
#start high at lag 1 (around 0.7) and trail off toward zero as the lags increase; the pace is, however,
#slower than in the previous series.
#The pacf has a significant spike at lags 1 and 2 but then it decays to 0.
#These graphs are compatible suggest that the underlying process could be an ARMA(1,1).

qqnorm(y2)
qqline(y2)
#plot shows that there might be problems with the normality assumption 

shapiro_result_2 <- shapiro.test(y2)
print(shapiro_result_2)
#H0: population is normally distributed
#from the result of the test we in fact reject the assumption of normality with respect to 
#the second time series

################################################
#log-scaled transformation

y2_transf<- diff(y2)
df2_transf <- data.frame(time = seq_along(y2), value = as.numeric(y2_transf))

fit2_transf <- lowess(y2_transf, f = 1/3)$y
df_sd_y2_transf<-data.frame(time = seq_along(y2), fit = fit2_transf, value=sqrt(lowess((y2_transf - fit2_transf)^2, 
                            f = 1/3)$y)) %>% mutate(ub2_transf = fit2_transf + 2*value, lb2_transf = fit2_transf - 2*value)

g2<-ggplot(df2_transf, aes(x = time, y = value)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_line(data=df_sd_y2_transf, aes(x = time, y = fit2_transf), color = "green", linewidth = 0.75)+
  geom_line(data=df_sd_y2_transf, aes(x = time, y = ub2_transf), color = "purple", linewidth = 0.75)+
  geom_line(data=df_sd_y2_transf, aes(x = time, y = lb2_transf), color = "purple", linewidth = 0.75)+ 
  theme_minimal() + labs(x = "Time Index", y = "y2_transf")
g2
