# Library
library(sf)
library(ggplot2)
library(tigris)
library(dplyr)
library(DHARMa)
library(ggspatial)
library(AER)
library(corrplot)
library(aod)
library(psych)
library(GPArotation)
library(clValid)
library(cluster)
library(factoextra)
library(tidyverse)
library(car)
library(statmod)
library(cplm)
library(COMPoissonReg)
library(Rfast2)
library(flexmix)
library(DGLMExtPois)
library(MASS)
library(esquisse)

# Input Data
data <- read.csv("D:\\SKRIPSI\\cobaa.csv")
summary(data)


# Eksplorasi Data
## Eksplorasi Peubah Respon
### Boxplot Peubah Respon
boxplot(data$Y, ylab = "Jumlah Kasus Kekerasan pada Perempuan")
### Histogram Peubah Respon
hist(data$Y, main="Histogram JKP", xlab = "Jumlah Kasus Kekerasan pada Perempuan", ylab = "Frequency")

## Scatter Plot Peubah Respon dengan Setiap Peubah Penjelas
### Peubah X1
ggplot(data, aes(x = Y, y = X1)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X1")
### Peubah X2
ggplot(data, aes(x = Y, y = X2)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X2")
### Peubah X3
ggplot(data, aes(x = Y, y = X3)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X3")
### Peubah X4
ggplot(data, aes(x = Y, y = X4)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X4")
### Peubah X5
ggplot(data, aes(x = Y, y = X5)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X5")
### Peubah X6
ggplot(data, aes(x = Y, y = X6)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X6")
### Peubah X7
ggplot(data, aes(x = Y, y = X7)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X7")
### Peubah X8
ggplot(data, aes(x = Y, y = X8)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X8")
### Peubah X9
ggplot(data, aes(x = Y, y = X9)) + 
  geom_point(color= "#cc7952") + 
  geom_smooth(method = "lm", color = "darkblue") +
  theme_minimal() + 
  labs(x = "Y", y = "X9")

## Matriks Korelasi 
corr_matrix <- round(cor(data),2)
corrplot(corr_matrix, 
         type="lower",
         method = "color", 
         tl.cex = 0.5, 
         tl.col = "black",
         addCoef.col = "#2F2F2F",
         addCoefasPercent = FALSE,
         number.cex = 0.5,
         diag = FALSE)
corrplot

# Pemodelan Regresi Linear
Model <- lm(Y~., data=data)
summary(Model)

# Pengecekkan Multikolinearitas
vif(Model)

# STEPWISE SELECTION
step_mo <- step(lm(Y~ .,data=data),direction="both")
summary(step_mo)

# Pembangunan Model Alternatif
## Matriks Kebaikan Model
rmse <- matrix(NA, ncol = 1, nrow = 4)
aic <- matrix(NA, ncol = 1, nrow = 4)
bic <- matrix(NA, ncol = 1, nrow = 4)

## Pemodelan Regresi Poisson
poisson_model <- glm(Y ~ X1+X2+X5+X9, data = data, family = poisson())
a <- summary(poisson_model)
a
bic[1,1] <- BIC(poisson_model)
aic [1,1] <- AIC(poisson_model)
err_po <- residuals(poisson_model)
mse_po <- mean(err_po^2)
rmse[1,1] <- sqrt(mse_po)

## Pengecekan Overdispersi
rasio <- a$deviance/a$df.residual
rasio
x = data[,-1]
overdispreg.test(data$Y, x)

## Pemodelan Regresi COM-Poisson
com <- glm.CMP(Y ~ X1+X2+X5+X9, formula.nu = Y~1, data = data)
b <- summary(com)
b
aic[2,1] <- com$aic
bic[2,1] <- com$bic
err_com <- residuals(com)
mse_com <- mean(err_com^2)
rmse[2,1] <- sqrt(mse_com)

## Pemodelan Poisson-Tweedie
pois.tw <- cpglm(Y ~ X1+X2+X5+X9, data = data)
c <- summary(pois.tw)
c
err_tw <- residuals(pois.tw)
mse_tw <- mean(err_tw^2)
rmse[3,1] <- sqrt(mse_tw)
aic[3,1] <- c$aic
LL <- (c$aic - 2*6)/-2
bic[3,1] <- -2*LL + 6*log(length(data$Y))

## Pemodelan Negatif Binomial`
nb_model <- glm.nb(Y ~ X1+X2+X5+X9, data = data)
d <- summary(nb_model)
d
aic[4,1] <- nb_model$aic
log_likelihood <- logLik(nb_model)
num_parameters <- length(coef(nb_model))
num_observations <- nrow(data)
BIC_NB <- -2 * log_likelihood + num_parameters * log(num_observations)
BIC_NB
bic[4,1] <- BIC_NB
err_nb <- residuals(nb_model)
mse_nb <- mean(err_nb^2)
rmse[4,1] <- sqrt(mse_nb)

# Evaluasi Model
evaluasi <- cbind(rmse, aic, bic)
colnames(evaluasi) <- c("RMSE", "AIC", "BIC")
rownames(evaluasi) <- c("Poisson", "COM-Pois", "Pois-Tw", "Negative Binomial")
evaluasi
