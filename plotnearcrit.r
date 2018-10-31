setwd('/Users/grayson/evosir')
rm(list = ls())
library(ggplot2)
simdata <- read.csv('data/exptime.csv')
simdata <- simdata[simdata$rho == 4,]
simdata <- simdata[simdata$n == 100000,]
simdata <- simdata[(simdata$lambda < 1.4) & (simdata$lambda > 1.1),]
ggplot() +
  geom_point(aes(simdata$lambda,1-simdata$S)) +
  geom_vline(xintercept = 1.25) +
  labs(y = "Size of Epidemic", x = expression(lambda)) +
  title("Constant Time Recovery Size of Epidemic vs lambda")
rm(list = ls())
simdata <- read.csv('data/consttime.csv')
simdata <- simdata[simdata$rho == 4,]
simdata <- simdata[simdata$n == 100000,]
simdata <- simdata[(simdata$lambda < 3) & (simdata$lambda > 0),]
ggplot() +
  geom_point(aes(simdata$lambda,1-simdata$S)) +
  geom_vline(xintercept = 1.0084) +
  labs(y = "Size of Epidemic", x = expression(lambda)) +
  title("Exponential Time Recovery Size of Epidemic vs lambda")