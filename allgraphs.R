rm(list = ls())
setwd('/Users/grayson/evosir')
exptime <- read.csv('data/exptime.csv')
lambdadata <- exptime[exptime$rho == 4,]
lambdadata <- lambdadata[lambdadata$n == 100000,]
#PROPORTION OF POPULATION INFECTED EXPONENTIAL TIME RHO = 4
theory <- read.csv('theory/exprho4size.csv')
lowerbound <- read.csv('theory/exprho4genfunc.csv')
ode <- read.csv('theory/exprho4improvedode.csv')
newode <- read.csv('theory/ode3withrewiring.csv')
newode <- newode[newode$rho == 4,]
newode <- head(newode,-2)
png(filename = 'plots/exptimerho4size.png')
plot(lambdadata$lambda,1-lambdadata$S, xlab = expression(lambda), ylab = '',main = 'Proportion of Population Infected in Epidemic \n (Exponential Infection Time, rho= 4)')
lines(theory$lambda,theory$t,lty = 2)
lines(lowerbound$lambda,1-lowerbound$z,lty = 3)
#lines(ode$lambda,ode$R,col = 'green')
#lines(newode$lambda,newode$R,col = 'purple')
abline(v = 1.25)
dev.off()
#PROBABILITY OF EPIDEMIC GRAPH EXPONENTIAL TIME RHO = 4
epidemics <- lambdadata[lambdadata$S < .70,]
theory <- read.csv('theory/exprho4prob.csv')
numsurvive <- table(epidemics$lambda)
num <- as.numeric(names(numsurvive))
totalsims <- table(lambdadata$lambda)
fracsurvive <- c()
for (i in names(numsurvive))
{
  fracsurvive <- c(fracsurvive,numsurvive[i]/totalsims[i])
}
png(filename = 'plots/exptimerho4prob.png')
plot(num,fracsurvive,main = 'Probability of Epidemic\n (Exponential Infection Time, rho= 4)',xlab = expression(lambda), ylab = '')
lines(theory$lambda,1-theory$u,lty = 2)
dev.off()
#PROPORTION OF POPULATION INFECTED EXPONENTIAL TIME RHO = 0
norewire <- exptime[exptime$rho == 0,]
theory <- read.csv('theory/exprho0size.csv')
lowerbound <- read.csv('theory/exprho0genfunc.csv')
png(filename = 'plots/exptimerho0size.png')
plot(norewire$lambda,1-norewire$S, xlab = expression(lambda), ylab = '',main = 'Proportion of Population Infected in Epidemic \n (Exponential Infection Time, rho= 0)')
lines(theory$lambda,theory$t,lty = 2)
lines(lowerbound$lambda,1-lowerbound$z,lty = 3)
dev.off()
#PROBABILITY OF EPIDEMIC EXPONENTIAL TIME RHO = 0
theory <- read.csv('theory/exprho0prob.csv')
epidemics <- norewire[norewire$S < .90,]
numsurvive <- table(epidemics$lambda)
num <- as.numeric(names(numsurvive))
totalsims <- table(norewire$lambda)
fracsurvive <- c()
for (i in names(numsurvive))
{
  fracsurvive <- c(fracsurvive,numsurvive[i]/totalsims[i])
}
png(filename = 'plots/exptimerho0prob.png')
plot(num,fracsurvive,main = 'Probability of Epidemic\n (Exponential Infection Time, rho= 0)',xlab = expression(lambda), ylab = '')
lines(theory$lambda,1-theory$u,lty = 2)
dev.off()
#PROPORTION OF POPULATION INFECTED EXPONENTIAL TIME LAMBDA = 1
rhodata <- exptime[exptime$n == 100000,]
rhodata <- rhodata[rhodata$lambda == 1,]
epidemics <- rhodata[rhodata$S < .90,]
theory <- read.csv('theory/explambda1size.csv')
newtheory <- read.csv('theory/ode3withrewiring.csv')
newtheory <- newtheory[newtheory$rho < 4,]
#newtheory <- tail(newode,2)
lowerbound <- read.csv('theory/explambda1genfunc.csv')
png(filename = 'plots/exptimelambda1size.png')
plot(rhodata$rho,1-rhodata$S, xlab = 'Rho', ylab = '',main = 'Proportion of Population Infected in Epidemic\n(Exponential Infection Time, lambda = 1)',ylim = c(0,1))
lines(theory$rho,theory$t,lty = 2)
lines(lowerbound$rho,1-lowerbound$z,lty = 3)
#lines(newtheory$rho,newtheory$R,col = 'purple')
dev.off()
numsurvive <- table(epidemics$rho)
num <- as.numeric(names(numsurvive))
totalsims <- table(rhodata$rho)
fracsurvive <- c()
for (i in names(numsurvive))
{
  fracsurvive <- c(fracsurvive,numsurvive[i]/totalsims[i])
}
#PROBABILITY OF EPIDEMIC GRAPH EXPONENTIAL TIME LAMBDA = 1
theory <- read.csv('theory/explambda1prob.csv')
png(filename = 'plots/exptimelambda1prob.png')
plot(num,fracsurvive,main = 'Probability of Epidemic\n (Exponential Infection Time, lambda = 1)',xlab = expression(rho),ylab = '')
lines(theory$rho,1-theory$u,lty = 2)
dev.off()
#CONSTANT TIME PROPORTION OF POPULATION INFECTED IN EPIDEMIC RHO = 4
consttime <- read.csv('data/consttime.csv')
theory <- read.csv('theory/constrho4.csv')
lowerbound <- read.csv('theory/constrho4lower.csv')
lambdadata <- consttime[consttime$rho == 4,]
lambdadata <- lambdadata[lambdadata$n == 100000,]
epidemics <- lambdadata[lambdadata$S < .90,]
png(filename = 'plots/consttimerho4size.png')
plot(lambdadata$lambda,1-lambdadata$S, xlab = expression(lambda), ylab = 'Proportion of Population Infected \n (Constant Time, rho = 4)',main = 'Proportion of Population Infected in Epidemic \n (Constant Infection Time, rho = 4)')
lines(theory$lambda,theory$t,lty = 2)
lines(lowerbound$lambda,1-lowerbound$z,lty = 3)
dev.off()
numsurvive <- table(epidemics$lambda)
num <- as.numeric(names(numsurvive))
totalsims <- table(lambdadata$lambda)
fracsurvive <- c()
for (i in names(numsurvive))
{
  fracsurvive <- c(fracsurvive,numsurvive[i]/totalsims[i])
}
#CONSTANT TIME PROBABILITY OF EPIDEMIC RHO = 4
theory <- read.csv('theory/constrho4lower.csv')
png(filename = 'plots/consttimerho4prob.png')
plot(num,fracsurvive,xlab = expression(lambda),ylab = 'Probability of Epidemic', main = 'Probability of Epidemic \n (Constant Time, rho = 4)')
lines(theory$lambda,1-theory$z,lty = 2)
dev.off()
#CONSTANT TIME PROPORTION OF POPULATION INFECTED IN EPIDEMIC LAMBDA = 1
rhodata <- consttime[consttime$lambda == 1,]
rhodata <- rhodata[rhodata$n == 100000,]
epidemics <- rhodata[rhodata$S < .90,]
theory <- read.csv('theory/constlambda1.csv')
lowerbound <- read.csv('theory/constlambda1lower.csv')
png(filename = 'plots/consttimelambda1size.png')
plot(rhodata$rho,1-rhodata$S, xlab = expression(rho), ylab = 'proportion of population who survive',main = 'Proportion of Population Infected in Epidemic \n and Probability of Epidemic \n (Constant Infection Time, lambda = 1)')
lines(theory$rho,theory$t,lty = 2)
lines(lowerbound$rho,1-lowerbound$z,lty = 3)
dev.off()
numsurvive <- table(epidemics$rho)
num <- as.numeric(names(numsurvive))
totalsims <- table(rhodata$rho)
fracsurvive <- c()
for (i in names(numsurvive))
{
  fracsurvive <- c(fracsurvive,numsurvive[i]/totalsims[i])
}
#CONSTANT TIME PROBABILITY OF EPIDEMIC LAMBDA = 1
theory <- read.csv('theory/constlambda1lower.csv')
png(filename = 'plots/consttimelambda1prob.png')
plot(num,fracsurvive,xlab = expression(rho),main = 'Probability of Survival \n (Constant Time, lambda = 1)', ylab = 'Probability of Survival')
lines(theory$rho,1-theory$z,lty = 2)
dev.off()
#PROPORTION OF POPULATION INFECTED CONSTANT TIME RHO = 0
norewire <- consttime[consttime$rho == 0,]
theory <- read.csv('theory/constrho0.csv')
lowerbound <- read.csv('theory/constrho0lower.csv')
png(filename = 'plots/consttimerho0size.png')
plot(norewire$lambda,1-norewire$S, xlab = expression(lambda), ylab = '',main = 'Proportion of Population Infected in Epidemic \n (Constant Infection Time, rho= 0)')
lines(theory$lambda,theory$t,lty = 2)
lines(lowerbound$lambda,1-lowerbound$z,lty = 3)
dev.off()
#REAL TIME EXPONENTIAL GRAPH
#Simulation
simulationdata <- read.csv('data/realexptime.csv')
png(filename = 'plots/realtimesimgraph.png')
plot(simulationdata$t,simulationdata$r/simulationdata$n,type = 'l',xlab = 'time', ylab = 'population proportion',main = 'Simulation Data',xlim = c(0,10),ylim = c(0,1))
lines(simulationdata$t,simulationdata$i/simulationdata$n, lty = 2)
lines(simulationdata$t,simulationdata$s/simulationdata$n, lty = 3)
dev.off()
#Simulation with rewiring
simulationdata <- read.csv('data/realexptimerewire.csv')
png(filename = 'plots/realtimesimgraphrewire.png')
normalizeddata <- data.frame(simulationdata$r/simulationdata$n,simulationdata$i/simulationdata$n,simulationdata$s/simulationdata$n,simulationdata$t)
colnames(normalizeddata) = c('r','i','s','t')
downsamplevector <- runif(nrow(normalizeddata)) < .001
normalizeddata <- normalizeddata[downsamplevector,]
#lo <- loess(r~t,data = normalizeddata,span = .001)
#plot(normalizeddata$t,predict(lo),type = 'l')
plot(normalizeddata$t,normalizeddata$s,type = 'l')
lines(simulationdata$t,simulationdata$i/simulationdata$n, lty = 3)
#lo <- loess(i~t,data = normalizeddata,span = .001)
#lines(normalizeddata$t,predict(lo),lty = 2)
#lines(normalizeddata$t,normalizeddata$r,lty = 2)
lines(simulationdata$t,simulationdata$r/simulationdata$n, lty = 2)
#lo <- loess(r~t,data = normalizeddata,span = .001)
#lines(normalizeddata$t,predict(lo),lty = 3)
#lines(simulationdata$t,(simulationdata$mu-2)/3,col = 'purple')
dev.off()
#old theory
png(filename = 'plots/realtimeoneedgetheory.png')
plot(theory$t,theory$R,type = 'l',xlab = 'time', ylab = 'population proportion', main = 'One Edge Only Diffeq',xlim = c(0,10),ylim = c(0,1))
lines(theory$t,theory$I,lty = 2)
lines(theory$t,theory$S,lty = 3)
dev.off()
#new theory
png(filename = 'plots/realtimemanyedgetheory.png')
theory <- read.csv('theory/realtimeexptheoryimproved.csv')
plot(theory$t,theory$R,type = 'l',xlab = 'time', ylab = 'population proportion', main = 'All Edges Diffeq',xlim = c(0,10),ylim = c(0,1))
lines(theory$t,theory$I,lty = 2)
lines(theory$t,theory$S,lty = 3)
dev.off()
#Miller 11 real time
png(filename = 'plots/realtimemiller11theory.png')
theory <- read.csv('theory/miller11.csv')
plot(theory$t,theory$R,type = 'l',xlab = 'time', ylab = 'population proportion', main = 'Miller 2011 Diffeq',xlim = c(0,10),ylim = c(0,1))
lines(theory$t,theory$I,lty = 2)
lines(theory$t,theory$S,lty = 3)
dev.off()
#ODE 3 real time
png(filename = 'plots/realtimeode3.png')
theory <- read.csv('theory/ode3.csv')
plot(theory$t,theory$R,type = 'l',xlab = 'time', ylab = 'population proportion', main = 'ODE 3',xlim = c(0,10),ylim = c(0,1))
lines(theory$t,theory$I,lty = 2)
lines(theory$t,theory$S,lty = 3)
dev.off()
#stargraph
png(filename = 'plots/stargraphrho.png')
theory <- read.csv('theory/consttimestar.csv')
simulationdata <- read.csv('data/consttimestarcomponents.csv')
lambdatheory <- theory[theory$rho==4,]
lambdadata <- simulationdata[simulationdata$rho == 4,]
lambdadata <- lambdadata[lambdadata$n == 100000,]
lambdatheory <- head(lambdatheory,-1)
plot(lambdadata$lambda,1-lambdadata$prop,xlab = expression(lambda),ylab = '')
lines(lambdatheory$lambda,lambdatheory$prop)
dev.off()
png(filename = 'plots/stargraphlambda.png')
rhotheory <- theory[theory$lambda == 1,]
rhodata <- simulationdata[simulationdata$lambda == 1,]
rhodata <- rhodata[rhodata$n == 100000,]
rhotheory <- rhotheory[-1,]
plot(rhodata$rho,1-rhodata$prop,xlab = expression(rho),ylab = '')
lines(rhotheory$rho,rhotheory$prop)
dev.off()
#Real time ODE with rewiring
png(filename = 'plots/realtimeoderewire.png')
theory <- read.csv('theory/realtimeoderewire.csv')
plot(theory$t,theory$R,type = 'l',xlab = 'time',ylab = 'population proportion', main = 'ODE WITH REWIRING',xlim = c(0,10),ylim = c(0,1))
lines(theory$t,theory$I,lty = 2)
lines(theory$t,theory$S,lty = 3)
#lines(theory$t,theory$mu-4.5,col = 'purple')
dev.off()
#simulation vs theoretical mu
png(filename = 'plots/theoreticalvssimulatedmu.png')
theory <- read.csv('theory/realtimeoderewire.csv')
simulationdata <- read.csv('data/realexptimerewire.csv')
plot(theory$t,theory$mu,type = 'l',ylim = c(0,6),xlim = c(0,6),xlab = 't',ylab = expression(mu),main = 'mu comparison over time')
lines(simulationdata$t,simulationdata$mu)
dev.off()
#no rewire mu
png(filename = 'plots/norewiremu.png')
simulationdata <- read.csv('data/realexptime.csv')
plot(simulationdata$t,simulationdata$mu,type = 'l',xlim = c(0,6),ylim = c(0,6))
dev.off()
#noodles
png(filename = 'plots/noodlen.png')
noodledata <- read.csv('data/evonoodle.csv')
noodledata <- noodledata[noodledata$rho == 10,]
plot(noodledata$n,noodledata$diam)
dev.off()
png(filename = 'plots/noodlerho.png')
noodledata <- read.csv('data/evonoodle.csv')
noodledata <- noodledata[noodledata$n == 1000,]
plot(noodledata$rho,noodledata$diam)
dev.off()
#NIV
png(filename = 'plots/newinfecteddeg.png')
simulationdata <- read.csv('data/realexptimetrackinfected.csv')
plot(simulationdata$t,simulationdata$newmu,type = 'l',xlim = c(0,6),ylim = c(0,10))
abline(h=5)
dev.off()
png(filename = 'plots/newinfecteddegnorewire.png')
simulationdata <- read.csv('data/realexptimetrackinfectednorewire.csv')
plot(simulationdata$t,simulationdata$newmu,type = 'l',xlim = c(0,6),ylim = c(0,10))
abline(h=5)
dev.off()
