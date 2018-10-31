#setwd('~/Desktop/evoSIR-DOmath')
setwd('~/evosir')
library(ggplot2)
rm(list = ls())
if(FALSE){
#DIFFEQ WITH REWIRING
#observed epidemic size over time with rewiring
theorydata <- read.csv('theory/ODE3.csv')
realtimeode3rewire <- ggplot(theorydata, aes(t)) + 
  stat_smooth(aes(y = S), method = "loess", size = 1, span=0.05, linetype="dotted", color="black") +
  stat_smooth(aes(y = I), method = "loess", size = 1, span=0.05, linetype="solid", color="black") +
  stat_smooth(aes(y = R), method = "loess", size = 1, span=0.05, linetype="dashed", color="black") +
  theme_bw() +
  xlim(0,6)
  labs(y = NULL)
ggsave("paper_plots/realtimeode3.png", plot = realtimeode3rewire, device="png")#,path = "~/Desktop/evoSIR-DOmath/paper_plots")
#theoretical epidemic size over time with rewiring
simulationdata <- read.csv('data/realexptime.csv')
realtimesimgraphrewire <- ggplot(simulationdata, aes(t)) + 
  stat_smooth(aes(y = s), method = "loess", size = 1, span=0.001, linetype="dotted", color="black") +
  stat_smooth(aes(y = i), method = "loess", size = 1, span=0.001, linetype="solid", color="black") +
  stat_smooth(aes(y = r), method = "loess", size = 1, span=0.001, linetype="dashed", color="black") +
  theme_bw() +
  xlim(0,6)
  labs(y = NULL)
ggsave("paper_plots/realtimesimgraph.png", plot = realtimesimgraphrewire, device="png")#,path = "~/Desktop/evoSIR-DOmath/paper_plots")
#END DIFFEQ WITH REWIRING
}
#observed epidemic size at crit val
simulationdata <- read.csv('data/critvalsim.csv')
simulationdata$s = simulationdata$s/simulationdata$n
simulationdata$i = simulationdata$i/simulationdata$n
simulationdata$r = simulationdata$r/simulationdata$n
realtimeode3rewire <- ggplot(simulationdata, aes(t)) + 
  stat_smooth(aes(y = s), method = "loess", size = 1, span=0.001, linetype="dotted", color="black") +
  stat_smooth(aes(y = i), method = "loess", size = 1, span=0.001, linetype="solid", color="black") +
  stat_smooth(aes(y = r), method = "loess", size = 1, span=0.001, linetype="dashed", color="black") +
  theme_bw() +
  xlim(0,30) +
  labs(y = NULL)
ggsave("paper_plots/critval.png", plot = realtimeode3rewire, device="png")#,path = "~/Desktop/evoSIR-DOmath/paper_plots")
#observed epidemic size over time without rewiring
theorydata <- read.csv('theory/ODE3.csv')
realtimeode3rewire <- ggplot(theorydata, aes(t)) + 
  stat_smooth(aes(y = S), method = "loess", size = 1, span=0.05, linetype="dotted", color="black") +
  stat_smooth(aes(y = I), method = "loess", size = 1, span=0.05, linetype="solid", color="black") +
  stat_smooth(aes(y = R), method = "loess", size = 1, span=0.05, linetype="dashed", color="black") +
  theme_bw() +
  xlim(0,8) +
  labs(y = NULL)
ggsave("paper_plots/realtimeode3.png", plot = realtimeode3rewire, device="png")#,path = "~/Desktop/evoSIR-DOmath/paper_plots")
#theoretical epidemic size over time
simulationdata <- read.csv('data/realexptime.csv')
simulationdata$s = simulationdata$s/simulationdata$n
simulationdata$i = simulationdata$i/simulationdata$n
simulationdata$r = simulationdata$r/simulationdata$n
realtimesimgraphrewire <- ggplot(simulationdata, aes(t)) + 
  stat_smooth(aes(y = s), method = "loess", size = 1, span=0.01, linetype="dotted", color="black") +
  stat_smooth(aes(y = i), method = "loess", size = 1, span=0.01, linetype="solid", color="black") +
  stat_smooth(aes(y = r), method = "loess", size = 1, span=0.01, linetype="dashed", color="black") +
  theme_bw() +
  xlim(0,8) +
  labs(y = NULL)
ggsave("paper_plots/realtimesimgraph.png", plot = realtimesimgraphrewire, device="png")#,path = "~/Desktop/evoSIR-DOmath/paper_plots")
#epidemic size vs lambda constant time
plot.new()
simdata <- read.csv('data/consttime.csv')
simdata <- simdata[simdata$rho == 4,]
simdata <- simdata[simdata$n == 100000,]
simdata <- simdata[(simdata$lambda < 1.4) & (simdata$lambda > .5),]
ode <- read.csv('theory/constrho4.csv')
lowerbound <- read.csv('theory/constrho4lower.csv')
consttimelambda <- ggplot() +
  geom_point(aes(simdata$lambda,1-simdata$S)) +
  geom_vline(xintercept = 1.0084) +
  labs(y = "Size of Epidemic", x = expression(lambda)) +
  title("Constant Time Recovery Size of Epidemic vs lambda") +
  geom_line(aes(lowerbound$lambda,1-lowerbound$z),linetype = "solid") +
  geom_line(aes(ode$lambda,ode$t), linetype = "dotted") +
  xlim(.5,1.4) +
  theme_bw()
ggsave("paper_plots/consttimelambda.png", plot = consttimelambda, device="png")#,path = "~/Desktop/evoSIR-DOmath/paper_plots")
#epidemic size vs lambda exponential time
simdata <- read.csv('data/exptime.csv')
simdata <- simdata[simdata$rho == 4,]
simdata <- simdata[simdata$n == 100000,]
simdata <- simdata[(simdata$lambda < 2) & (simdata$lambda > 0.0),]
lowerbound <- read.csv('theory/exprho4genfunc.csv')
ode <- read.csv('theory/exprho4size.csv')
newode <- read.csv('theory/ode3withrewiring.csv')
newode <- newode[newode$rho == 4,]
#newode <- head(newode,-2)
exptimelambda <- ggplot() +
  geom_point(aes(simdata$lambda,1-simdata$S)) +
  geom_vline(xintercept = 1.25) +
  labs(y = "Size of Epidemic", x = expression(lambda)) +
  title("Exponential Time Recovery Size of Epidemic vs lambda") +
  geom_line(aes(lowerbound$lambda,1-lowerbound$z),linetype = "solid") +
  geom_line(aes(ode$lambda,ode$t), linetype = "dotted") +
  geom_line(aes(newode$lambda,newode$R), linetype = "dashed") +
  theme_bw()+
  xlim(0,2)
ggsave("paper_plots/exptimelambda.png", plot = exptimelambda, device="png")#,path = "~/Desktop/evoSIR-DOmath/paper_plots")
