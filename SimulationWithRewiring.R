setwd('~/Desktop/evoSIR-DOmath')
library(ggplot2)
simulationdata <- read.csv('data/realexptimerewire.csv')
realtimesimgraphrewire <- ggplot(simulationdata, aes(t)) + 
  stat_smooth(aes(y = s), method = "loess", size = 1.5, span=0.001, linetype="dotted", color="black") +
  stat_smooth(aes(y = i), method = "loess", size = 1.5, span=0.001, linetype="solid", color="black") +
  stat_smooth(aes(y = r), method = "loess", size = 1.5, span=0.001, linetype="dashed", color="black") +
  theme_bw() +
  labs(y = NULL)
ggsave("realtimesimgraphrewire.png", plot = realtimesimgraphrewire, device="png",path = "~/Desktop/evoSIR-DOmath/paper_plots")