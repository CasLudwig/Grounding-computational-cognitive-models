# C.J.H. Ludwig, June 2024
# Reward rate visualisation from expanded judgement task (to generate Figure C1 in the paper).

library(ggplot2)
library(dplyr)
library(viridis)
library(latex2exp)

rm(list=ls())
dev.off()

# Read in all the data we need
asymptoticRR <- read.csv("rewardDynamicProgramming.csv")
w5RR <- read.csv("rewardSimulation.csv")
rr_v_w <- read.csv("rewardWindowSize.csv")

# To ensure consistency with Figure 5 in the main text, convert degrees to radians
asymptoticRR$theta1 <- asymptoticRR$theta1 * pi/180
w5RR$theta1 <- w5RR$theta1 *pi/180

colScaleLims <- c(-0.065, 0.065)
xyPeak <- c(asymptoticRR$theta1[which.max(asymptoticRR$objective)], 
            asymptoticRR$theta2[which.max(asymptoticRR$objective)])

# Panel A: Asymptotic reward rate as a function of threshold height and gradient, derived through dynamic programming (cf. Malhotra et al., 2018).
A <- ggplot() +
  geom_raster(aes(x = theta1, y = theta2, fill = objective), data = asymptoticRR, interpolate = FALSE) +
  scale_fill_viridis(na.value="transparent", limits = colScaleLims, breaks=c(-0.06, -0.03, 0, 0.03, 0.06)) +
  geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyPeak[1], theta2 = xyPeak[2]), shape = 3, colour = "black", size = 2) +
  scale_x_continuous(breaks = seq(-1, 0, by = 0.25), limits = c(-1.07, 0.11)) +
  scale_y_continuous(limits = c(-1, 28), breaks = seq(0, 25, by = 5)) +
  labs(x = TeX("$\\beta$"), y = TeX("$\\alpha$"), fill = "reward rate") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill = guide_colourbar(position = "inside")) +
  theme(legend.position.inside = c(0.205, 0.68)) +
  theme(legend.key.height = unit(0.3,"cm")) +
  theme(legend.margin = margin(l = 0, r = 0, t = 0, b = 4.5)) +
  theme(legend.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 10)) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))
A

# Panel B: Reward rate as a function of threshold height and gradient, when estimated over a small window of 5 trials.
# xyNoisyPeak <- c(w5RR$theta1[which.max(w5RR$objective)], 
#                  w5RR$theta2[which.max(w5RR$objective)])
B <- ggplot() +
  geom_raster(aes(x = theta1, y = theta2, fill = objective), data = w5RR, interpolate = FALSE) +
  scale_fill_viridis(na.value="transparent", limits = colScaleLims, breaks=c(-0.06, -0.03, 0, 0.03, 0.06)) +
  geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyPeak[1], theta2 = xyPeak[2]), shape = 3, colour = "black", size = 2) +
  # geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyNoisyPeak[1], theta2 = xyNoisyPeak[2]), shape = 17, colour = "black", size = 2) +
  scale_x_continuous(breaks = seq(-1, 0, by = 0.25), limits = c(-1.07, 0.11)) +
  scale_y_continuous(limits = c(-1, 28), breaks = seq(0, 25, by = 5)) +
  labs(x = TeX("$\\beta$"), y = TeX("$\\alpha$"), fill = "reward rate") +
  theme_bw() +
  guides(size = guide_legend(position = "inside"), fill = "none") +
  theme(legend.position.inside = c(0.15, 0.725)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))
B

# Panel C: Distribution of reward rates as a function of the size of the window.
rr_v_w$wsize <- as.factor(rr_v_w$wsize)
summaryStats <- rr_v_w %>%
  group_by(., wsize) %>%
  summarise(., meanRR = mean(objective),
            medianRR = median(objective),
            q25 = quantile(objective, 0.25),
            q75 = quantile(objective, 0.75))
C <- ggplot() + 
  geom_violin(data = rr_v_w, aes(x = wsize, y = objective), trim = FALSE, scale = "area", width = 1, fill = "grey90", colour = "black", alpha = 1) +
  geom_hline(yintercept = max(asymptoticRR$objective), linetype = "dashed") +
  geom_linerange(data = summaryStats, aes(x = wsize, ymin = q25, ymax = q75), linewidth = 1, colour = "black")+ # IQR
  geom_point(data = summaryStats, aes(x = wsize, y = medianRR), size = 3, shape = 15, colour= "grey50", position = position_dodge(width = 1), alpha = 1)+
  scale_x_discrete(breaks = c(5, 10, 25, 50), labels = c(5, 10, 25, 50)) +
  labs(x = 'window size', y='reward rate') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))
C

combinedPlot <- ggpubr::ggarrange(A, B, C,
                                  labels = c("A", "B", "C"), 
                                  ncol = 3,
                                  nrow = 1)
combinedPlot
ggsave("rewardRates.eps", dpi = 300, width = 21, height = 7, units = "cm")
  