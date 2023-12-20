# C.J.H. Ludwig, August-September 2022

# Script that illustrates some crude optimisation procedures that a hypothetical agent might use in navigating a 2D objective function.

library(tidyverse)
library(viridis)
library(latex2exp)
library(gganimate)

rm(list=ls())

# Load a script that contains the expression of a 2D objective, along with some function(s) to search this objective.
source("objectiveSearch.R")

###Define and plot a 2D objective function###

# Set up a grid of our two parameters
range1 <- c(-3,3)
range2 <- c(1e-6,5)
theta1 <- seq(range1[1],range1[2],length.out=20)
theta2 <- seq(range2[1],range2[2],length.out=20)
thetaGrid <- expand.grid(theta1, theta2)
colnames(thetaGrid) <- c("theta1", "theta2")
muLoc <- c(0.5,0.5)
width <- c(1.5,1)
rho <- 0.75
peak <- 100

# Compute the objective
thetaGrid$objective <- obj2D(x=thetaGrid$theta1, y=thetaGrid$theta2, mu=muLoc, sigma=width, corrxy=rho, amp=peak)
# Compute the location of the peak - because of the non-linear correlation between the two parameters, this is not necessarily 'muLoc'.
peakLoc <- optim(muLoc, obj2Derror, mu = muLoc, sigma = width, corrxy = rho, amp = peak)
xyPeak <- peakLoc$par

# Plot the objective as a heatmap
ggplot() +
  geom_raster(aes(x = theta1, y = theta2, fill = objective), data = thetaGrid, interpolate = TRUE) +
  scale_fill_viridis(na.value="transparent", limits = c(0,100), breaks=c(0,25,50,75,100)) +
  geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyPeak[1], theta2 = xyPeak[2]), shape = 3, colour = "black", size = 3) +
  # geom_point(aes(x = theta1, y = theta2), data = xyGrid, size = 2, colour = "grey", alpha = 0.5) +
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5)) +
  labs(x = TeX("$\\theta_g^1$"), y = TeX("$\\theta_g^2$")) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))

###Simulate an agent who navigates the objective using a simple grid search###
sigma_0 <- 25
totalN <- 500
windowSize <- 5
nParmLevels <- c(4,3) # number of levels over which the two parameters can be varied
gridSearchResult <- gridSearch(sd0 = sigma_0, steps <- nParmLevels, 
                    xrange = c(min(theta1), max(theta1)),
                    yrange = c(min(theta2), max(theta2)),
                    N = totalN, W = windowSize)
searchGrid <- filter(gridSearchResult$grid, visited == TRUE)
searchPath <- data.frame(theta1 = gridSearchResult$grid$theta1[gridSearchResult$locs],
                         theta2 = gridSearchResult$grid$theta2[gridSearchResult$locs])
searchText <- data.frame(theta1 = searchPath$theta1[c(1,nrow(searchPath))],
                         theta2 = searchPath$theta2[c(1,nrow(searchPath))],
                         annoText = c("start", "end"))

# Visualise the search
ggplot() +
  geom_path(aes(x = theta1, y = theta2), data = searchPath, size = 2, colour = "grey", alpha = 0.5) +
  geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyPeak[1], theta2 = xyPeak[2]), shape = 3, colour = "black", size = 3) +
  geom_point(aes(x = theta1, y = theta2, colour = meanObjective, size = nSamples), data = searchGrid) +
  geom_text(aes(x = theta1, y = theta2, label = annoText), data = searchText, nudge_x = c(0.25,0), nudge_y = c(0.25, -0.25)) +
  scale_size_area(name = "samples", limits = c(0,250)) +
  scale_colour_viridis(na.value="transparent", limits = c(0,100), breaks=c(0,25,50,75,100)) +
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5)) +
  labs(x = TeX("$\\psi_g^1$"), y = TeX("$\\psi_g^2$"), colour = "objective value", size = "samples") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))

###Simulate an agent who can only take steps of a certain size, but is able to use crude gradient information.###
sigma_0 <- 50
totalN <- 500
stepSD <- c(1,1) # Parameters of the multivariate normal that controls the steps taken in the parameter space
nSteps <- 500 # Maximum number of steps we allow
windowSize <- 10
satisficeThreshold <- obj2D(x=xyPeak[1], y=xyPeak[2], mu=muLoc, sigma=width, corrxy=rho, amp=peak) * 0.75 # Stopping criterion for the search
greedySearchResult <- greedySearch(sd0 = sigma_0, stepParms = stepSD,
                               nMaxSteps = nSteps,
                               xrange = c(min(theta1), max(theta1)),
                               yrange = c(min(theta2), max(theta2)),
                               N = totalN, W = windowSize, 
                               satisfice = satisficeThreshold)
greedySearchResult$sample_window <- seq(1:nrow(greedySearchResult))

searchText <- data.frame(theta1 = greedySearchResult$theta1[c(1,nrow(greedySearchResult))],
                         theta2 = greedySearchResult$theta2[c(1,nrow(greedySearchResult))],
                         annoText = c("start", "end"))
theta1Lims <- c(min(c(-3, min(greedySearchResult$theta1) - 0.1)), max(c(3, max(greedySearchResult$theta1) + 0.1)))
theta2Lims <- c(min(c(0, min(greedySearchResult$theta2) - 0.1)), max(c(5, max(greedySearchResult$theta2) + 0.1)))
sizeScaleLims <- c(1, max(c(100,max(greedySearchResult$cumulativeSamples))))
colScaleLims <- c(min(c(0,min(greedySearchResult$estObjective))), 
                  max(c(100, max(greedySearchResult$estObjective))))

# Visualise the search (static)
# You need to adjust the number of points shown in the size scaling legend.
# You also need to adjust the position of the 'start' and 'end' annotations (in geom_text). Indeed, comment this line out for generating the animation further down.
staticSearchPlot <- ggplot() +
  geom_path(aes(x = theta1, y = theta2), data = greedySearchResult, size = 2, colour = "grey", alpha = 1) +
  geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyPeak[1], theta2 = xyPeak[2]), shape = 3, colour = "black", size = 3) +
  geom_text(aes(x = theta1, y = theta2, label = annoText), data = searchText, nudge_x = c(0.3,-0.35), nudge_y = c(0,0)) +
  geom_point(aes(x = theta1, y = theta2, colour = estObjective, size = cumulativeSamples), data = greedySearchResult) +
  scale_size_area(name = "samples", limits = sizeScaleLims, breaks=c(5, 50, 100, 200, 275)) +
  scale_colour_viridis(na.value="transparent", limits = colScaleLims, breaks=c(0,25,50,75,100)) +
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3), limits=theta1Lims) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits=theta2Lims) +
  labs(x = TeX("$\\psi_1$"), y = TeX("$\\psi_2$"), color = "estimate") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))
staticSearchPlot

ggsave("staticSearch_window10_noise50_satisfice75.eps", plot = staticSearchPlot, height = 97, units = "mm", dpi = 300)
#ggsave("staticSearch_window5_noise25_2.eps", plot = staticSearchPlot, height = 89, units = "mm", dpi = 300)

# Visualise the search with an animation
dynSearchPlot <- staticSearchPlot +
  transition_reveal(sample_window)
animGreedySearch <- animate(dynSearchPlot, fps=5)
anim_save("animGreedySearch_window10_noise50_satisfice75.gif")


