# C.J.H. Ludwig, March - June 2024
# Simple demo of an agent varying two parameters in order to maximise an objective.
# In this demo, the objective is arbitrary and the agent can only get uncertain estimates of the objective for any combination of parameters.
# For the sake of illustration, the parameters that are varied are simply the mean and standard deviation of a Gaussian distribution.
# At each location in the parameter space, the agent generates observations: random draws from a Gaussian distribution.
# We then estimate the parameters of a Gaussian model fit to the agent's data.
# We either estimate these with a static model in which the parameters remain constant over time, or we use a simple bootstrap particle filter that estimates time-varying parameters sequentially.
# These simulations and model fits are used to generate Figure 4 in the paper.

library(tidyverse)
library(viridis)
library(latex2exp)
library(gganimate)
library(rstan)

rm(list=ls())
dev.off()

# Load a script that contains the expression of a 2D objective, along with some function(s) to search this objective.
source("objectiveFunctions.R")
source("samplingFunctions.R")
source("timeVarNormal_bootstrapPF.R")

set.seed(1042024)

###Set up the 2D parameter space and compute objective function###
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

thetaGrid$objective <- obj2D_uniModal(x=thetaGrid$theta1, y=thetaGrid$theta2, mu=muLoc, sigma=width, corrxy=rho, amp=peak)
# Compute the location of the peak - because of the non-linear correlation between the two parameters, this is not necessarily 'muLoc'.
peakLoc <- optim(muLoc, obj2D_uniModal_sqErr, mu = muLoc, sigma = width, corrxy = rho, amp = peak)
xyPeak <- peakLoc$par

###Let an agent explore this parameter space###
totalN <- 100
sigma0 <- 25 # Noise around the objective estimates
windowSize <- 5 # Window over which the agent estimates the objective
stepSDs <- c(1,1) # Size of steps in 2D parameter space
satisficeThreshold <- 1 # Stopping criterion
# Given these parameters, run local sampling and comparison algorithm.
simulateTrajectory <- LoSaCo1(sd0 = sigma0,
                              stepParms = stepSDs,
                              nMaxSteps = totalN,
                              xrange = c(min(theta1), max(theta1)),
                              yrange = c(min(theta2), max(theta2)),
                              N = totalN, 
                              W = windowSize, 
                              satisfice = satisficeThreshold,
                              muLoc = muLoc,
                              width = width,
                              rho = rho,
                              peak = peak)
simulateTrajectory$sample_window <- seq(1:nrow(simulateTrajectory))
# Compute weighted average location
wMean <- c(weighted.mean(simulateTrajectory$theta1, simulateTrajectory$nSamples),
           weighted.mean(simulateTrajectory$theta2, simulateTrajectory$nSamples))

# Visualise the objective and the search trajectory
# You need to adjust the number of points shown in the size scaling legend.
theta1Lims <- c(min(c(-3, min(simulateTrajectory$theta1) - 0.1)), max(c(3, max(simulateTrajectory$theta1) + 0.1)))
theta2Lims <- c(min(c(0, min(simulateTrajectory$theta2) - 0.1)), max(c(5, max(simulateTrajectory$theta2) + 0.1)))
sizeScaleLims <- c(1, max(c(100,max(simulateTrajectory$cumulativeSamples))))
colScaleLims <- c(min(c(0,min(simulateTrajectory$estObjective))), 
                  max(c(100, max(simulateTrajectory$estObjective))))

###Generate "observable data" from this trajectory.###
obsData <- simulateTrajectory %>%
  slice(rep(1:nrow(simulateTrajectory), nSamples)) %>%
  select(theta1, theta2, meanObjective, sample_window)
obsData$trialID <- 1:nrow(obsData)
# Now for each trial, generate a draw from a Gaussian with mu = theta1, sd = theta2
for (i in 1:nrow(obsData)){
  obsData$y[i] <- rnorm(1, mean = obsData$theta1[i], sd = obsData$theta2[i])
}

###Estimate parameters ignoring the temporal dynamics.###
# This is over the top, but we use Stan to sample from the marginal posteriors for mu and sigma, with priors that are the same as we use in the particle filter below.
stanData <- list(N = totalN, y = obsData$y)
stanParms <- stan("postMeanSD.stan", data = stanData)
print(fitSummary <- data.frame(summary(stanParms, pars = c("mu", "sigma"))$summary))
dfPostSamples <- as.data.frame(stanParms, pars = c("mu","sigma")) # Extract samples in a dataframe
# Manipulate the posterior samples and summary statistics to get in a useful format for plotting
# For fitSummary we need to do some renaming to facilitate subsequent plotting
fitSummary <- select(fitSummary, c("mean", "X2.5.", "X25.", "X50.", "X75.", "X97.5."))
colnames(fitSummary) <- c("mean", "q_025", "q_25", "q_50", "q_75", "q_975")
fitSummary$parameter <- rownames(fitSummary)
fitSummary$trialID <- max(obsData$trialID) + 10 # add a "virtual" trial ID
dfPostSamples$trialID <- max(obsData$trialID) + 10 # add a "virtual" trial ID
dfPostSamples <- tidyr::pivot_longer(dfPostSamples,
                                   cols = c("mu", "sigma"), 
                                   names_to = "parameter", 
                                   values_to = "sampleVal")

### Estimate the trajectory based on the observable data ###

# Define the log-likelihood for a single observation
logLikOneTrial <- function(params, obs, returnLogLik = FALSE){
  dnorm(obs, mean = params[1], sd = params[2], log = returnLogLik)
}

# Set some of the hyperparameters of the particle filter algorithm: parameters of the truncated multivariate Gaussian prior and transition distributions. We specify the prior means (mu0), the prior variances (sigma0) and the variances of the transition distribution (sigmaQ).
hParms <- list(mu0 = c(0, 2.5), # prior means 
               sigma0 = matrix(c(1^2, 0, 0, 1^2), ncol = 2), # 2 x 2 covariance matrix for initial particle values; note sds on the diagonal are squared to give variance.
               sigmaQ = matrix(c(0.5^2, 0, 0, 0.5^2), ncol = 2), # 2 x 2 covariance matrix for transition distribution; note sds on the diagonal are squared to give variance.
               lims = matrix(c(-3, 3, 1e-6, 5), ncol = 2) # 2 x 2 matrix with lower and upper bounds in the two rows.
)
nPart <- 2000

# Bootstrap particle filter parameter estimation for dataset
pfOut <- pfRun(data = obsData, n = nPart, hyperParms = hParms)

# Match the true data generating parameters for each trial with the inferred states.
pfSummary <- summaryStats(pfOut) # Get useful summary statistics for each trial
parmDF <- obsData %>%
  rename(., mu = theta1, sigma = theta2) %>% # Rename generating parameters
  tidyr::pivot_longer(., cols = c("mu", "sigma"), 
                      names_to = "parameter", 
                      values_to = "trueVal") %>% # From wide to long
  left_join(., pfSummary, by=c("trialID", "parameter")) # Join the two dataframes, based on trialID and parameter columns

# Quantify the goodness-of-fit: RMSE
parmDF$staticPostMean <- rep(fitSummary$mean, times = stanData$N)
parmDF$staticPostMedian <- rep(fitSummary$q_50, times = stanData$N)
parmDF$staticSquareErrMean <- sqrt((parmDF$staticPostMean - parmDF$trueVal)^2)
parmDF$staticSquareErrMedian <- sqrt((parmDF$staticPostMedian - parmDF$trueVal)^2)
parmDF$pfSquareErrMean <- sqrt((parmDF$post_mean - parmDF$trueVal)^2)
parmDF$pfSquareErrMedian <- sqrt((parmDF$q_5 - parmDF$trueVal)^2)
errorSummary <- parmDF %>%
  group_by(., parameter) %>%
  summarise(., rmse_StaticPostMean = mean(staticSquareErrMean, na.rm = TRUE),
            rmse_StaticPostMedian = mean(staticSquareErrMedian, na.rm = TRUE),
            rmse_pfPostMean = mean(pfSquareErrMean, na.rm = TRUE), 
            rmse_pfPostMedian = mean(pfSquareErrMedian, na.rm = TRUE))
print(errorSummary)

### Visualisations
# May need to tweak the positioning of certain non-key elements, such as text annotations and legends.
# The strip text in panel D is changed into Greek symbols outside of this script (could not get this changed in R).

# If you prefer: use a colour-blind palette.
# cbbPalette <- c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour-deficient-friendly palette; first four entries: black, orange, blue, green

# Objective
objPlot <- ggplot() +
  geom_raster(aes(x = theta1, y = theta2, fill = objective), data = thetaGrid, interpolate = TRUE) +
  scale_fill_viridis(na.value="transparent", limits = colScaleLims, breaks=c(0,25,50,75,100)) +
  geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyPeak[1], theta2 = xyPeak[2]), shape = 3, colour = "black", size = 3) +
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5)) +
  labs(x = TeX("$\\psi_1$"), y = TeX("$\\psi_2$")) +
  theme_bw() +
  guides(fill = guide_colourbar(position = "inside")) +
  theme(legend.position.inside = c(0.17, 0.625)) +
  # theme(aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))
objPlot

# Search trajectory through the parameter space
staticSearchPlot <- ggplot() +
  geom_path(aes(x = theta1, y = theta2), data = simulateTrajectory, linewidth = 2, colour = "grey", alpha = 1) + # trajectory
  geom_point(aes(x = theta1, y = theta2, colour = meanObjective, size = cumulativeSamples), data = simulateTrajectory) + # sampled points
  geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyPeak[1], theta2 = xyPeak[2]), shape = 3, colour = "black", size = 3) + # peak of the objective
  geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = wMean[1], theta2 = wMean[2]), shape = 17, colour = "black", size = 3) + # weighted average of the visited locations
  # Play around with the position of the annotations by tweaking the offsets
  annotate("text", x = simulateTrajectory$theta1[1] + 0.1, y = simulateTrajectory$theta2[1] + 0.25, label = "start") +
  annotate("text", x = simulateTrajectory$theta1[nrow(simulateTrajectory)] + 0, y = simulateTrajectory$theta2[nrow(simulateTrajectory)] + 0.25, label = "end") +
  scale_size_area(name = "samples", limits = sizeScaleLims, breaks=c(5, 25, 50)) +
  scale_colour_viridis(na.value="transparent", limits = colScaleLims, breaks=c(0,25,50,75,100)) +
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3), limits=theta1Lims) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits=theta2Lims) +
  labs(x = TeX("$\\mu$"), y = TeX("$\\sigma$")) +
  theme_bw() +
  guides(size = guide_legend(position = "inside"), colour = "none") +
  theme(legend.position.inside = c(0.125, 0.75)) +
  # theme(aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold")) +
  # Show the mean posterior estimates in the plot of the search trajectory
  geom_point(aes(x = mean[1], y = mean[2]), data = fitSummary, shape = 2, colour = "black", size = 3) + # weighted average of the visited locations
  geom_errorbar(aes(x = mean[1], ymin = q_025[2], ymax = q_975[2]), data = fitSummary, colour = "black", size = 0.5, width = 0.1) +
  geom_errorbarh(aes(y = mean[2], xmin = q_025[1], xmax = q_975[1]), data = fitSummary, colour = "black", size = 0.5, height = 0.1)
staticSearchPlot

# Time series of data and the marginal distribution (mixture of Gaussians)
dataPlot <- ggplot(data = obsData) +
  geom_point(aes(x = trialID, y = y), size = 2, shape = 16, colour = "black") +
  geom_violin(aes(x = max(trialID) + 20, y = y), trim = FALSE, scale = "area", width = 15, fill = "grey90", colour = "black", alpha = 1) +
  labs(x = "trial", y = "y") +
  scale_x_continuous(breaks = seq(0, 100, by = 25)) +
  scale_y_continuous(limits = c(-10, 10), breaks = c(-10, -5, 0, 5, 10)) +
  theme_bw() +
  # theme(aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))
dataPlot

# Parameter estimates from particle filter and the static model fit
parmPlot <- ggplot(data = parmDF) +
  geom_ribbon(aes(x = trialID, ymin = q_025, ymax = q_975), colour = "grey90", fill = "grey90", alpha = 1) + # 95% credible interval
  geom_line(aes(x = trialID, y = trueVal), size = 1, colour = "grey50", show.legend = TRUE) +
  geom_line(aes(x = trialID, y = post_mean), colour = "black", alpha = 1, size = 1, show.legend = TRUE) + # Posterior mean
  # To the right, show the estimates from a static model (posterior parameters from Stan)
  geom_violin(data = dfPostSamples, aes(x = trialID, y = sampleVal), trim = FALSE, scale = "area", width = 15, fill = "grey90", colour = "black", alpha = 1)+
  geom_errorbar(data = fitSummary, aes(x = trialID, ymin = q_25, ymax = q_75), size = 1, colour = "black")+ # IQR
  geom_point(data = fitSummary, aes(x = trialID, y = mean), size = 2, shape = 15, colour= "grey50", position = position_dodge(width = 1), alpha = 1)+
  facet_grid(parameter~., scales = "free_y", labeller = as_labeller(c(`mu` = TeX(r'($\mu$)'), `sigma` = TeX(r'($\sigma$)')))) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, 50, 75, 100)) +
  labs(x = 'trial', y='parameter value') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text=element_text(size=14,face="bold"))+
  theme(axis.text.x  = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.title.x.bottom = element_text(size=14,face="bold"))+
  theme(axis.title.y.left = element_text(size=14,face="bold"))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=14, face="bold"))
parmPlot

# Combine plots into a multi-panel figure.
combinedPlot <- ggpubr::ggarrange(objPlot, staticSearchPlot, dataPlot, parmPlot,
                  labels = c("A", "B", "C", "D"), 
                  ncol = 2,
                  nrow = 2)
combinedPlot
ggsave("exploreObj2D_Gauss.eps", dpi = 300, width = 21, height = 16, units = "cm")

# Add estimated trajectory (posterior mean) to the search plot. 
# This looks ugly, even after some filtering (median filter with same window size as the data generating agent)
# filtered_mu <- robfilter::rm.filter(parmDF$post_mean[parmDF$parameter == "mu"], width = windowSize, online = FALSE)
# filtered_sigma <- robfilter::rm.filter(parmDF$post_mean[parmDF$parameter == "sigma"], width = windowSize, online = FALSE)
# filteredParms <- data.frame(trialID = seq(1,totalN),
#                             mu = filtered_mu$level,
#                             sigma = filtered_sigma$level)
# colnames(filteredParms) <- c("trialID", "mu", "sigma")
# 
# staticSearchPlot <- staticSearchPlot +
#   geom_line(aes(x = mu, y = sigma),
#             data = filteredParms,
#             colour = "grey",
#             size = 0.5,
#             alpha = 0.25) +
#   annotate("text", x = filteredParms$mu[1], y = filteredParms$sigma[1], label = "pf1") +
#   annotate("text", x = filteredParms$mu[totalN], y = filteredParms$sigma[totalN], label = "pf100")
# staticSearchPlot


