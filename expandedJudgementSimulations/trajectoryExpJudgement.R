# C.J.H. Ludwig, June 2024
# Visualisation of trajectories in parameter space during a simulated expanded judgement task (to generate Figure 5 in the paper).
# This script calls a function that estimates the parameters for a given simulated agent (using a static model and a particle filter)
# and returns a plot of (i) the trajectories of three agents in a 2D parameter space, and (ii) the parameter estimates for each agent separately.

library(ggplot2)
library(dplyr)
library(viridis)
library(latex2exp)

rm(list=ls())

# Navigate to the appropriate folder and load the function that computes the log-likelihood for a given (stimulus, action) sequence.
setwd("")

source("decisionBoundEstimation.R")

# Load data generated from simulated agents.
# This dataset contains data and trajectories from three different types of agents:
#   - A stationary agent who does not move at all
#   - An agent who explores a little
#   - An agent who explores extensively
fName <- "data/simulatedAgents.RData"
load(fName)
selectAgents <- c(1, 2, 3)
agentID <- factor(c("static", "limited", "extensive"), levels = c("static", "limited", "extensive"), ordered = TRUE)

# For each participant, estimate parameters and generate individual visualisations
decisionBounds <- vector(mode = "list", length = length(selectAgents))
for (i in 1:length(selectAgents)){
  decisionBounds[[i]] <- decisionBoundEstimation(pID = agentID[i], traject = simulatedAgents[[i]]$visitData, allData = simulatedAgents[[i]]$allData)
}

# Combine plots into a multi-panel figure for an easy overview, but we don't actually save this.
combinedPlot <- ggpubr::ggarrange(decisionBounds[[1]]$searchPlot,
                                  decisionBounds[[2]]$searchPlot,
                                  decisionBounds[[3]]$searchPlot,
                                  decisionBounds[[1]]$parmPlot,
                                  decisionBounds[[2]]$parmPlot,
                                  decisionBounds[[3]]$parmPlot,
                                  labels = c("A", "B", "C"),
                                  ncol = 3,
                                  nrow = 2)
combinedPlot

### Re-generate visualisations so we can show the data in more compact form. ###

## Combined trajectory plot ##

# Extract the "visit" and "trajectory" data into dataframes
visitData <- bind_rows(decisionBounds[[1]]$uniqueLocs, decisionBounds[[2]]$uniqueLocs, decisionBounds[[3]]$uniqueLocs, .id = "agentID")
visitData$agentID <- agentID[as.numeric(visitData$agentID)]
trajectData <- bind_rows(decisionBounds[[1]]$trajectory, decisionBounds[[2]]$trajectory, decisionBounds[[3]]$trajectory, .id = "agentID")
trajectData$agentID <- factor(agentID[as.numeric(trajectData$agentID)], levels = agentID, ordered = TRUE)
meanLocs <- data.frame(agentID = c("static", "limited", "extensive"),
                       meanX = c(decisionBounds[[1]]$meanLocs[1],
                                 decisionBounds[[2]]$meanLocs[1],
                                 decisionBounds[[3]]$meanLocs[1]),
                       meanY = c(decisionBounds[[1]]$meanLocs[2],
                                 decisionBounds[[2]]$meanLocs[2],
                                 decisionBounds[[3]]$meanLocs[2])
                      )
meanLocsDyn <- filter(meanLocs, agentID != "static") # In case we only want to show the average position for the dynamic agents

# Some hard-coded parameters that we know are suitable for the agents we've selected to illustrate.
xyPeak <- c(-2.5*pi/180, 5) # Gradient and intercept with the "true" long-run maximum reward rate.
# Make sure the size and colour scales accommodate the range across agents
sizeScaleLims <- c(1,20)
sizeScaleBreaks <- c(1, 5, 10, 20)
colScaleLims <- c(-0.01, 0.041)
colScaleBreaks <- seq(-0.01, 0.04, by = 0.01)
# Set up dataframe with "start" and "end" annotations for dynamic agents. You
# will need to adjust the offsets with some trial-and-error. Values chosen here
# work for these particular agents and the image size we're saving (although even then, I ended up tweaking them offline).
txtOffsets <- data.frame(agentID = rep(c("limited", "extensive"), times = 2),
                         x = c(decisionBounds[[2]]$trajectory$gradient[1] + 0, 
                               decisionBounds[[3]]$trajectory$gradient[1] + 0,
                               decisionBounds[[2]]$trajectory$gradient[nrow(decisionBounds[[2]]$trajectory)] + 0, 
                               decisionBounds[[3]]$trajectory$gradient[nrow(decisionBounds[[3]]$trajectory)] + 0),
                         y = c(decisionBounds[[2]]$trajectory$intercept[1] - 0.8, 
                               decisionBounds[[3]]$trajectory$intercept[1] + 0.8, 
                               decisionBounds[[2]]$trajectory$intercept[nrow(decisionBounds[[2]]$trajectory)] + 1, 
                               decisionBounds[[3]]$trajectory$intercept[nrow(decisionBounds[[3]]$trajectory)] + 1),
                         txt = rep(c("start", "end"), each = 2)
                         )

staticSearchCombined <- ggplot() +
  geom_path(aes(x = gradient, y = intercept, group = agentID), data = trajectData, linewidth = 2, colour = "grey", alpha = 1) + # trajectory
  geom_point(aes(x = gradient, y = intercept, shape = agentID, colour = meanRewardRate, size = totalVisits), data = visitData) + # sampled points
  geom_point(aes(x = meanX, y = meanY, shape = agentID), data = meanLocs, colour = "black", size = 3) + # weighted average of the visited locations
  scale_shape(solid = TRUE) +
  geom_point(aes(x = xyPeak[1], y = xyPeak[2]), shape = 3, colour = "black", size = 3) + # peak of the objective
  # Play around with the position of the annotations by tweaking the offsets
  geom_text(aes(x = x, y = y, label = txt), data = txtOffsets) +
  scale_size_area(name = "# visits", limits = sizeScaleLims, breaks = sizeScaleBreaks) +
  scale_colour_viridis(na.value="transparent", limits = colScaleLims, breaks = colScaleBreaks) +
  scale_x_continuous(breaks = seq(-1, 0, by = 0.25), limits = c(-1.2, 0.2)) +
  scale_y_continuous(breaks = seq(0, 25, by = 5), limits = c(0, 26)) +
  labs(x = TeX("$\\beta$"), y = TeX("$\\alpha$"), colour = "reward rate", shape = "agent") +
  theme_bw() +
  # guides(size = guide_legend(position = "inside"), colour = guide_colorbar(position = "inside"), shape = guide_legend(position = "inside")) +
  # theme(legend.position.inside = c(0.85, 0.75)) +
  # theme(aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
  theme(axis.title.y.left = element_text(size=14,face="bold"))
staticSearchCombined
## Parameter estimates ##

# Extract the individual parameter estimates into dataframes.
pfEst <- bind_rows(decisionBounds[[1]]$pfParms, decisionBounds[[2]]$pfParms, decisionBounds[[3]]$pfParms, .id = "agentID")
pfEst$agentID <- agentID[as.numeric(pfEst$agentID)]
staticPostSamples <- bind_rows(decisionBounds[[1]]$staticSamples, decisionBounds[[2]]$staticSamples, decisionBounds[[3]]$staticSamples, .id = "agentID")
staticPostSamples$agentID <- agentID[as.numeric(staticPostSamples$agentID)]
staticParmSummary <- bind_rows(decisionBounds[[1]]$staticParms, decisionBounds[[2]]$staticParms, decisionBounds[[3]]$staticParms, .id = "agentID")
staticParmSummary$agentID <- agentID[as.numeric(staticParmSummary$agentID)]
# Add an offset on the trial axis to show the violins
staticPostSamples$trial <- staticPostSamples$trial + 15
staticParmSummary$trial <- staticParmSummary$trial + 15

# Given the facet_grid construction below, it's difficult to set the limits of the individual axes, so we do the following dirty trick.
dfBlank <- data.frame(agentID = rep(agentID, each = 6),
                      parameter = rep(factor(c("alpha", "beta", "eta")), times = 6),
                      x = rep(100, times = 18),
                      y = rep(c(0.1, -60*pi/180, 0.01, 25, 10*pi/180, 5), times = 3)
                      )

# Comparison between true and estimated parameters
estParmCombined <- ggplot(data = pfEst) +
  geom_ribbon(aes(x = trial, ymin = q_025, ymax = q_975), colour = "grey90", fill = "grey90", alpha = 1) + # 95% credible interval
  geom_line(aes(x = trial, y = trueVal), linewidth = 1, colour = "grey50", show.legend = TRUE) +
  geom_line(aes(x = trial, y = post_mean), colour = "black", linewidth = 1, alpha = 1, show.legend = TRUE) + # Posterior mean
  # To the right, show the estimates from a static model (posterior parameters from Stan)
  geom_violin(data = staticPostSamples, aes(x = trial, y = sampleVal), trim = FALSE, scale = "width", fill = "grey90", colour = "black", alpha = 1)+
  geom_errorbar(data = staticParmSummary, aes(x = trial, ymin = q_25, ymax = q_75), size = 1, width = 3, colour = "black")+ # IQR
  geom_point(data = staticParmSummary, aes(x = trial, y = mean), size = 2, shape = 15, colour= "grey50", position = position_dodge(width = 1), alpha = 1)+
  geom_blank(data = dfBlank, aes(x = x, y = y)) +
  facet_grid(parameter~agentID, scales = "free") +
  # scale_x_continuous(breaks = c(0, 50, 100, 150, 200), labels = c(0, 50, 100, 150, 200)) +
  labs(x = 'trial', y='parameter value') +
  theme_bw()+
  theme(strip.text=element_text(size=14,face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.title.x.bottom = element_text(size=14,face="bold"))+
  theme(axis.title.y.left = element_text(size=14,face="bold"))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=14, face="bold"))
estParmCombined

## RMSEs ##
dfRMSE <- bind_rows(decisionBounds[[1]]$RMSE, decisionBounds[[2]]$RMSE, decisionBounds[[3]]$RMSE, .id = "agentID")
dfRMSE$agentID <- agentID[as.numeric(dfRMSE$agentID)]
dfRMSE <- select(dfRMSE, -c(rmse_StaticPostMedian, rmse_pfPostMedian))
dfRMSE <- tidyr::pivot_longer(dfRMSE, cols = starts_with("rmse_"), 
                           values_to = "RMSE",
                           names_to = "comparison",
                           names_prefix = "rmse_")
dfRMSE$comparison <- factor(dfRMSE$comparison, levels = c("StaticPostMean", "pfPostMean"),
                            labels = c("static", "dynamic"),
                            ordered = TRUE)
dfRMSE$lnRMSE <- log(dfRMSE$RMSE)
yScaleLims <- c(.002, 5.1) 
yScaleBreaks <- c(0.01, 0.1, 1, 5)

rmseCombined <- ggplot(data = dfRMSE) +
  geom_point(aes(x = parameter, y = RMSE, shape = agentID, colour = comparison, group = agentID), position = position_dodge(width = 0.5), size = 4) +
  scale_colour_grey() +
  scale_y_log10(limits = yScaleLims, breaks = yScaleBreaks) +
  labs(x = 'parameter', y='RMSE', colour = "", shape = "") +
  guides(shape = "none") +
  # theme(legend.position.inside = c(1, 1)) +
  theme_bw()+
  # theme(aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text=element_text(size=14,face="bold"))+
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.title.x.bottom = element_text(size=14,face="bold"))+
  theme(axis.title.y.left = element_text(size=14,face="bold"))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=14, face="bold"))
rmseCombined

# Combine plots into a multi-panel figure.
# This is really ugly and hacky, but does the job: we combine the top row first and then put the parameter estimates underneath.
topRow <- ggpubr::ggarrange(staticSearchCombined, rmseCombined,
                            labels = c("A", "C"),
                            ncol = 2,
                            nrow = 1)
combinedPlot <- ggpubr::ggarrange(topRow, estParmCombined,
                                  labels = c("", "B"),
                                  ncol = 1,
                                  nrow = 2)
combinedPlot

# With this combined plot, the only things to tweak offline are:
#   - 'start' and 'end' annotations in panel A (though we could do this in txtOffsets).
#   - Make the weighted average locations in panel A hollow.
#   - Legend positioning (I actually save a larger version of the plot and grab the legends from there, then paste into the smaller figure).
#   - Parameter names as Greek symbols (axis labels in panel B; row labels in panel C).
#   - Remove the right most x-axis tick mark and numerical label
# The following size just about fits on A4 when rotated 90 degs. You need to do the legend positioning offline though.
ggsave("expJudgeEstimates.eps", dpi = 300, width = 26, height = 18.5, units = "cm")
