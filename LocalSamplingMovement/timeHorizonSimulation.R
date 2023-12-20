# C.J.H. Ludwig, July 2023.
# Function that runs the time horizon simulations, combining the different arguments (e.g. different noise levels and window sizes).

timeHorizonSimulation <- function(simulationN = NULL, 
                                  totalN = NULL, 
                                  sigma0 = NULL,
                                  windowSize = NULL,
                                  satisficeThreshold = NULL){
  
  # Set up the 2D parameter space and compute objective function. These settings are the same as in Ludwig et al.
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
  
  # Visualise the objective.
  plotObj <- ggplot() +
    geom_raster(aes(x = theta1, y = theta2, fill = objective), data = thetaGrid, interpolate = TRUE) +
    scale_fill_viridis(na.value="transparent", limits = c(0,100), breaks=c(0,25,50,75,100)) +
    geom_point(aes(x = theta1, y = theta2), data = data.frame(theta1 = xyPeak[1], theta2 = xyPeak[2]), shape = 3, colour = "black", size = 3) +
    # geom_point(aes(x = theta1, y = theta2), data = xyGrid, size = 2, colour = "grey", alpha = 0.5) +
    scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3)) +
    scale_y_continuous(breaks = c(0,1,2,3,4,5)) +
    labs(x = TeX("$\\theta_1$"), y = TeX("$\\theta_2$")) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.x  = element_text(size=12)) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
    theme(axis.title.y.left = element_text(size=14,face="bold"))
  
  
  # Now run a number of agents, defined by: different noise levels; different sampling windows; satisficing vs. maximising; the time horizon (totalN). The last variable is key here. All the other variables are a (judiciously chosen) sub-set of the values reported in Ludwig et al.
  ## Start by setting up a dataframe that holds the simulation parameters and results.
  simulationRun <- seq(1,simulationN)
  simulationResults <- expand.grid(windowSize, sigma0, satisficeThreshold, totalN,  simulationRun)
  colnames(simulationResults) <- c("windowSize", "sd0", "satisficeThreshold", "N", "run")
  simulationResults$stopRule[simulationResults$satisficeThreshold == 100] <- "maximising"
  simulationResults$stopRule[simulationResults$satisficeThreshold < 100] <- "satisficing"
  simulationResults$stepSD1 <- 1
  simulationResults$stepSD2 <- 1
  simulationResults$maxSteps <- simulationResults$N # Allow for W = 1, i.e. maximum number of steps equals the number of trials
  simulationResults$totalObjective <- 0
  simulationResults$nUniqueLocs <- 0
  simulationResults$spatialSpread <- 0
  
  ## Now just loop through the dataframe and run a simulation for each row. I know the for-loop is inefficient, but it's clear and time is not of the essence here!
  for(i in 1:nrow(simulationResults)){
    cat(i, "out of", nrow(simulationResults), "\n")
    tmpResult <- LoSaCo1(sd0 = simulationResults$sd0[i], 
                         stepParms = c(simulationResults$stepSD1[i],
                         simulationResults$stepSD2[i]),
                         nMaxSteps = simulationResults$maxSteps[i],
                         xrange = c(min(theta1), max(theta1)),
                         yrange = c(min(theta2), max(theta2)),
                         N = simulationResults$N[i], 
                         W = simulationResults$windowSize[i], 
                         satisfice = simulationResults$satisficeThreshold[i],
                         muLoc = muLoc,
                         width = width,
                         rho = rho,
                         peak = peak)
    # How many unique locations visited?
    nUniqueLocs <- sum(tmpResult$cumulativeSamples == tmpResult$nSamples)
    # Compute the weighted average (x,y) location over the whole trajectory
    wAvLoc <- c(sum(tmpResult$theta1*tmpResult$nSamples)/simulationResults$N[i], sum(tmpResult$theta2*tmpResult$nSamples)/simulationResults$N[i])
    # What is the spatial spread of the search, relative to this mean?
    distanceFromAvLoc <- sqrt((tmpResult$theta1 - wAvLoc[1])^2 + (tmpResult$theta2 - wAvLoc[2])^2) # First compute the distances from the weighted average location
    distanceSD <- sum(distanceFromAvLoc*tmpResult$nSamples)/(simulationResults$N[i]*(nrow(tmpResult)-1)/nrow(tmpResult))
    
    
    simulationResults$totalObjective[i] <- tmpResult$cumulativeObjective[nrow(tmpResult)]
    simulationResults$nUniqueLocs[i] <- nUniqueLocs
    simulationResults$spatialSpread[i] <- distanceSD
  }

  # Summarise the simulation results. We have three key measures: the mean gain, the number of unique locations visited and the spatial spread of the trajectory. For each of these, compute the mean, median and CI.
  simulationSummary <- simulationResults %>%
    group_by(., sd0, windowSize, stopRule, N) %>%
    summarise(., meanTotalObjective = mean(totalObjective), 
              medianTotalObjective = median(totalObjective), 
              objectiveCI = qt(0.975, df = simulationN-1) * sd(totalObjective)/sqrt(simulationN),
              meanUniqueLocs = mean(nUniqueLocs),
              medianUniqueLocs = median(nUniqueLocs),
              nLocsCI = qt(0.975, df = simulationN-1) * sd(nUniqueLocs)/sqrt(simulationN),
              meanSpread = mean(spatialSpread),
              medianSpread = median(spatialSpread),
              spreadCI = qt(0.975, df = simulationN-1) * sd(spatialSpread)/sqrt(simulationN)
              )
  
  # Show results in a way to highlight the comparison between maximising and satisficing
  # Plot the mean gain
  cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  plotGain <- ggplot(data = simulationSummary, aes(x = windowSize, y = meanTotalObjective, colour = stopRule)) +
    geom_point(size = 4, pch = 16, position=position_dodge(width=0.5)) +
    geom_line(size = 2) +
    geom_errorbar(aes(ymin = meanTotalObjective - objectiveCI, ymax = meanTotalObjective + objectiveCI), position=position_dodge(width=0.5)) +
    facet_grid(as.factor(sd0)~as.factor(N)) +
    scale_x_continuous(limits = c(0,51), breaks = c(0,5,10,25,50)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0,25,50,75,100), oob = scales::oob_keep) +
    scale_color_manual(values = cbbPalette[1:length(satisficeThreshold)]) +
    labs(x = TeX("window size"), y = TeX("gain per trial"), color = "stopping rule") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(panel.grid.minor = element_blank()) +
    theme(strip.text = element_text(size = 14, face = "bold")) +
    theme(axis.text.x  = element_text(size=12)) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
    theme(axis.title.y.left = element_text(size=14,face="bold"))
  
  # Plot the mean number of unique locations visited
  plotLocs <- ggplot(data = simulationSummary, aes(x = windowSize, y = meanUniqueLocs, colour = stopRule)) +
    geom_point(size = 4, pch = 16, position=position_dodge(width=0.5)) +
    geom_line(size = 2) +
    geom_errorbar(aes(ymin = meanUniqueLocs - nLocsCI, ymax = meanUniqueLocs + nLocsCI), position=position_dodge(width=0.5)) +
    facet_grid(as.factor(sd0)~as.factor(N)) +
    scale_x_continuous(limits = c(0,51), breaks = c(0,5,10,25,50)) +
    scale_y_continuous(limits = c(1,2000), breaks = c(1, 10, 100, 1000), trans = "log", oob = scales::oob_keep) +
    scale_color_manual(values = cbbPalette[1:length(satisficeThreshold)]) +
    labs(x = TeX("window size"), y = TeX("# unique locations"), color = "stopping rule") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(panel.grid.minor = element_blank()) +
    theme(strip.text = element_text(size = 14, face = "bold")) +
    theme(axis.text.x  = element_text(size=12)) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
    theme(axis.title.y.left = element_text(size=14,face="bold"))
  
  # Finally, plot the mean spatial spread
  plotSpread <- ggplot(data = simulationSummary, aes(x = windowSize, y = meanSpread, colour = stopRule)) +
    geom_point(size = 4, pch = 16, position=position_dodge(width=0.5)) +
    geom_line(size = 2) +
    geom_errorbar(aes(ymin = meanSpread - spreadCI, ymax = meanSpread + spreadCI), position=position_dodge(width=0.5)) +
    facet_grid(as.factor(sd0)~as.factor(N)) +
    scale_x_continuous(limits = c(0,51), breaks = c(0,5,10,25,50)) +
    scale_y_continuous(limits = c(0,20), breaks = c(0,5,10,15,20), oob = scales::oob_keep) +
    scale_color_manual(values = cbbPalette[1:length(satisficeThreshold)]) +
    labs(x = TeX("window size"), y = TeX("spatial spread"), color = "stopping rule") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(panel.grid.minor = element_blank()) +
    theme(strip.text = element_text(size = 14, face = "bold")) +
    theme(axis.text.x  = element_text(size=12)) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
    theme(axis.title.y.left = element_text(size=14,face="bold"))
  
  return(list(simulationSummary = simulationSummary, simulationResults = simulationResults, pObj = plotObj, pGain = plotGain, pLocs = plotLocs, pSpread = plotSpread))
}