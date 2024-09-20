# C.J.H. Ludwig, June 2024
# Function to estimate and visualise trajectories in parameter space during a
# simulated expanded judgement paradigm (cf. Malhotra et al., 2017).
# Returns a list with 
#   - Static parameter estimates
#   - Particle filter parameter estimates
#   - RMSE from static and particle filter estimates
#   - 2D trajectory plot (intercept, gradient)
#   - 3 x 1 plot of true vs estimated parameters

library(ggplot2)
library(dplyr)
library(viridis)
library(latex2exp)

decisionBoundEstimation <- function(pID = NULL, traject = NULL, allData = NULL){
  
  ### Some initial data manipulation and checks ###
  allData$participant <- pID
  allData$noise <- 1 # Add column with the (fixed) noise parameter that was used to generate the data
  
  # Get the accuracy for each trial and get rid of any time-outs 
  # (however, we don't want datasets with time-outs, because that messes up the particle filter estimates)
  trials <- unique(allData$trial)
  nTrials <- length(trials)
  accuracy <- rep(NA, nTrials)
  rejectTrials <- 0
  for(i in 1:nTrials){
    trialData <- filter(allData, trial == trials[i]) # select evidence and action sequences for the current trial
    finalAction <- trialData$action[nrow(trialData)]
    if(finalAction == 1){
      allData$keepTrial[allData$trial == trials[i]] <- 1 # Keep the trial if the final action was 'go'
      if (trialData$evi[trialData$action==1]>0){
        accuracy[i] <- 1
      } else {
        accuracy[i] <- 0
      }
    } else { # If the final action was not 'go' just get rid of the trial to avoid problems down the line
      allData$keepTrial[allData$trial == trials[i]] <- 0
      rejectTrials <- rejectTrials + 1
    }
  }
  cat(rejectTrials, "trials rejected\n") # Show if we've had to reject any trials
  allData <- filter(allData, keepTrial==1)
  nTrials <- nTrials - rejectTrials
  accuracy <- accuracy[!is.na(accuracy)]
  
  ### Particle filter to identify parameter dynamics ###
  
  # Set some of the hyperparameters of the algorithm: parameters of the
  # truncated multivariate Gaussian prior and transition distributions. We
  # specify the prior means (mu0), the prior variances (sigma0) and the
  # variances of the transition distribution (sigmaQ). Note that we do not
  # expect the decision noise to vary over time very much, as we assume this is
  # not a parameter that is under (much) strategic control.
  hParms <- list(mu0 = c(5, 0, 2), # prior means 
                 sigma0 = matrix(c(5^2, 0, 0, 0, 0.4^2, 0, 0, 0, 3^2), ncol = 3), # 3 x 3 covariance matrix; note sds on the diagonal are squared to give variance.
                 sigmaQ = matrix(c(1^2, 0, 0, 0, 0.05^2, 0, 0, 0, 0.2^2), ncol = 3), # 3 x 3 covariance matrix; note sds on the diagonal are squared to give variance.
                 lims = matrix(c(0.1, 25, -60 * pi/180, 10 * pi/180, 0.01, 5), ncol = 3) # 2 x 2 matrix with lower and upper bounds in the two rows.
  )
  nPart <- 2000
  
  # Bootstrap particle filter parameter estimation for dataset
  source("likelihoodFunctions.R") # Functions for computing likelihood and doing MCMC fit.
  source("bootstrapPF.R") # Functions for bootstrap particle filter sequential estimation of parameters with each new trial.
  pfOut <- pfRun(data = allData, n = nPart, hyperParms = hParms)
  
  # Match the true data generating parameters for each trial with the inferred states.
  pfSummary <- summaryStats(pfOut) # Get useful summary statistics for each trial
  parmDF <- allData %>%
    filter(., time == 1) %>% # Select the data generating parameters for each trial
    select(., participant, trial, int, grad, noise) %>% # We don't need various columns
    dplyr::rename(., alpha = int, beta = grad, eta = noise) %>% # Rename generating parameters
    tidyr::pivot_longer(., cols = c("alpha", "beta", "eta"), 
                        names_to = "parameter", 
                        values_to = "trueVal") %>% # From wide to long
    left_join(., pfSummary, by=c("trial", "parameter")) # Join the two dataframes, based on trial and parameter columns
  
  ### Static parameter estimates (decisionBound_Stan.R) ###
  source("decisionBound_Stan.R")
  stanFit <- fitStanModel(allData = allData, nTrials = nTrials, accuracy = accuracy)
  stanParmSummary <- as.data.frame(stanFit$fitSummary)
  stanSamples <- stanFit$dfPostSamples
  
  # Manipulate the posterior samples and summary statistics to get in a useful format for plotting
  # For stanParmSummary we need to do some renaming to facilitate subsequent plotting
  stanParmSummary <- select(stanParmSummary, c("mean", "2.5%", "25%", "50%", "75%", "97.5%"))
  colnames(stanParmSummary) <- c("mean", "q_025", "q_25", "q_50", "q_75", "q_975")
  stanParmSummary$parameter <- rownames(stanParmSummary)
  stanParmSummary$trial <- max(allData$trial) + 10 # add a "virtual" trial ID
  stanSamples$trial <- max(allData$trial) + 10 # add a "virtual" trial ID
  stanSamples <- tidyr::pivot_longer(stanSamples,
                                     cols = c("alpha", "beta", "eta"), 
                                     names_to = "parameter", 
                                     values_to = "sampleVal")
  
  ### Quantify the goodness-of-fit: RMSE ###
  parmDF$staticPostMean <- rep(stanParmSummary$mean, times = nTrials)
  parmDF$staticPostMedian <- rep(stanParmSummary$q_50, times = nTrials)
  parmDF$staticSquareErrMean <- sqrt((parmDF$staticPostMean - parmDF$trueVal)^2)
  parmDF$staticSquareErrMedian <- sqrt((parmDF$staticPostMedian - parmDF$trueVal)^2)
  parmDF$pfSquareErrMean <- sqrt((parmDF$post_mean - parmDF$trueVal)^2)
  parmDF$pfSquareErrMedian <- sqrt((parmDF$q_5 - parmDF$trueVal)^2)
  errorSummary <- parmDF %>%
    group_by(., participant, parameter) %>%
    summarise(., rmse_StaticPostMean = mean(staticSquareErrMean, na.rm = TRUE),
              rmse_StaticPostMedian = mean(staticSquareErrMedian, na.rm = TRUE),
              rmse_pfPostMean = mean(pfSquareErrMean, na.rm = TRUE), 
              rmse_pfPostMedian = mean(pfSquareErrMedian, na.rm = TRUE))
  
  ### Visualisations ###
  
  # Search trajectory through the parameter space
  # Computed weighted average of all the locations visited
  uniqueLocInfo <- traject %>% group_by(., intercept, gradient) %>% 
    summarise(totalSamples = cumulativeSamples[max(visitCount)],
              totalVisits = max(visitCount),
              meanRewardRate = mean(meanObjective))
  wMean <- c(weighted.mean(uniqueLocInfo$gradient, uniqueLocInfo$totalSamples),
             weighted.mean(uniqueLocInfo$intercept, uniqueLocInfo$totalSamples))
  # Some hard-coded parameters that we know are suitable for the agents we've selected to illustrate.
  xyPeak <- c(-2.5*pi/180, 5) # Gradient and intercept with the "true" long-run maximum reward rate.
  # Make sure the size and colour scales accommodate the range across agents
  #sizeScaleLims <- c(1,35)
  #colScaleLims <- c(-0.06, 0.06)
  #legendPos <- c(0.1, 0.7)
  txtOffsets <- data.frame(ID = c("static", "limited", "extensive"),
                           xOffStart = c(0, 0, 0),
                           yOffStart = c(0, 0, 0),
                           xOffEnd = c(0, 0, 0),
                           yOffEnd = c(0, 0, 0))

  staticSearchPlot <- ggplot() +
    geom_path(aes(x = gradient, y = intercept), data = traject, linewidth = 2, colour = "grey", alpha = 1) + # trajectory
    geom_point(aes(x = gradient, y = intercept, colour = meanRewardRate, size = totalVisits), data = uniqueLocInfo) + # sampled points
    geom_point(aes(x = xyPeak[1], y = xyPeak[2]), shape = 3, colour = "black", size = 3) + # peak of the objective
    geom_point(aes(x = wMean[1], y = wMean[2]), shape = 17, colour = "black", size = 3) + # weighted average of the visited locations
    # Play around with the position of the annotations by tweaking the offsets
    annotate("text", x = traject$gradient[1] + txtOffsets$xOffStart[txtOffsets$ID == pID], y = traject$intercept[1] + txtOffsets$yOffStart[txtOffsets$ID == pID], label = "start") +
    annotate("text", x = traject$gradient[nrow(traject)] + txtOffsets$xOffEnd[txtOffsets$ID == pID], y = traject$intercept[nrow(traject)] + txtOffsets$yOffEnd[txtOffsets$ID == pID], label = "end") +
    scale_size_area(name = "# visits") + # limits = sizeScaleLims, breaks=c(5, 25, 125, 500)) +
    scale_colour_viridis(na.value="transparent") + # limits = colScaleLims, breaks=c(-0.06, -0.03, 0, 0.03, 0.06)) +
    scale_x_continuous(breaks = seq(-1, 0, by = 0.25), limits = c(-1.2, 0.2)) +
    scale_y_continuous(breaks = seq(0, 25, by = 5), limits = c(0, 26)) +
    labs(x = TeX("$\\beta$"), y = TeX("$\\alpha$"), colour = "reward rate") +
    theme_bw() +
    #guides(size = guide_legend(position = "inside"), colour = guide_colorbar(position = "inside")) +
    #theme(legend.position.inside = legendPos) +
    # theme(aspect.ratio = 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.x  = element_text(size=12)) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.title.x.bottom = element_text(size=14,face="bold")) +
    theme(axis.title.y.left = element_text(size=14,face="bold"))
    # Optional: show the mean posterior estimates in the plot of the search trajectory
    # geom_point(aes(x = stanParmSummary$mean[2], y = stanParmSummary$mean[1]), shape = 2, colour = "black", size = 3) + # weighted average of the visited locations
    # geom_errorbar(aes(x = stanParmSummary$mean[2], ymin = stanParmSummary$q_025[1], ymax = stanParmSummary$q_975[1]), colour = "black", size = 0.5, width = 0.1) +
    # geom_errorbarh(aes(y = stanParmSummary$mean[1], xmin = stanParmSummary$q_025[2], xmax = stanParmSummary$q_975[2]), colour = "black", size = 0.5, height = 0.1)
  staticSearchPlot
  
  # Comparison between true and estimated parameters
  estParmPlot <- ggplot(data = parmDF) +
    geom_ribbon(aes(x = trial, ymin = q_025, ymax = q_975), colour = "grey90", fill = "grey90", alpha = 1) + # 95% credible interval
    geom_line(aes(x = trial, y = trueVal), linewidth = 1, colour = "grey50", show.legend = TRUE) +
    geom_line(aes(x = trial, y = post_mean), colour = "black", linewidth = 1, alpha = 1, show.legend = TRUE) + # Posterior mean
    # To the right, show the estimates from a static model (posterior parameters from Stan)
    geom_violin(data = stanSamples, aes(x = trial, y = sampleVal), trim = FALSE, scale = "area", width = 15, fill = "grey90", colour = "black", alpha = 1)+
    geom_errorbar(data = stanParmSummary, aes(x = trial, ymin = q_25, ymax = q_75), size = 1, colour = "black")+ # IQR
    geom_point(data = stanParmSummary, aes(x = trial, y = mean), size = 2, shape = 15, colour= "grey50", position = position_dodge(width = 1), alpha = 1)+
    facet_grid(parameter~., scales = "free_y") +
    # scale_x_continuous(breaks = c(0, 50, 100, 150, 200), labels = c(0, 50, 100, 150, 200)) +
    labs(x = 'trial', y='parameter value') +
    theme_bw()+
    theme(strip.text=element_text(size=14,face="bold"))+
    theme(axis.text.x  = element_text(size=12))+
    theme(axis.text.y = element_text(size=12))+
    theme(axis.title.x.bottom = element_text(size=14,face="bold"))+
    theme(axis.title.y.left = element_text(size=14,face="bold"))+
    theme(legend.text = element_text(size=12))+
    theme(legend.title = element_text(size=14, face="bold"))
estParmPlot

# Return a list with 
#   - Static parameter estimates
#   - Particle filter parameter estimates
#   - RMSE from static and particle filter estimates
#   - 2D trajectory plot (intercept, gradient)
#   - 3 x 1 plot of true vs estimated parameters
return(list(staticParms = stanParmSummary,
            staticSamples = stanSamples, 
            pfParms = parmDF,
            RMSE = errorSummary,
            trajectory = traject,
            uniqueLocs = uniqueLocInfo,
            meanLocs = wMean,
            searchPlot = staticSearchPlot,
            parmPlot = estParmPlot)
       )
}