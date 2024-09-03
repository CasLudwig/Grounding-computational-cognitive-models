# C.J.H. Ludwig, July 2023
# Edited by E. Stuchly in January 2024
# Local sampling and comparison routine.
# In this version of the algorithm, we assume metric knowledge of the objective, prior knowledge of the peak, a bivariate Gaussian jump distribution, and a linear (clipped) jump probability function.
# arguments: standard distribution of the "jump" probability;
# time limit (number of time/evidence samples) for which the trajectory/simulated experiment will run
# window size;
# upper and lower range of intercept values;
# upper and lower range of gradient values;
# gaussian noise around the current state of accumulated evidence;
# satisficing constant (1 = finding the objective peak reward rate);
# the value of the objective peak reward rate;

LoSaCo_decBound <- function(stepParms = c(3,pi/18), nMaxSteps = 1200, W = 5, intrange = c(0,25), gradrange = c(-pi/3,pi/36), evinoise = 1, satisfice = 1, peak = 0.05){
  satisfice_peak = satisfice*peak #first, compute what outcome is deemed as "acceptable" by the agent, in units of reward rate
  # Set up a dataframe that's going to contain the results of the parameter search, where each row corresponds to a sampling period with current parameter values for W trials
  paramsGrid <- data.frame(intercept = 0,
                       gradient = 0,
                       meanObjective = 0,
                       nTrials = W,
                       cumulativeSamples = 0,
                       jumpP = 0,
                       visitCount = 0)
  
  ### Sample objective information from the starting location
  winCount <- 1 # Pointer to the row in paramsGrid that contains the current location
  # Start at a random (uniform) location inside the parameter space/grid
  paramsGrid$intercept[winCount] <- round(runif(1, intrange[1], intrange[2]),5)
  paramsGrid$gradient[winCount] <- round(runif(1, gradrange[1], gradrange[2]),5)
  
  rrate_estim <- rr_trialwin(window=paramsGrid$nTrials, MaxTrialsamples=MaxTrialsamples, intercept = paramsGrid$intercept, gradient = paramsGrid$gradient, evinoise = evinoise, drift = drift, iti_cor = iti_cor, iti_pen = iti_pen, monpen = monpen, returnTrials=1,0) #runs the reward rate estimation, based on a window of trials with a particular set of threshold parameters (intercept, gradient, evinoise)
  paramsGrid$meanObjective[winCount] <- rrate_estim$rr_mean #mean of the reward rate estimate
  paramsGrid$cumulativeSamples[winCount] <- rrate_estim$nsampleswin #keep track of how many (time/evidence) samples were presented in the window of trials
  
  rrate_estim$data_window$win <- winCount #put the window index into the data frame
  data_gen <- rrate_estim$data_window #create a dataframe that will contain trial data from all the windows

  # Now that we've sampled the objective, set a "jump" probability for the next sample. If the estimated objective has met the stopping (satisficing) criterion, we do not jump. We could of course choose a small non-zero value here so that there is always some (small) probability of further exploration.
  jumpProb <- 1 - (paramsGrid$meanObjective[winCount]/satisfice_peak) # Jump probability approaches 0 as we get close to our stopping criterion
  if(jumpProb < 0){jumpProb <- 0} # This will happen if the stopping criterion has been met
  if(jumpProb > 1){jumpProb <- 1} # This will happen if the objective estimate is negative, in which case, we definitely want to jump
  paramsGrid$jumpP[winCount] <- jumpProb
  paramsGrid$visitCount[winCount] <- 1 #indicate that this particular spot in parameter space has been visited
  
  # Set the parameters for the second sampling period
  paramsGrid[winCount+1,] <- 0
  paramsGrid$intercept[winCount+1] <- paramsGrid$intercept[winCount]
  paramsGrid$gradient[winCount+1] <- paramsGrid$gradient[winCount]
  paramsGrid$nTrials[winCount+1] <- W
  returnToPrevious <- FALSE
  stepCount = paramsGrid$cumulativeSamples[winCount]
  
  ### Now enter a loop in which the agent performs a greedy search for all subsequent sampling periods
  while(stepCount < nMaxSteps){
    
    winCount = winCount+1 #increase the window counter
    #append another empty row to the paramsGrid dataframe, to contain data from the next trial window
    paramsGrid[winCount+1,] <- 0
    paramsGrid$nTrials[winCount] <- W
    paramsGrid$visitCount[winCount] <- 1
    
    # See if we need to take a step to a different spot in the parameter space. If the most recently sampled value was worse than the previous one, we want to either return to the previous location and sample there, or take a step elsewhere from that previous location (rather than the from the most recently sampled location).
    takeStep <- rbinom(1,1,jumpProb)
    
    #if the agent is "static" (0 jump distribution SD), prevent them from attempting to perform a jump
    if (stepParms[1]==0 & stepParms[2]==0){
      takeStep=0
    }
    
    if(takeStep == 0){ # If we are staying at the current location
      rrate_estim <- rr_trialwin(window=paramsGrid$nTrials[winCount], MaxTrialsamples=MaxTrialsamples, intercept = paramsGrid$intercept[winCount], gradient = paramsGrid$gradient[winCount], evinoise = evinoise, drift = drift, iti_cor = iti_cor, iti_pen = iti_pen, monpen = monpen, returnTrials=1,0) #runs the reward rate estimation, based on a window of trials with a particular set of threshold parameters (intercept, gradient, evinoise) and task parameters
      paramsGrid$meanObjective[winCount] <- rrate_estim$rr_mean #reward rate estimate from the current sampling period
      paramsGrid$cumulativeSamples[winCount] <- paramsGrid$cumulativeSamples[winCount-1] + rrate_estim$nsampleswin #keep track of how many total time/evidence samples have been taken since the beginning
      stepCount = paramsGrid$cumulativeSamples[winCount]
      rrate_estim$data_window$win <- winCount #put the window (i.e., the sampling period) index into the data frame
      data_gen <- rbind(data_gen,rrate_estim$data_window) #append the latest simulated trial data into the dataframe
      
      # if the next sampling period continues from the previously visited location, increase its "visitation counter"; otherwise remain at the current location and increase its visit counter
      if(returnToPrevious){
        paramsGrid$visitCount[winCount] <- paramsGrid$visitCount[winCount-2] + 1
      }else{
        paramsGrid$visitCount[winCount] <- paramsGrid$visitCount[winCount-1] + 1
      }
      
      # Set the start location for the next sampling epoch and compute the jump probability
      if(stepCount < nMaxSteps){
        paramsGrid$intercept[winCount+1] <- paramsGrid$intercept[winCount]
        paramsGrid$gradient[winCount+1] <- paramsGrid$gradient[winCount]
        # Compute the jump probability for the next step
        jumpProb <- 1 - (paramsGrid$meanObjective[winCount]/satisfice_peak) # Jump probability approaches 0 as we get close to our stopping criterion
        if(jumpProb < 0){jumpProb <- 0} # This will happen if the stopping criterion has been met
        if(jumpProb > 1){jumpProb <- 1} # This will happen if the objective estimate is negative, in which case, we definitely want to jump
      }else{
        jumpProb <- 0
      }
      paramsGrid$jumpP[winCount] <- jumpProb
      returnToPrevious <- FALSE
      
    }else{ # when takeStep == 1; make a "jump" to elsewhere in the parameter space, centered around the current location
      paramsNew <- round(rnorm(2,mean = c(paramsGrid$intercept[winCount], paramsGrid$gradient[winCount]), sd = stepParms),5)
      #introduce a loop that keeps proposing new parameters until their values fall within the accepted range
      while (paramsNew[1]<intrange[1] | paramsNew[1]>intrange[2] | paramsNew[2]<gradrange[1] | paramsNew[2]>gradrange[2]){
        paramsNew <- round(rnorm(2,mean = c(paramsGrid$intercept[winCount], paramsGrid$gradient[winCount]), sd = stepParms),5)
      }
      # Sample the objective at the new location
      paramsGrid$intercept[winCount] <- paramsNew[1]
      paramsGrid$gradient[winCount] <- paramsNew[2]
      
      rrate_estim <- rr_trialwin(window=paramsGrid$nTrials[winCount], MaxTrialsamples=MaxTrialsamples, intercept = paramsGrid$intercept[winCount], gradient = paramsGrid$gradient[winCount], evinoise = evinoise, drift = drift, iti_cor = iti_cor, iti_pen = iti_pen, monpen = monpen, returnTrials=1,0) #runs the reward rate estimation, based on a window of trials with a particular set of threshold parameters (intercept, gradient, evinoise)
      paramsGrid$meanObjective[winCount] <- rrate_estim$rr_mean #reward rate estimate from the current sampling period
      paramsGrid$cumulativeSamples[winCount] <- paramsGrid$cumulativeSamples[winCount-1] + rrate_estim$nsampleswin #keep track of how many (time/evidence) samples were presented in the window of trials
      stepCount = paramsGrid$cumulativeSamples[winCount]
      rrate_estim$data_window$win <- winCount #put the window (i.e., the sampling period) index into the data frame
      data_gen <- rbind(data_gen,rrate_estim$data_window) #append the latest simulated trial data to the dataframe
      
      # Now that we've sampled the objective, test whether the estimated value from the current sampling period is higher/lower than the previously sampled value. If it is lower, we want to continue our search from the previous location. If it's better, we want to continue our search from the current location.
      if(stepCount < nMaxSteps){
        if(paramsGrid$meanObjective[winCount] < paramsGrid$meanObjective[winCount-1]){ # Current objective is worse than the previous one, so we'll continue from the previous location/make a new jump from that previous location. In this case, we need not compute a jump probability, because it would be the same as the last time we sampled the previous location.
          returnToPrevious <- TRUE
          paramsGrid$intercept[winCount+1] <- paramsGrid$intercept[winCount-1]
          paramsGrid$gradient[winCount+1] <- paramsGrid$gradient[winCount-1]
        }else{ # Current objective is better than the previous one, so we'll continue from the current location. Recompute the jump probability.
          returnToPrevious <- FALSE
          paramsGrid$intercept[winCount+1] <- paramsGrid$intercept[winCount]
          paramsGrid$gradient[winCount+1] <- paramsGrid$gradient[winCount]
          # Compute the jump probability for the next step
          jumpProb <- 1 - (paramsGrid$meanObjective[winCount]/satisfice_peak) # Jump probability approaches 0 as we get close to our stopping criterion
          if(jumpProb < 0){jumpProb <- 0} # This will happen if the stopping criterion has been met
          if(jumpProb > 1){jumpProb <- 1} # This will happen if the objective estimate is negative, in which case, we definitely want to jump
        }
      }
      else{ #if the maximum number of allowed samples (per experiment) has been exceeded
        jumpProb <- 0
      }
      paramsGrid$jumpP[winCount] <- jumpProb
      
    } # if takeStep == 0/1
  } # for winCount = 1:nMaxSteps loop

  # Compute the cumulative *actual* reward rate values collected (per sampling period, so that this number lives on the same scale as the objective function)
  paramsGrid$cumulativeObjective <- cumsum(paramsGrid$meanObjective * paramsGrid$nTrials)/cumsum(paramsGrid$nTrials)
  
  data_gen$trialid <- data_gen$trial
  data_gen$trial <- data_gen$trial + (data_gen$win-1)*W
  
  paramsGrid <- filter(paramsGrid,cumulativeSamples<nMaxSteps & cumulativeSamples!=0) #filter out the parameters grid to only include rows from complete windows (i.e., those where the nMaxSteps was not exceeded)
  return(list(paramsGrid,data_gen)) # ,params_overall
}
