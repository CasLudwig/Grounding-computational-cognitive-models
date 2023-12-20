# C.J.H. Ludwig, July 2023
# Local sampling and comparison routine.

LoSaCo1 <- function(sd0 = 1, stepParms = c(1,1), nMaxSteps = NULL, xrange = NULL, yrange = NULL, N = NULL, W = NULL, satisfice = 100, muLoc = c(0,0), width = c(1,1), rho = 0, peak = 1){
  # In this version of the algorithm, we assume metric knowledge of the objective, prior knowledge of the peak, a bivariate Gaussian jump distribution, and a linear (clipped) jump probability function.
  nMaxSteps <- min(floor(N/W), nMaxSteps) # Check whether the supplied maximum number of steps is compatible with the window size.
  nRemainTrials <- N - nMaxSteps*W # How many trials are left after the last sampling window (> 0 only if N is not an integer multiple of W).
  # If there are any remaining trials, add an extra step, but with a smaller window of trials (i.e. lastW < W).
  if(nRemainTrials > 0){
    nMaxSteps <- nMaxSteps + 1
    lastW <- nRemainTrials
  }
  
  # Set up a dataframe that's going to contain the results of our search
  xyGrid <- data.frame(theta1 = rep(0,nMaxSteps),
                       theta2 = rep(0,nMaxSteps),
                       meanObjective = rep(0,nMaxSteps),
                       sdObjective = rep(0,nMaxSteps),
                       estObjective = rep(0,nMaxSteps),
                       nSamples = rep(W,nMaxSteps),
                       cumulativeSamples = rep(0,nMaxSteps),
                       jumpP = rep(0,nMaxSteps))
  if(nRemainTrials > 0){
    xyGrid$nSamples[nMaxSteps] <- lastW
  }
  
  # Sample objective information from the starting location
  stepCount <- 1 # Pointer to the row in xyGrid that contains the current location
  # Start at a random (uniform) location inside the grid
  xyGrid$theta1[stepCount] <- runif(1, xrange[1], xrange[2]) 
  xyGrid$theta2[stepCount] <- runif(1, yrange[1], yrange[2])
  xyGrid$meanObjective[stepCount] <- obj2D_uniModal(x=xyGrid$theta1[stepCount], y=xyGrid$theta2[stepCount], mu=muLoc, sigma=width, corrxy=rho, amp=peak) # Mean objective value
  xyGrid$sdObjective[stepCount] <- sd0/sqrt(xyGrid$nSamples[stepCount]) # The standard deviation around the mean objective value---depends on the number of trials in a window (larger windows result in more precise estimates).
  xyGrid$estObjective[stepCount] <- rnorm(1, xyGrid$meanObjective[stepCount], xyGrid$sdObjective[stepCount])
  xyGrid$cumulativeSamples[stepCount] <- xyGrid$cumulativeSamples[stepCount] + xyGrid$nSamples[stepCount]
  
  # Now that we've sampled the objective, set a "jump" probability for the next sample. If the estimated objective has met the stopping (satisficing) criterion, we do not jump. We could of course choose a small non-zero value here so that there is always some (small) probability of further exploration.
  jumpProb <- 1 - (xyGrid$estObjective[stepCount]/satisfice) # Jump probability approaches 0 as we get close to our stopping criterion
  if(jumpProb < 0){jumpProb <- 0} # This will happen if the stopping criterion has been met
  if(jumpProb > 1){jumpProb <- 1} # This will happen if the objective estimate is negative, in which case, we definitely want to jump
  xyGrid$jumpP[stepCount] <- jumpProb
  
  # Set the start location for the second sample
  xyGrid$theta1[stepCount+1] <- xyGrid$theta1[stepCount]
  xyGrid$theta2[stepCount+1] <- xyGrid$theta2[stepCount]
  returnToPrevious <- FALSE
  
  # Now enter a loop in which the agent performs a greedy search
  for(stepCount in 2:nMaxSteps){
    # See if we need to take a step from the current best value. If the most recently sampled value was worse than the best, we want to either return to the previous best location or take a step from that best location (rather than the most recently sampled location).
    takeStep <- rbinom(1,1,jumpProb)
    
    if(takeStep == 0){ # Stay at the current location
      xyGrid$meanObjective[stepCount] <- obj2D_uniModal(x=xyGrid$theta1[stepCount], y=xyGrid$theta2[stepCount], mu=muLoc, sigma=width, corrxy=rho, amp=peak) # Mean objective value
      xyGrid$sdObjective[stepCount] <- sd0/sqrt(xyGrid$nSamples[stepCount]) # The standard deviation around the mean objective value — depends on the number of trials in a window (larger windows result in more precise estimates).
      xyGrid$estObjective[stepCount] <- rnorm(1, xyGrid$meanObjective[stepCount], xyGrid$sdObjective[stepCount])
      
      # Keep track of the total number of samples in a given location. We haven't moved on this step, but we may either stayed in the most recent location or we have have returned to a previous location.
      if(returnToPrevious){
        xyGrid$cumulativeSamples[stepCount] <- xyGrid$cumulativeSamples[stepCount-2] + xyGrid$nSamples[stepCount]
      }else{
        xyGrid$cumulativeSamples[stepCount] <- xyGrid$cumulativeSamples[stepCount-1] + xyGrid$nSamples[stepCount]
      }
      
      # Set the start location for the next sampling epoch and compute the jump probability
      if(stepCount < nMaxSteps){
        xyGrid$theta1[stepCount+1] <- xyGrid$theta1[stepCount]
        xyGrid$theta2[stepCount+1] <- xyGrid$theta2[stepCount]
        # Compute the jump probability for the next step
        jumpProb <- 1 - (xyGrid$estObjective[stepCount]/satisfice) # Jump probability approaches 0 as we get close to our stopping criterion
        if(jumpProb < 0){jumpProb <- 0} # This will happen if the stopping criterion has been met
        if(jumpProb > 1){jumpProb <- 1} # This will happen if the objective estimate is negative, in which case, we definitely want to jump
      }else{
        jumpProb <- 0
      }
      xyGrid$jumpP[stepCount] <- jumpProb
      returnToPrevious <- FALSE
      
    }else{ # takeStep == 1; select a location randomly from the current location
      
      xyNew <- rnorm(2,mean = c(xyGrid$theta1[stepCount], xyGrid$theta2[stepCount]), sd = stepParms)
      # Sample the objective at the new location
      xyGrid$theta1[stepCount] <- xyNew[1]
      xyGrid$theta2[stepCount] <- xyNew[2]
      xyGrid$meanObjective[stepCount] <- obj2D_uniModal(x=xyNew[1], y=xyNew[2], mu=muLoc, sigma=width, corrxy=rho, amp=peak) # Mean objective value
      xyGrid$sdObjective[stepCount] <- sd0/sqrt(xyGrid$nSamples[stepCount]) # The standard deviation around the mean objective value---depends on the number of trials in a window (larger windows result in more precise estimates).
      xyGrid$estObjective[stepCount] <- rnorm(1, xyGrid$meanObjective[stepCount], xyGrid$sdObjective[stepCount])
      xyGrid$cumulativeSamples[stepCount] <- xyGrid$nSamples[stepCount]
      # Now that we've sampled the objective, test whether the estimated value is higher/lower than the previously sampled value. If it is worse, we want to continue our search from the previous location. If it's better, we want to continue our search from the current location.
      if(stepCount < nMaxSteps){
        if(xyGrid$estObjective[stepCount] < xyGrid$estObjective[stepCount-1]){ # Current objective is worse than the previous one, so we'll continue from the previous location. In this case, we need not compute a jump probability, because it would be the same as the last time we sampled the previous location.
          returnToPrevious <- TRUE
          xyGrid$theta1[stepCount+1] <- xyGrid$theta1[stepCount-1]
          xyGrid$theta2[stepCount+1] <- xyGrid$theta2[stepCount-1]
        }else{ # Current objective is better than the previous one, so we'll continue from the current location. Recompute the jump probability.
          returnToPrevious <- FALSE
          xyGrid$theta1[stepCount+1] <- xyGrid$theta1[stepCount]
          xyGrid$theta2[stepCount+1] <- xyGrid$theta2[stepCount]
          # Compute the jump probability for the next step
          jumpProb <- 1 - (xyGrid$estObjective[stepCount]/satisfice) # Jump probability approaches 0 as we get close to our stopping criterion
          if(jumpProb < 0){jumpProb <- 0} # This will happen if the stopping criterion has been met
          if(jumpProb > 1){jumpProb <- 1} # This will happen if the objective estimate is negative, in which case, we definitely want to jump
        }
      }
      else{
        jumpProb <- 0
      }
      xyGrid$jumpP[stepCount] <- jumpProb
      
    } # if takeStep == 0/1
    
  } # for stepCount = 1:nMaxSteps loop
  
  # Compute the cumulative *actual* objective values collected (per trial, so that this number lives on the same scale as the objective function)
  xyGrid$cumulativeObjective <- cumsum(xyGrid$meanObjective * xyGrid$nSamples)/cumsum(xyGrid$nSamples)
  
  return(xyGrid)
}

# Set up a small function that turns a continuous (noisy) objective estimate into a discrete, categorical label. Input vector b specifies the category boundaries relative to the peak (i.e. as proportions).
categoricalObjective <- function(x, b = c(0.2, 0.4, 0.6, 0.8), peak = 1){
  tmp1 <- x - b*peak
  tmp2 <- which(tmp1 <= 0)
  if(any(tmp2)){
    catObj <- tmp2[1]
  }else{catObj <- (length(b) +1)}
  return(catObj)
}

LoSaCo2 <- function(sd0 = 1, stepParms = c(1,1), nMaxSteps = NULL, xrange = NULL, yrange = NULL, N = NULL, W = NULL, muLoc = c(0,0), width = c(1,1), rho = 0, peak = 1){
  # In this version of the algorithm, we assume coarse, categorical knowledge of the objective with linearly spaced bins over the objective scale and a linear mapping to the jump probabilities.
  nMaxSteps <- min(floor(N/W), nMaxSteps) # Check whether the supplied maximum number of steps is compatible with the window size.
  nRemainTrials <- N - nMaxSteps*W # How many trials are left after the last sampling window (> 0 only if N is not an integer multiple of W).
  # If there are any remaining trials, add an extra step, but with a smaller window of trials (i.e. lastW < W).
  if(nRemainTrials > 0){
    nMaxSteps <- nMaxSteps + 1
    lastW <- nRemainTrials
  }
  
  # Specify the edges for the categorical objective values (relative to the peak, i.e. as proportions).
  # Note that any values below 0.2 (including negative values due to noise) are assigned label 1; any values above 0.8 (including values above the peak due to noise) are assigned label 5.
  binEdges <- c(0.2, 0.4, 0.6, 0.8)
  maxCat <- length(binEdges) + 1 # best possible value we can get
  piStep <- 1/length(binEdges) # size of the steps by which the jump probability decreases over the categorical objective value scale
  
  # Set up a dataframe that's going to contain the results of our search
  xyGrid <- data.frame(theta1 = rep(0,nMaxSteps),
                       theta2 = rep(0,nMaxSteps),
                       meanObjective = rep(0,nMaxSteps),
                       sdObjective = rep(0,nMaxSteps),
                       estObjective = rep(0,nMaxSteps),
                       nSamples = rep(W,nMaxSteps),
                       cumulativeSamples = rep(0,nMaxSteps),
                       jumpP = rep(0,nMaxSteps))
  if(nRemainTrials > 0){
    xyGrid$nSamples[nMaxSteps] <- lastW
  }
  
  # Sample objective information from the starting location
  stepCount <- 1 # Pointer to the row in xyGrid that contains the current location
  # Start at a random (uniform) location inside the grid
  xyGrid$theta1[stepCount] <- runif(1, xrange[1], xrange[2]) 
  xyGrid$theta2[stepCount] <- runif(1, yrange[1], yrange[2])
  xyGrid$meanObjective[stepCount] <- obj2D_uniModal(x=xyGrid$theta1[stepCount], y=xyGrid$theta2[stepCount], mu=muLoc, sigma=width, corrxy=rho, amp=peak) # Mean objective value
  xyGrid$sdObjective[stepCount] <- sd0/sqrt(xyGrid$nSamples[stepCount]) # The standard deviation around the mean objective value---depends on the number of trials in a window (larger windows result in more precise estimates).
  xyGrid$estObjective[stepCount] <- categoricalObjective(rnorm(1, xyGrid$meanObjective[stepCount], xyGrid$sdObjective[stepCount]), b = binEdges, peak = peak)
  xyGrid$cumulativeSamples[stepCount] <- xyGrid$cumulativeSamples[stepCount] + xyGrid$nSamples[stepCount]
  
  # Now that we've sampled the objective, set a "jump" probability for the next sample. This is simply a linearly decreasing function of the categorical objective estimate, with slope piStep
  jumpProb <- 1 - (xyGrid$estObjective[stepCount] - 1) * piStep # So jumpProb = 1 when c = 1; jumpProb = 0 when c = 5.
  xyGrid$jumpP[stepCount] <- jumpProb
  
  # Set the start location for the second sample
  xyGrid$theta1[stepCount+1] <- xyGrid$theta1[stepCount]
  xyGrid$theta2[stepCount+1] <- xyGrid$theta2[stepCount]
  returnToPrevious <- FALSE
  
  # Now enter a loop in which the agent performs a greedy search
  for(stepCount in 2:nMaxSteps){
    # See if we need to take a step from the current best value. If the most recently sampled value was worse than the best, we want to either return to the previous best location or take a step from that best location (rather than the most recently sampled location).
    takeStep <- rbinom(1,1,jumpProb)
    
    if(takeStep == 0){ # Stay at the current location
      xyGrid$meanObjective[stepCount] <- obj2D_uniModal(x=xyGrid$theta1[stepCount], y=xyGrid$theta2[stepCount], mu=muLoc, sigma=width, corrxy=rho, amp=peak) # Mean objective value
      xyGrid$sdObjective[stepCount] <- sd0/sqrt(xyGrid$nSamples[stepCount]) # The standard deviation around the mean objective value — depends on the number of trials in a window (larger windows result in more precise estimates).
      xyGrid$estObjective[stepCount] <- categoricalObjective(rnorm(1, xyGrid$meanObjective[stepCount], xyGrid$sdObjective[stepCount]), b = binEdges, peak = peak)
      
      # Keep track of the total number of samples in a given location. We haven't moved on this step, but we may either stayed in the most recent location or we have have returned to a previous location.
      if(returnToPrevious){
        xyGrid$cumulativeSamples[stepCount] <- xyGrid$cumulativeSamples[stepCount-2] + xyGrid$nSamples[stepCount]
      }else{
        xyGrid$cumulativeSamples[stepCount] <- xyGrid$cumulativeSamples[stepCount-1] + xyGrid$nSamples[stepCount]
      }
      
      # Set the start location for the next sampling epoch and compute the jump probability
      if(stepCount < nMaxSteps){
        xyGrid$theta1[stepCount+1] <- xyGrid$theta1[stepCount]
        xyGrid$theta2[stepCount+1] <- xyGrid$theta2[stepCount]
        # Compute the jump probability for the next step
        jumpProb <- 1 - (xyGrid$estObjective[stepCount] - 1) * piStep
      }
      xyGrid$jumpP[stepCount] <- jumpProb
      returnToPrevious <- FALSE
      
    }else{ # takeStep == 1; select a location randomly from the current location
      
      xyNew <- rnorm(2,mean = c(xyGrid$theta1[stepCount], xyGrid$theta2[stepCount]), sd = stepParms)
      # Sample the objective at the new location
      xyGrid$theta1[stepCount] <- xyNew[1]
      xyGrid$theta2[stepCount] <- xyNew[2]
      xyGrid$meanObjective[stepCount] <- obj2D_uniModal(x=xyNew[1], y=xyNew[2], mu=muLoc, sigma=width, corrxy=rho, amp=peak) # Mean objective value
      xyGrid$sdObjective[stepCount] <- sd0/sqrt(xyGrid$nSamples[stepCount]) # The standard deviation around the mean objective value---depends on the number of trials in a window (larger windows result in more precise estimates).
      xyGrid$estObjective[stepCount] <- categoricalObjective(rnorm(1, xyGrid$meanObjective[stepCount], xyGrid$sdObjective[stepCount]), b = binEdges, peak = peak)
      xyGrid$cumulativeSamples[stepCount] <- xyGrid$nSamples[stepCount]
      # Now that we've sampled the objective, test whether the estimated value is higher/lower than the previously sampled value. If it is worse, we want to continue our search from the previous location. If it's better, we want to continue our search from the current location.
      if(stepCount < nMaxSteps){
        if(xyGrid$estObjective[stepCount] < xyGrid$estObjective[stepCount-1]){ # Current objective is worse than the previous one, so we'll continue from the previous location. In this case, we need not compute a jump probability, because it would be the same as the last time we sampled the previous location.
          returnToPrevious <- TRUE
          xyGrid$theta1[stepCount+1] <- xyGrid$theta1[stepCount-1]
          xyGrid$theta2[stepCount+1] <- xyGrid$theta2[stepCount-1]
        }else{ # Current objective is better than the previous one, so we'll continue from the current location. Recompute the jump probability.
          returnToPrevious <- FALSE
          xyGrid$theta1[stepCount+1] <- xyGrid$theta1[stepCount]
          xyGrid$theta2[stepCount+1] <- xyGrid$theta2[stepCount]
          # Compute the jump probability for the next step
          jumpProb <- 1 - (xyGrid$estObjective[stepCount] - 1) * piStep
        }
      }
      else{
        jumpProb <- 0
      }
      xyGrid$jumpP[stepCount] <- jumpProb
      
    } # if takeStep == 0/1
    
  } # for stepCount = 1:nMaxSteps loop
  
  # Compute the cumulative *actual* objective values collected (per trial, so that this number lives on the same scale as the objective function)
  xyGrid$cumulativeObjective <- cumsum(xyGrid$meanObjective * xyGrid$nSamples)/cumsum(xyGrid$nSamples)
  
  return(xyGrid)
}