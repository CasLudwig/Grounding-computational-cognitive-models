# C.J.H. Ludwig, January 2024
# Functions to compute likelihood for single and multiple trials. These
# functions are the engines for parameter estimation through MCMC and particle
# filtering. The data come from an expanded judgment paradigm so that for each
# trial we have a sequence of evidence samples, 'wait' decisions and a final
# 'go' decision. The model assumes a decision threshold with a certain height
# (alpha) and gradient (beta) and perfect integration of evidence. At each
# time-point, the agent generates a probabilistic wait/go decision depending on
# how close the evidence is to the decision threshold. The extent to which the
# agent behaves in a probabilistic or deterministic manner is controlled by a
# noise parameter (eta). So there are three parameters (alpha, beta and eta).

# Log-likelihood for a single trial. This is for simulated data where a correct response is made when the cumulative evidence is positive. By default it returns the log-likelihood, but if returnLogLik == FALSE, the function returns the trial data augmented with various model-generated quantities (boundary values, action probabilities).
logLikOneTrial <- function(params, trialData, returnLogLik = TRUE){

  alpha <- params[1]
  beta <- params[2]
  eta <- params[3]
  
  # Add a row at the beginning of trial data to represent the start state (time = 0, evidence = 0). No action is taken in this state (at least, not in these simulated data). In other words, the agent will always wait for at least one sample. However, the model will assign a probability of going/waiting in this state.
  trialData <- rbind(trialData[1,], trialData)
  trialData$action[1] <- 0
  trialData$time[1] <- 0
  trialData$evi[1] <- 0
  
  nSamples <- nrow(trialData) # Number of samples observed
  if(trialData$action[nSamples]==1){ #if the trial isn't timed-out
    if (trialData$evi[trialData$action==1]>0){
      correct <- 1
    } else {
      correct <- 0
    }

    # Compute upper decision boundary at each time point (including 0, before any evidence is observed).
    trialData$boundary <- alpha + trialData$time * tan(beta)
    trialData$boundary[trialData$boundary < 0] <- 0
    
    # Now we are in a position to compute probabilities for each of the three actions at each time point.
    trialData$p_up <- 1 - pnorm(trialData$boundary, trialData$evi, eta)
    trialData$p_down <- pnorm(-trialData$boundary, trialData$evi, eta)
    trialData$p_wait <- rep(1,nSamples) - trialData$p_up - trialData$p_down
    
    # Likelihood of each wait/go action.
    if(correct == 1){
      trialData$lik <- (rep(1,nSamples) - trialData$action) * trialData$p_wait + trialData$action * trialData$p_up
    }else{
      trialData$lik <- (rep(1,nSamples) - trialData$action) * trialData$p_wait + trialData$action * trialData$p_down
    }
    # Deal with the usual difficulties (probabilities of 0 and 1); these will only occur if the evidence is >8*s away from the nearest boundary. We only care about these if these values happen for the likelihood of the observed actions.
    trialData$lik[trialData$lik <= 0] <- 1e-12
    trialData$lik[trialData$lik >= 1] <- 1 - 1e-12
    
    # Log-likelihood for the entire observed action sequence
    logLikSeq <- sum(log(trialData$lik))
    
  } else { # Time out; I'm setting these to NA and we should just avoid
    logLikSeq <- NA
  }
  
  # Check what to return
  if(returnLogLik){
    return(logLikSeq)
  } else{return(trialData)} # If we don't want the log-likelihood, return the data for the trial (useful for checking the likelihood calculation for an individual trial).
}


# Log-likelihood for an entire dataset. This is for simulated data where a correct response is made when the cumulative evidence is positive. By default it returns the (negative) log-likelihood.
logLikOverall <- function(params, data){
  
  trials <- unique(data$trial)
  nTrials <- length(trials)
  
  logLiks <- rep(NA, nTrials)
  
  # Compute the likelihood for each trial
  for(i in 1:nTrials){
    trialData <- filter(data, trial == trials[i]) # select evidence and action sequences for the current trial
    logLiks[i] <- logLikOneTrial(params = params, trialData = trialData) # return log-likelihood for current trial
  }
  
  return(-1*sum(logLiks)) # Overall negative log-likelihood - sum across trials
}

