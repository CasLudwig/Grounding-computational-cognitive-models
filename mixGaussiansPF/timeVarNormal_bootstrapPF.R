# Suite of functions to conduct particle filter estimation of dynamic parameter trajectories. Some functions copied/modified from Speekenbrink (2016, JMathPsych; Supplementary data in Appendix A).

# Code here is tailored for a (toy) problem where we estimate the mean and sd of a Gaussian distribution, where both of these parameters change over time.

library(Hmisc) # For weighted summary stats
library(tmvtnorm) # For truncated multivariate normal

ESS <- function(weights) {
  # Effective sample size
  # weights: a vector with normalised weights.
  # output: effective sample size
  1/sum(weights^2)
}


resample_systematic <- function(weights) {
  # Systematic resampling
  # weights: a vector of length N with (unnormalised) importance weights
  # output: a vector of length N with indices of the replicated particles
  N <- length(weights)
  weights <- weights/sum(weights) # normalize weights
  csum <- cumsum(weights)
  u1 <- runif(1,min=0,max=1/N) # draw a single uniform number
  u <- (1:N - 1)/N + u1
  idx <- vector("integer",length=length(weights)) # allocate a vector for the results
  j <- 1
  for(i in 1:N) {
    while (u[i] > csum[j]) {
      j <- j + 1
    }
    idx[i] <- j
  }
  return(idx)
}


initParticles <- function(nPart = 2000, mu, sigma, lims){
  # Initialise particles and weights.
  # Parameters m and s drawn from a truncated multivariate Gaussian.
  # Return a list with the particle values and the weights.
  particles <- rtmvnorm(n = nPart,
                        mean = mu, 
                        sigma = sigma,
                        lower = lims[1,],
                        upper = lims[2,])
  weights <- rep(1/nPart, nPart)
  return(list(particles = particles, weights = weights))
}


propParticles <- function(nPart = 2000, preVals, sigma, lims){
  # Propagate the particles and the weights, given some transition distribution.
  # Parameters m and s drawn from a truncated multivariate Gaussian.
  # preVals is a nPart x 2 matrix and for each row of parameters, we will draw a new pair centred on the previous values.
  newVals <- matrix(rep(0,nPart * ncol(sigma)), ncol = ncol(sigma))
  for (i in 1:nPart){
    newVals[i,] <- rtmvnorm(1, mean = preVals[i,], sigma = sigma, lower = lims[1,], upper = lims[2,])
  }
  return(newVals)
}


# Takes the list output from pfRun and returns some useful summary statistics of the posterior distributions.
summaryStats <- function(pf){
  nTrials <- dim(pf$particles)[1]-1 # number of trials
  nPart <- dim(pf$particles)[2] # number of particles
  nVars <- dim(pf$particles)[3] # number of parameters
  
  # Set up a dataframe to hold the summary stats
  pfSummary <- data.frame(trialID = rep(1:nTrials, each = nVars),
                         parameter = factor(rep(c("mu", "sigma"), times = nTrials)),
                         q_025 = rep(0, nTrials * nVars),
                         q_05 = rep(0, nTrials * nVars),
                         q_25 = rep(0, nTrials * nVars),
                         q_5 = rep(0, nTrials * nVars),
                         q_75 = rep(0, nTrials * nVars),
                         q_95 = rep(0, nTrials * nVars),
                         q_975 = rep(0, nTrials * nVars),
                         post_mean = rep(0, nTrials * nVars))
  rowid <- 1 # pfSummary row index
  
  # Loop over trials and extract the relevant measures
  for(i in 1:nTrials){
    tmpPart <- pf$particles[i+1,,] # Note +1, because the trial dimension contains the initial set of particles (drawn from prior)
    tmpW <- pf$weights[i+1,]
    
    # Get the various statistics of interest
    # Get some useful quantiles (using apply to loop over parameters)
    tmpQuants <- apply(tmpPart,2,wtd.quantile, weights = tmpW * nPart, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)) 
    tmpMean <- apply(tmpPart,2,wtd.mean, weights = tmpW * nPart) # posterior mean

    # Write to dataframe
    rids <- c(rowid, rowid+nVars-1)
    pfSummary[rids[1]:rids[2], 3:ncol(pfSummary)] <- cbind(t(tmpQuants), tmpMean)

    # update row index for the next trial
    rowid <- rowid + nVars
  }
  pfSummary$logLik <- rep(pf$logLiks[2:(nTrials+1)], each = nVars)
  return(pfSummary)
}


# Function that runs a bootstrap particle filter for parameter estimation on a single dataset.
# Prior (hyper) parameters are set within this function, but you can supply these as a list if you want to. By default, we use n = 2000 particles.
pfRun <- function(data, n = 2000, hyperParms = NULL){

  trialIDs <- unique(data$trialID)
  nTrials <- length(trialIDs)
  resampleCrit <- 0.5
  
  # If hyper-parameters are not specified, we set them here. 
  # The prior for (m, s) is a truncated multivariate Gaussian, so we need to set the means, the variance-covariance matrix and the truncation limits. Gradients are in radians. 
  # The list also includes the variance-covariance of the transition distribution Q. We assume a multivariate Gaussian for the transitions, but with much smaller variances, so that the jumps are generally quite small.
  if(is.null(hyperParms)){
    hyperParms <- list(mu0 = c(0, 1), # prior means 
                       sigma0 = matrix(c(1^2, 0, 0, 0.5^2), ncol = 2), # 2 x 2 covariance matrix for initial particle values; note sds on the diagonal are squared to give variance.
                       sigmaQ = matrix(c(0.25^2, 0, 0, 0.1^2), ncol = 2), # 2 x 2 covariance matrix for transition distribution; note sds on the diagonal are squared to give variance.
                       lims = matrix(c(-3, 3, 1e-6, 5), ncol = 2) # 2 x 2 matrix with lower and upper bounds in the two rows.
                      )
  }
  
  # Initialise the particles and weights.
  pf <- initParticles(nPart = n, mu = hyperParms$mu0, sigma = hyperParms$sigma0, lims = hyperParms$lims)
  
  # Arrays to store the particles, weights, and log-likelihoods. Note we use nTrials+1 to store the initial sets at time 0 (i.e. drawn from prior).
  particleArray <- array(0.0, dim=c(nTrials+1, n, ncol(pf$particles)))
  weightArray <- matrix(1/n, nrow=nTrials+1, ncol=n)
  logLiks <- rep(NA,nTrials+1)
  
  # Store initial set of particles and weights (drawn from prior)
  particleArray[1,,] <- pf$particles
  weightArray[1,] <- pf$weights
  
  # Trial-wise inference
  for(i in 1:nTrials){

    print(i)
    # Propagate the particles through the transition distribution.
    pf$particles <- propParticles(nPart = n, 
                                  preVals = pf$particles, 
                                  sigma = hyperParms$sigmaQ, 
                                  lims = hyperParms$lims)
    
    # Loop over the matrix of particles and get the log-likelihood for each particle.
    ll <- apply(pf$particles, 1, logLikOneTrial, obs = data$y[i], returnLogLik = TRUE)
    likelihoods <- exp(ll) # go from log-likelihood to likelihood
    # Compute expectation of the log-likelihood: log of the weighted average likelihood.
    logLiks[i+1] <- log(sum(likelihoods*pf$weights))
    
    # Update weights, so that they are proportional to the likelihoods (if the weights were uniform to begin with).
    pf$weights <- pf$weights*likelihoods
    pf$weights <- pf$weights/sum(pf$weights) # normalise to 1
    
    # Resample if there are too many low weights (weight degeneracy)
    if(ESS(pf$weights) < resampleCrit*n) {
      # Draw the indices of the replicated particles with systematic resampling
      idx <- resample_systematic(pf$weights)
      # Select the replicated particles
      pf$particles <- pf$particles[idx,]
      # Uniform weights after resampling
      pf$weights <- rep(1/n,n)
    }
    
    # Store particles and their weights
    particleArray[i+1,,1:2] <- pf$particles
    weightArray[i+1,] <- pf$weights
  }
  return(list(particles = particleArray, weights = weightArray, logLiks = logLiks))
}
