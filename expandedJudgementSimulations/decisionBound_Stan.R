# C.J.H. Ludwig, January 2024
# Model for the expanded judgement task of Malhotra et al. (2017). We estimate parameters for simulated data. The parameters change over time, but here we just estimate a single set of fixed parameters for the overall dataset (i.e. ignoring the temporal variation).
# The model assumes perfect integration of evidence to a decision boundary that can vary in height and slope. After each sample, a probabilistic decision is made whether the 'wait' or 'go' (Gaussian decision noise).
# Parameter sampling done with Stan.

library(rstan)
library(coda)
library(bayesplot)

fitStanModel <- function(allData = NULL, nTrials = NULL, accuracy = NULL){
  # For Stan:
  ## For execution on a local, multicore CPU with excess RAM call:
  options(mc.cores = parallel::detectCores())
  ## To avoid recompilation of unchanged Stan programs, call:
  rstan_options(auto_write = TRUE, javascript = FALSE)
  
  # Put the data in a list for Stan
  allDataStan <- list(
    N = nTrials, # total number of trials (after rejecting any time outs)
    accuracy = accuracy, # accuracy of the final action of each trial
    nTotalSamples = nrow(allData), 
    trial = allData$trial,
    timeInTrial = allData$time, 
    cumEvidence = allData$evi,
    actions = allData$action
  )
  
  # Fit model
  staticModel <- stan_model(file = "decisionBoundStan.stan")
  staticModelStan <- sampling(staticModel,
                              data = allDataStan,  # named list of data
                              chains = 4,             # number of Markov chains
                              warmup = 1000,          # number of warm-up iterations per chain
                              iter = 6000,            # total number of iterations per chain
                              cores = 4,              # number of cores (could use one per chain)
                              thin = 5,
                              refresh = 1             # progress shown
                              # control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 15)
  )
  
  # Let's check it out
  # color_scheme_set("brightblue")
  # print(fitSummary <- summary(staticModelStan, pars = c("alpha", "beta", "eta"))$summary)
  # dfPostSamples <- as.array(staticModelStan, pars = c("alpha","beta","eta")) # Extract samples
  # mcmc_pairs(dfPostSamples)
  # mcmc_trace(dfPostSamples)
  # print(R <- rhat(staticModelStan, pars = c("alpha","beta","eta")))
  # mcmc_rhat(R) + yaxis_text(hjust = 1)
  # mcmc_neff(neff_ratio(staticModelStan), size = 2)
  # mcmc_acf(dfPostSamples)
  # mcmc_dens(
  #   dfPostSamples, 
  #   pars = c("alpha", "beta", "eta"),
  #   facet_args = list(nrow=1,ncol=3,scales="free")
  # )
  
  fitSummary <- summary(staticModelStan, pars = c("alpha", "beta", "eta"))$summary
  dfPostSamples <- as.data.frame(staticModelStan, pars = c("alpha","beta","eta")) # Extract samples again, but put them in a dataframe, as this is easier to work with
  
  return(list(dfPostSamples = dfPostSamples, fitSummary = fitSummary))
         
}

### OPTIONAL: CHECK OUT PRIORS ###
# priorModel <- stan_model(file = "stanPriorSamples.stan")
# priorVals <- sampling(priorModel,
#                       chains = 1,
#                       warmup = 1000,
#                       iter = 6000,
#                       cores = 4,
#                       refresh = 1,
#                       control = list(adapt_delta = 0.9, stepsize = 0.5, max_treedepth = 20)
#                       )
# print(priorSummary <- summary(priorVals, pars = c("alpha", "beta", "eta"))$summary)
# dfPriorSamples <- as.array(priorVals, pars = c("alpha","beta","eta")) # Extract samples
# mcmc_dens(
#   dfPriorSamples,
#   pars = c("alpha", "beta", "eta"),
#   facet_args = list(nrow=1,ncol=3,scales="free")
# )
# save.image("data/staticParmsStan.RData")


# For comparison, let's do a MLE fit
# library(optimx)
# source("likelihoodFunctions.R")
# nStartVals <- 10
# startVals <- matrix(c(runif(nStartVals, 0.1, 20),
#                     runif(nStartVals, -0.7, 0.175),
#                     runif(nStartVals, 0.01, 5)), nrow = nStartVals)
# mleFit <- multistart(startVals, 
#                      logLikOverall,
#                      lower = c(0.1, -0.7, 0.01),
#                      upper = c(20, 0.175, 5),
#                      method = "L-BFGS-B",
#                      data = allData)
# mleFit <- mleFit[order(mleFit$value),]
# bestParms <- mleFit[1,1:3]
# colnames(bestParms) <- c("alpha", "beta", "eta")