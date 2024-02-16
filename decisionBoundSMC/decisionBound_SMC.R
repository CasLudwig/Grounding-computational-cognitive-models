# C.J.H. Ludwig, January/February 2024
# Particle filtering approach for estimating time-varying parameters in an expanded judgement paradigm (cf. Malhotra et al., 2017). Here the data are generated through simulation and we are looking at identifying time-varying parameters.

library(ggplot2)
library(dplyr)
library(tmvtnorm) # For truncated multivariate normal

rm(list=ls())

# Navigate to the appropriate folder and load the function that computes the log-likelihood for a given (stimulus, action) sequence.
dirName <- "~/OneDrive - University of Bristol/Data/localSamplingComparison/decisionBoundApplication"
setwd(dirName)
source("likelihoodFunctions.R") # Functions for computing likelihood and doing MCMC fit.
source("bootstrapPF.R") # Functions for bootstrap particle filter sequential estimation of parameters with each new trial.

fName <- "data/toyDynamicData.csv"
allData <- read.csv(fName)

### Particle filter to identify parameter dynamics ###

# Set some of the hyperparameters of the algorithm: parameters of the truncated multivariate Gaussian prior and transition distributions. We specify the prior means (mu0), the prior variances (sigma0) and the variances of the transition distribution (sigmaQ). Note that we do not expect the decision noise to vary over time very much, as we assume this is not a parameter that is under (much) strategic control.
hParms <- list(mu0 = c(5, 0, 1), # prior means 
                   sigma0 = matrix(c(3^2, 0, 0, 0, 0.2^2, 0, 0, 0, 1^2), ncol = 3), # 3 x 3 covariance matrix; note sds on the diagonal are squared to give variance.
                   sigmaQ = matrix(c(2^2, 0, 0, 0, 0.2^2, 0, 0, 0, 0.5^2), ncol = 3), # 3 x 3 covariance matrix; note sds on the diagonal are squared to give variance.
                   lims = matrix(c(0.1, 20, -40 * pi/180, 10 * pi/180, 0.01, 5), ncol = 3) # 2 x 2 matrix with lower and upper bounds in the two rows.
)
nPart <- 2000

# Bootstrap particle filter parameter estimation for dataset
pfOut <- pfRun(data = allData, n = nPart, hyperParms = hParms)

# Match the true data generating parameters for each trial with the inferred states.
pfSummary <- summaryStats(pfOut) # Get useful summary statistics for each trial
allData$noise <- 2.5 # Add column with the (fixed) noise parameter used to generate the data
parmDF <- allData %>%
  filter(., time == 1) %>% # Select the data generating parameters for each trial
  select(., participant, trialID, int, grad, noise) %>% # We don't need various columns
  rename(., alpha = int, beta = grad, eta = noise) %>% # Rename generating parameters
  tidyr::pivot_longer(., cols = c("alpha", "beta", "eta"), 
                      names_to = "parameter", 
                      values_to = "trueVal") %>% # From wide to long
  left_join(., pfSummary, by=c("trialID", "parameter")) # Join the two dataframes, based on trialID and parameter columns

# Get the parameters estimated from a static model (loaded from a different datafile).
attach('data/staticParmsStan.RData'); 
stanParmSummary <- as.data.frame(fitSummary)
stanSamples <- dfPostSamples
detach('file:data/staticParmsStan.RData')

# Manipulate the posterior samples and summary statistics to get in a useful format for plotting
# For stanParmSummary we need to do some renaming to facilitate subsequent plotting
stanParmSummary <- select(stanParmSummary, c("mean", "2.5%", "25%", "50%", "75%", "97.5%"))
colnames(stanParmSummary) <- c("mean", "q_025", "q_25", "q_50", "q_75", "q_975")
stanParmSummary$parameter <- rownames(stanParmSummary)
stanParmSummary$trialID <- max(allData$trialID) + 10 # add a "virtual" trial ID
stanSamples$trialID <- max(allData$trialID) + 10 # add a "virtual" trial ID
stanSamples <- tidyr::pivot_longer(stanSamples,
                                   cols = c("alpha", "beta", "eta"), 
                                   names_to = "parameter", 
                                   values_to = "sampleVal")

# Visualise the comparison
cbbPalette <- c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour-deficient-friendly palette; first four entries: black, orange, blue, green
ggplot(data = parmDF) +
  geom_line(aes(x = trialID, y = trueVal), size = 2, colour = cbbPalette[1], show.legend = TRUE) +
  geom_ribbon(aes(x = trialID, ymin = q_025, ymax = q_975), colour = cbbPalette[2], alpha = 0.25) + # 95% credible interval
  geom_line(aes(x = trialID, y = post_mean), colour = cbbPalette[2], size = 2, show.legend = TRUE) + # Posterior mean
  # To the right, show the estimates from a static model (posterior parameters from Stan)
  geom_violin(data = stanSamples, aes(x = trialID, y = sampleVal), trim = FALSE, scale = "area", width = 10, fill = cbbPalette[3], colour = cbbPalette[3], alpha = 0.25)+
  geom_errorbar(data = stanParmSummary, aes(x = trialID, ymin = q_25, ymax = q_75), size = 1, colour = cbbPalette[1])+ # IQR
  geom_point(data = stanParmSummary, aes(x = trialID, y = mean), size = 2, shape = 15, colour= cbbPalette[3],position = position_dodge(width = 1), alpha = 0.75)+
  facet_grid(parameter~., scales = "free_y") +
  scale_x_continuous(breaks = c(0, 25, 50, 75), labels = c(0, 25, 50, 75)) +
  labs(x = 'trial', y='parameter value') +
  theme_bw()+
  theme(strip.text=element_text(size=14,face="bold"))+
  theme(axis.text.x  = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.title.x.bottom = element_text(size=14,face="bold"))+
  theme(axis.title.y.left = element_text(size=14,face="bold"))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=14, face="bold"))

  

### Single level MCMC for entire dataset with BayesianTools ###

# Parameter estimation for the overall dataset, using BayesianTools to perform MCMC sampling. Plot mean and CI of posterior densities as as a single data point + CI to the right of the dynamic estimates (i.e. at nTrials + 10 or something like that).

#library(BayesianTools) # For MCMC
#library(coda)
#library(bayesplot)

# data <- allData

# Specify prior parameters
# mu0 <- c(5, 0, 1) # means 
# sigma0 <- matrix(c(3^2, 0, 0, 0, 0.2^2, 0, 0, 0, 1^2), ncol = 3) # 3 x 3 covariance matrix; note sds on the diagonal are squared to give variance.
# lims <- matrix(c(0.1, 20, -40 * pi/180, 10 * pi/180, 0.01, 5), ncol = 3)
# prior <- createTruncatedNormalPrior(mean = mu0, sd = sqrt(diag(sigma0)), lower = lims[1,], upper = lims[2,])

# Specify parameters of MCMC sampler
# iter <- 30000
# burnin <- 5000
# nChains <- 3
# nCores<-nChains
# startVals <- prior$sampler(n = nChains)
# parmNames <- c("alpha", "beta", "eta")

# bayesSetup <- createBayesianSetup(likelihood = logLikOverall,
                                  prior = prior,
                                  parallel = nChains,
                                  names=parmNames)

# bayesSettings = list(iterations = iter, burnin = burnin) # startValue = startVals, #eps = 1^-12, e = 0.05; those are settings described in published papers on DREAM(zs), but have just stuck with BT defaults here.

# Perform sampling
# mcmcOut <- runMCMC(bayesianSetup = bayesSetup,
#                    sampler = "DREAMzs",
#                    settings = bayesSettings)

# Generate some diagnostics
# summary(mcmcOut$chain)
# print(ess <- lapply(mcmcOut$chain,effectiveSize)) # Effective sample size/chain; this can give some idea of the thinning required
# gelman.diag(mcmcOut$chain, multivariate = FALSE, autoburnin = FALSE)
# mcmc_trace(mcmcOut$chain, pars=parmNames)
# correlationPlot(mcmcOut$chain, whichParameters = parmNames)
# mcmc_acf(mcmcOut$chain, pars=parmNames)
# mcmc_dens(mcmcOut$chain, pars=parmNames)