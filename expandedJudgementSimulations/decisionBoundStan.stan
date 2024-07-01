// C.J.H. Ludwig, January 2024

// Model for the expanded judgement task of Malhotra et al. (2017). 
// We estimate parameters for simulated data. 
// The parameters may change over time, but here we just estimate a single set 
// of fixed parameters for the overall dataset (i.e. ignoring the temporal variation).
// The model assumes perfect integration of evidence to a decision boundary that 
// can vary in height and slope. After each sample, a probabilistic decision is 
// made whether the 'wait' or 'go' (Gaussian decision noise).

// #include /pre/license.stan

data {
  int<lower=1> N; // total number of trials
  int<lower=0, upper=1> accuracy[N]; // vector indicating whether the final 'go' action on a trial was correct (1) or not (0)
  int<lower=1> nTotalSamples; // total number of samples taken across all trials
  int<lower=1> trial[nTotalSamples]; // vector with trial numbers
  int<lower=0> timeInTrial[nTotalSamples]; // sample number within a trial
  int cumEvidence[nTotalSamples]; // cumulative evidence within a trial
  int<lower=0, upper=1> actions[nTotalSamples]; // wait (0) or go (1) action
}

parameters {
  // Single-level parameters, specified as standard normal random deviates.
  real alpha_norm; // intercept
  real beta_norm;  // gradient
  real eta_norm; // decision noise
}

transformed parameters {
  // Parameters are specified (above) as normal random deviates (the single-level '._norm' variables).
  // They are transformed here in that they are put on their "natural", meaningful scale.

  // single-level (transformed) parameters - these are the ones of interest.
  real alpha;
  real beta;
  real eta;
  
  alpha = exp(alpha_norm + 1.5);
  beta = beta_norm*0.5;
  eta = exp(eta_norm + 1);
}

model {
  
  vector[3] p; // probability of going up (correct), going down, waiting
  real logLik; // log-likelihood of each action
  real upBoundary;
  real loBoundary;
  int trialIndex;
  int prevTrialIndex;

  // parameters from standard normal distribution
  target += normal_lpdf(alpha_norm | 0, 1);
  target += normal_lpdf(beta_norm | 0, 1);
  target += normal_lpdf(eta_norm | 0, 1);
  
  // Loop over all samples: we do not need to differentiate between trials
  trialIndex = 1;
  prevTrialIndex = trial[1];
  for (i in 1:nTotalSamples) {

    // Check whether we have moved on to a new trial or are still in the current trial. We need this to deal with any rejected trials (e.g. when trialID array jumps from n to n+2).
    if (trial[i] > prevTrialIndex){
      trialIndex = trialIndex + 1;
      prevTrialIndex = trial[i];
    }
    // Given the time within the trial, what are the lower and upper decision boundaries?
    upBoundary = alpha + timeInTrial[i] * tan(beta);
    if (upBoundary < 0){upBoundary = 0;}
    loBoundary = -1 * upBoundary;
    
    // Given current evidence sample, boundary and decision noise, compute probability of each of the three possible actions
    p[1] = 1 - normal_cdf(upBoundary, cumEvidence[i], eta); // probability of going up (correct if actions = 1, because in these simulations where the evidence was pointing up)
    p[2] = normal_cdf(loBoundary, cumEvidence[i], eta); // probability of going down
    p[3] = 1 - p[1] - p[2]; // probability of waiting
    
    // Compute the log-likelihood of the action
    if (actions[i] == 0){ // wait
      logLik = p[3];
    }else{
      if (accuracy[trialIndex] == 1){ // go, correct (up)
        logLik = p[1];
      }else{logLik = p[2];} // go, incorrect (down)
    }

    target += log(logLik);
  }
}

