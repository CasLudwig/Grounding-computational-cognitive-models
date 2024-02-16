// C.J.H. Ludwig, January 2024

// Model for the expanded judgement task of Malhotra et al. (2017). We estimate parameters for simulated data. The parameters change over time, but here we just estimate a single set of fixed parameters for the overall dataset (i.e. ignoring the temporal variation).
// The model assumes perfect integration of evidence to a decision boundary that can vary in height and slope. After each sample, a probabilistic decision is made whether the 'wait' or 'go' (Gaussian decision noise).

// #include /pre/license.stan

data {
  int<lower=1> N; // total number of trials
  int<lower=0, upper=1> accuracy[N]; // vector indicating whether the final 'go' action on a trial was correct (1) or not (0)
  int<lower=1> nTotalSamples; // total number of samples taken across all trials
  int<lower=1> trialID[nTotalSamples]; // vector with trial numbers
  int<lower=0> timeInTrial[nTotalSamples]; // sample number within a trial
  int cumEvidence[nTotalSamples]; // cumulative evidence within a trial
  int <lower=0, upper=1> actions[nTotalSamples]; // wait (0) or go (1) action
}

parameters {
  // Single-level parameters, specified as normal random deviates.
  real alpha_norm; // intercept
  real beta_norm;  // gradient
  real eta_norm; // decision noise
}

transformed parameters {
  // Parameters are specified (above) as normal random deviates (the single-level '._norm' variables).
  // Pass these through a cumulative normal to put on a scale from 0-1. Then rescale to put
  // on the desired, bounded scale. *These* are the parameters of interest.

  // single-level parameters
  real<lower=0.1, upper=20> alpha;
  real<lower=-0.7, upper=0.175> beta;
  real<lower=0.01, upper=5> eta;

  // Here we transform the single-level parameters to their "natural", meaningful and bounded scale.
  // Make sure that the bounds here match those in the declaration of the parameters above (and below, when specifying priors).
  alpha = 0.1 + Phi_approx(alpha_norm) * (20 - 0.1);
  beta = -0.7 + Phi_approx(beta_norm) * (0.175 + 0.7);
  eta = 0.01 + Phi_approx(eta_norm) * (5 - 0.01);
}

model {

  vector[3] mu; // approximate prior means for alpha, beta and eta
  vector[3] p; // probability of going up (correct), going down, waiting
  real logLik; // log-likelihood of each action
  real upBoundary;
  real loBoundary;

  // Approximate prior means for alpha, beta, eta
  mu[1] = 5;
  mu[2] = 0;
  mu[3] = 1;
  
  // Somewhat informed priors for the single-level parameters, approximately centred on mu. Because of the transformation onto a bounded scale above, these priors will not be normal, but have some degree of skew.
  // The is ugly, but basically, we compute the random normal deviate that corresponds to the specified mu value, given the limits we impose on the parameters (inverse cumulative normal on (mu - l)/(u - l), where l and u are the lower and upper limits of the parameters on their natural scale. The sds are chosen by trial and error to get close to values of (3, 0.2 and 1).
  // The resulting distributions do not have an obvious parametric description, but they are generally unimodal and skewed, but with a good amount of mass over a reasonable parameter range.
  target += normal_lpdf(alpha_norm | inv_Phi((mu[1]-0.1)/(20-0.1)), 0.5);
  target += normal_lpdf(beta_norm | inv_Phi((mu[2]+0.7)/(0.175+0.7)), 0.75);
  target += normal_lpdf(eta_norm | inv_Phi((mu[3]-0.01)/(5-0.01)), 0.75);

  // Loop over all samples: we do not need to differentiate between trials
  for (i in 1:nTotalSamples) {

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
      if (accuracy[trialID[i]] == 1){ // go, correct (up)
        logLik = p[1];
      }else{logLik = p[2];} // go, incorrect (down)
    }
    if(logLik <= 0){logLik = 1e-12;} // avoid numerical problems with p <= 0
    if(logLik >= 1){logLik = 1 - 1e12;} // avoid numerical problems with p >= 1
    
    target += log(logLik);
  }
}
