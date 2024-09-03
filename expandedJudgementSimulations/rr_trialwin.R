#a function that generates a number of simulated "trials" for particular values of intercept and gradient, calculates the average reward rate and saves trial-level information
#arguments:  window size;
# the number of samples after which a trial "times out" (i.e., a within-trial response deadline);
# current intercept value;
# current gradient value; 
# gaussian noise around the current state of accumulated evidence; 
# drift rate of the task; 
# inter-trial interval for correct responses; 
# inter-trial interval for incorrect responses; 
# monetary penalty for incorrect responses; 
# a dummy variable indicating whether the information about individual trials should be returned (1) or not (0)
# pre-produced evidence samples sequence (set to 0 if the samples can vary for each combination of intercept*gradient)

rr_trialwin <- function(window=10, MaxTrialsamples=50, intercept = 10, gradient = 0, evinoise = 0, drift = 0.2, iti_cor = 15, iti_pen = 50, monpen = 0, returnTrials = 1,givensamples=0){
  evidence <- matrix(0,nrow=window,ncol=MaxTrialsamples)
  evidencenoisy <- matrix(0,nrow=window,ncol=MaxTrialsamples)
  time <- seq(1,MaxTrialsamples)
  action <- numeric(MaxTrialsamples)
  boundary <- numeric(MaxTrialsamples)
  
  itis <- rep(iti_pen,window) #assume ITIs for an incorrect response by default; overwritten later on in case of correct responses
  rewards <- rep(monpen,window) #assume a penalty (incorrect response) by default; overwritten later on in case of correct responses
  
  for (i in 1:window) {
    
    #calculate the value of the response boundary at each time step
    for (k in 1:MaxTrialsamples){
      boundary[k] <- intercept+k*tan(gradient) #generate boundary values for the current intercept and gradient
    }
    boundary[boundary<0] <- 0 #to ensure the threshold collapses to 0 and does not go to negative values
    
    if (length(givensamples)==1){ #if no evidence samples were pre-generated, generate them on the fly
      evisamples <- sample(x=c(1, -1),size=MaxTrialsamples,replace=TRUE, prob=c(0.5+drift, 0.5-drift))
    } else { #accept the pre-generated samples, if provided
      evisamples <- givensamples[i,]
    }
    evidence[i,] <- cumsum(evisamples) #transform the raw evidence samples into an evidence sequence
    evidencenoisy[i,] <- evidence[i,] + rnorm(n=MaxTrialsamples,mean=0,sd=evinoise) #apply representational Gaussian noise to the evidence sequence at each time-step independently
    p <-  which(abs(evidencenoisy[i,])>=boundary)[1] #the time point when evidence crossed a threshold
    
    #the following loop simply determines whether the response was correct/incorrect, assigns the corresponding reward/penalty and then creates all the relevant information about the trial into sample-by-sample vectors
    if (is.na(p)==TRUE){ # if p>MaxTrialsamples; the decision 'times out' and is coded as incorrect
      if (i == 1){ # first trial
       e <- evidence[i,1:(MaxTrialsamples)]
       en <- evidencenoisy[i,1:(MaxTrialsamples)]
       t <- time
       a <- action
       thresholdval <- boundary
       int <-  rep(intercept,MaxTrialsamples)
       grad <-  rep(gradient,MaxTrialsamples)
       trial <- rep(i,MaxTrialsamples)
      } else { # subsequent trials
       e <- c(e,evidence[i,1:(MaxTrialsamples)])
       en <- c(en,evidencenoisy[i,1:(MaxTrialsamples)])
       t <- c(t,time)
       a <- c(a,action)
       thresholdval <- c(thresholdval,boundary)
       int <-  c(int,rep(intercept,MaxTrialsamples))
       grad <-  c(grad,rep(gradient,MaxTrialsamples))
       trial <- c(trial,rep(i,MaxTrialsamples))
      }
    } else { # if p is valid (not a timed-out trial)
      if (i == 1){ #first trial
        if (evidencenoisy[i,p]>0){
          evisign <- 1 #positive value = correct response
        } else if (evidencenoisy[i,p]==0){
          evisign <- sample(x=c(1,-1), size=1)
        } else if (evidencenoisy[i,p]<0){ # flip the sign of evidence if an incorrect/negative response was made
          evisign <- -1 #negative value = incorrect response
        }
        e <- evidence[i,1:(p)]*evisign
        en <- evidencenoisy[i,1:(p)]*evisign
        t <- time(1:p)
        a <- action[1:p]
        a[length(a)] <- 1
        thresholdval <- boundary[1:p]
        int <-  rep(intercept,p)
        grad <-  rep(gradient,p)
        trial <- rep(i,p)
      } else { #subsequent trials
        if (evidencenoisy[i,p]>0){
          evisign <- 1 #positive value = correct response
        } else if (evidencenoisy[i,p]==0){
          evisign <- sample(x=c(1,-1), size=1)
        } else if (evidencenoisy[i,p]<0){ # flip the sign of evidence if an incorrect/negative response was made
          evisign <- -1 #negative value = incorrect response
        }
        e <- c(e,(evidence[i,1:(p)])*evisign)
        en <- c(en,(evidencenoisy[i,1:(p)])*evisign)
        t <- c(t,time[1:p])
        a <- c(a,action[1:p])
        a[length(a)] <- 1
        thresholdval <- c(thresholdval,boundary[1:p])
        int <-  c(int,rep(intercept,p))
        grad <-  c(grad,rep(gradient,p))
        trial <- c(trial,rep(i,p))
      }
      if(evisign==1){ #overwrite the rewards/penalties in case the response was correct (positive)
        rewards[i] <- 1
        itis[i] <- iti_cor
      }
    }
  }
  #put together data from all the simulated trials
  time <- t
  evi <- e
  evinoisy <- en
  action <- a
  data_window <- data.frame(trial,action,time,evi,evinoisy,int,grad,thresholdval)
  
  rr_mean <- sum(rewards)/(length(data_window$trial) + sum(itis)) #compute the average reward rate over the simulated trials
    
  if (returnTrials!=1){ #only return the average reward rate
    return(rr_mean=rr_mean)
  } else { #return the average reward rate, the total time (in samples) spent on these simulations, and the full dataframe with trial-level information
    return(list(rr_mean=rr_mean, nsampleswin=length(t), data_window=data_window)) # Return the objective value at the (x,y) position(s)
  }
}
