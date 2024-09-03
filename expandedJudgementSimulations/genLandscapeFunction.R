#a function that produces a reward rate landscape, plots it and returns both the matrix with reward rate information, as well as the plot
#arguments: range of intercept values; 
# range of gradient values; 
# drift rate of the task; 
# inter-trial interval for correct responses; 
# inter-trial interval for incorrect responses; 
# monetary penalty for incorrect responses; 
# window size; 
# gaussian noise around the current state of accumulated evidence; 
# pre-produced evidence samples sequence (set to 0 if the samples can vary for each combination of intercept*gradient)
gen_landscape <- function(intrange,gradrange,drift,iti_cor,iti_pen,monpen,window,evinoise,givensamples=0){
  
  #prepare a sequence of intercept and gradient values; call another function to run simulations for each parameter combination and store the result into a matrix 
  intercepts = seq(intrange[1],intrange[2],by=1)
  gradients = seq(gradrange[2],gradrange[1],by=-pi/72)
  rr_obj <- matrix(0,ncol=length(gradients),nrow=length(intercepts))
  for(i in 1:length(intercepts)){
    for(j in 1:length(gradients)){
      rr_obj[i,j] <- rr_trialwin(window=window, MaxTrialsamples=MaxTrialsamples, intercept = intercepts[i], gradient = gradients[j], evinoise = evinoise, drift = drift, iti_cor = iti_cor, iti_pen = iti_pen, monpen = monpen,returnTrials = 0,givensamples)
    }
  }
  
  #"tidy up" the matrix into a data-frame, find out the maximum reward rate and create a plot
  rr_obj <- rr_obj[nrow(rr_obj):1,]
  rownames(rr_obj) <- intercepts[length(intercepts):1]
  colnames(rr_obj) <- round(180*gradients/pi,3)
  rr_landscape <- melt(rr_obj)
  colnames(rr_landscape) <- c("Intercept", "Gradient", "RR_estim")
  
  max_index <- which(rr_landscape$RR_estim==max(rr_landscape$RR_estim))#get the maximum RR_estim
  max_rr <- c(rr_landscape$Intercept[max_index],rr_landscape$Gradient[max_index]) #get the intercept and gradient associated with highest RR
  
  plot <- ggplot(rr_landscape, aes(x = Gradient, y = Intercept)) +
    geom_raster(aes(x= Gradient, y = Intercept, fill = RR_estim)) +
    geom_point(x=-max_rr[2],y=max_rr[1],fill='black',size=5, shape=17) +
    scale_fill_viridis(na.value="transparent", limits = c(0,max(rr_landscape$RR_estim)), breaks=c(0,0.01,0.02,0.03,0.04)) +
    scale_x_reverse(breaks = c(5,0,-10,-20,-30,-40,-50,-60)) +
    scale_y_continuous(breaks = c(0,5,10,15,20,25)) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  #return the output: 1 - plot, 2 - maximum reward rate value, 3 - the "raw" data-frame containing reward rates for each parameter combination
  return(list(landscape_plot = plot,max_rr,rr_landscape))
}