# C.J.H. Ludwig, July 2023
# Some functions for generating 2D objective functions with some arbitrary, but interesting profiles.

# Define 2D objective function with some (non-linear) correlation between the two parameters.
obj2D_uniModal <- function(x=NULL, y=NULL, mu=c(0,0), sigma=c(1,1), corrxy=0, amp=1){
  x <- x
  y <- log(y)
  z <- amp * exp(-(1/(2*(1-corrxy^2))) * (((x-mu[1])/sigma[1])^2 -
                                            (2*corrxy) * ((x-mu[1])/sigma[1]) * ((y-mu[2])/sigma[2]) +
                                            ((y-mu[2])/sigma[2])^2))
  if(any(is.na(z))){
    z[which(is.na(z))] <- 0
  }
  return(z) # Return the objective value at the (x,y) position(s)
}

# Almost the same function, but now we return the squared difference between the (x,y) objective value and the peak. This is useful when we want to find the location of the peak in an automated way (e.g. through optim).
obj2D_uniModal_sqErr <- function(xy=NULL, mu=c(0,0), sigma=c(1,1), corrxy=0, amp=1){
  x <- xy[1]
  y <- log(xy[2])
  z <- amp * exp(-(1/(2*(1-corrxy^2))) * (((x-mu[1])/sigma[1])^2 -
                                            (2*corrxy) * ((x-mu[1])/sigma[1]) * ((y-mu[2])/sigma[2]) +
                                            ((y-mu[2])/sigma[2])^2))
  if(any(is.na(z))){
    z[which(is.na(z))] <- 0
  }
  return((amp-z)^2) # Return the squared difference between the z at (x,y) and the peak of the objective function.
}

