# Function to calculate Lotka R, this is a simplified version of a more complex function (the more complex function allows for life history estimates of fecund and mort)
# In this case we are feeding it the mortality and fecundity directly, effectively for this function we know the components of the Leslie matrix and are just
# doing the optimization.

# Function Arguments
#
# mort:     Mortality for each age class, you could use natural mortality (M) or total mortality (Z=F+M).  The function is expecting natural mortality to be instantaneous
# fecund:   Fecundity by age class, this is the average number of offspring produced by a female in an age class
# ages:     The ages vector, ideally this goes from 1 (or even 0) to the maximum age, should have same number of age classes as the mortality and fecundity vectors(matrices)
 
simple.lotka.r<-function(yrs = 1,mort = NULL,fecund = NULL,ages = NULL) 
{  
  
  
  # The Euler-Lotka Function to optimize on to find the intrinsic rate of growth for your population
  eulerlotka <- function(r) (sum(lx * mx * exp(-r * ages)) - 1)^2
  # Set up a dataframe
  res <- data.frame(years = NA,r = NA)
  
  n.yrs <- length(yrs)
  mx.tmp <- NULL
  lx.tmp <- NULL
  
  # Now run this through all the yrs 
  for(j in 1:(n.yrs))
  { 
    #browser() 
    # Get the r estimate for each year if we are running through a bunch of years
    # Just use the first year if we don't have same length nat.mort as years
    if(!is.null(nrow(mort))) lx <- 1-exp(-mort[j,]) else lx <- 1-exp(-mort)
    # Convert to survivorship vector
    si <- 1-unlist(lx)
    # Set the first age class to be 1
    lx<-1
    # And get cumulative survivorship for the stock
    for(s in 2:(length(si))) lx[s]<-lx[s-1]*si[s-1]
    #browser()
    if(!is.null(nrow(fecund))) mx <- fecund[j,] else mx <- fecund
    # Now we are cooking!
    junk<-optimize(lower=-0.99,upper=5,f = eulerlotka)
    #browser()
    res[j,] <- c(yrs[j],junk$minimum)
    mx.tmp[[j]] <- mx
    lx.tmp[[j]] <- lx
    #browser()
  } # end for(j in 1:(n.yrs-1))
  
  mx.mat <- do.call("rbind",mx.tmp)
  lx.mat <- do.call('rbind',lx.tmp)
  
  return(list(res=res,mx=mx,lx=lx,mort = mort,ages=ages))
}
