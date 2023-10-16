# Function to calculate Lotka R, DK revised this with several options of how we want to calculate

# Function Arguments
# nat.mort:    Natural mortality, 
# fecund:      Fecundity
# ages:

simple.lotka.r<-function(nat.mort = NULL,fecund = NULL,ages = NULL) 
{
  
  # The Euler-Lotka Function to optimize on to find the intrinsic rate of growth for your population
  eulerlotka <- function(r) (sum(lx * mx * exp(-r * ages)) - 1)^2
  # Set up a dataframe
  res <- data.frame(r = NA)
  
  
  
  #################################### Calculate Lotka SECTION #################################### Calculate Lotka SECTION #################################### Calculate Lotka SECTION
  mx.tmp <- NULL
  lx.tmp <- NULL
  
  # Now run this through all the yrs 

    # Get the r estimate for each year if we are running through a bunch of years
    # Just use the first year if we don't have same length nat.mort as years
    if(!is.null(nrow(nat.mort))) lx <- 1-exp(-nat.mort) else lx <- 1-exp(-nat.mort)
    # Convert to survivorship vector
    si <- 1-lx 
    # Set the first age class to be 1
    lx<-1
    # And get cumulative survivorship for the stock
    for(s in 2:(length(si))) lx[s]<-lx[s-1]*si[s-1]
    mx <- fecund
    # Now we are cooking!
    junk<-optimize(lower=-0.999999,upper=15,f = eulerlotka)
    res <- c(junk$minimum)

  return(list(res=res,mx=mx,lx=lx,nat.mort = nat.mort,ages=ages))
}