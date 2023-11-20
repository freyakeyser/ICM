# Function to calculate Lotka R, DK revised this with several options of how we want to calculate

# Function Arguments
# Z:    Total mortality, 
# fecund:      Fecundity
# ages:

simple.lotka.r<-function(Z = NULL,fecund = NULL,ages = NULL) 
{  
  
  # The Euler-Lotka Function to optimize on to find the intrinsic rate of growth for your population
  eulerlotka <- function(r) (sum(lx * mx * exp(-r * ages)) - 1)^2
  # Set up a dataframe
  res <- data.frame(r = NA)
  
  
  
  #################################### Calculate Lotka SECTION #################################### Calculate Lotka SECTION #################################### Calculate Lotka SECTION
  # Now run this through all the yrs 

    # Get the r estimate for each year if we are running through a bunch of years
    # Just use the first year if we don't have same length Z as years
    if(!is.null(nrow(Z))) lx <- 1-exp(-Z) else lx <- 1-exp(-Z)
    # Convert to survivorship vector
    si <- 1-lx 
    # Set the first age class to be 1
    lx<-1
    # And get cumulative survivorship for the stock
    for(s in 2:(length(si))) lx[s]<-lx[s-1]*si[s-1]
    mx <- fecund
    # Now we are cooking!
    #browser()
    junk<-optimize(lower=-5,upper=5,f = eulerlotka)
    res <- c(junk$minimum)

  return(list(res=res,mx=mx,lx=lx,Z = Z,ages=ages))
}