#this function finds the value of f (of the vulnerable, not total, population taken each year) note, because the age of selectivity is so low, this is similar to the total.
#this is based on an equilibrium age structure assumption selectivity is knife-edged at sel and 0 after maturation this calculates u only for the known removals component

# age.mat   The age at maturity estimated in parms_calc using life history information
# max.age   The longevity of the species estimated in parms_calc using life history information
# nat.mort:    Natural mortality, probably will be 0.2 for most stocks
# juv.mult:    A multiplier to allow for juvenile mortality if that's of interest.
# sel       The selectivity of the stock.  Currently this is set up as a single number and it is tweaked in the code, we'll need to make this more complex.
# removals  Removals from the fishery in a given year.
# N         Population size in a given year


u.calc<-function(age.mat,max.age,nat.mort,juv.mult,sel,removals,N)
{
  # Function to minimze
  minimise<-function(f)
  {
    # The usual set up to get the survival for each age, this accounts for natural mortality and fishing mortality
    f.vec<-c(rep(0,sel+1),rep(f,age.mat-sel),rep(0,max.age-age.mat))
    # Survivorship
    si<-exp(-(M+f.vec))
    # Not sure
    lx.f<-0
    #instanateous exploitation rate
    u<-1-exp(-f)
    #lx is the proportion alive at the start of each age class
    lx.f[1]<-1
    for(i in 2:(length(si)+1)) lx.f[i]<-lx.f[i-1]*si[i-1]

    #Add in number below the age of selectivity, calculated using lx,
    #This is the number of individuals in the population that are vulnerable to the fishery
    vulnerable<-N*(sum(lx.f[(.sel+1):length(lx.f)])/sum(lx.f))
    #   print(vulnerable)
    # Here is the catch given the exploitation rate
    pred.C<-vulnerable*u

    #   print(pred.C)
    # Now we need to get removals to equal the predicted removals so we can get our estimate of the exploitation rate
    sumsq <- (removals-pred.C)^2
    return(sumsq)
  }

  # Get natural mortality vector
  
  M<-c(nat.mort*juv.mult,rep(nat.mort,max.age-1))

  # Solve to minimize difference and get an estimate of the exploitation rate
  junk<-nlminb(start = 0.1,obj = minimise, lower=0.0001,upper=10)
  # Extract the exploitation rate.
  f.out<-junk$par
  u.out<-1-exp(-f.out) # Proportion removed

  return(u.out)
}