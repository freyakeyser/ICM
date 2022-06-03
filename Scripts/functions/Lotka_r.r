# Function to calculate Lotka R

# Function Arguements
# age.mat:     Age at maturity from params_calc function
# n.offspring: Number of offspring a female produces on average, for fish we'll need this for each mature age class, so should be a vector
# gest:        The gestation period, for fish well make this 1 (i.e. this will effectively enable an individual to reproduce once a year
# nat.mort:    Natural mortality, probably will be 0.2 for most stocks
# juv.mult:    A multiplier to allow for juvenile mortality if that's of interest.
# max.age:     Maximum age
# repo.lag:    Lag in years between giving birth and conceiving again.  Likely 0 for fish.
# u:           Exploitation rate (annual not instantaneous), currently set up to be 1 value.
# sel:         Selectivity, currently set up to be 1 value


lotka.r<-function(age.mat,n.offspring, gest = 1,nat.mort,juv.mult =1 ,max.age,repo.lag=0,u,sel)
{
 
  # Function to solve for lotka r.
  minimise<-function(start.r)
  {
    # Get a starting estimate of lotka's r
    rx<-exp(-start.r*x)
    # Multiply the repro contribution of each age by current value of rx and sum those up, as this approaches the correct value this approaches 1
    lotka<-sum(lxmx*rx)
    # Then minimize the sum of squares difference of this lotka term from 1 to determine what r is
    sumsq <- sum((lotka-1)^2)
    return(sumsq)
  }
  
  # The gestation period
  #gest.lag<-gest+repo.lag
  # A vector of something....
  #temp<-c(0,c(rep(0,age.mat+gest-1),n.offspring,rep(c(rep(0,gest.lag-1),offspring),100))) 
  # A fecundity vector 
  #mx<-temp[1:(max.age+1)]
  mx <- n.offspring
  r.M <- nat.mort
  #r.M<-c(nat.mort*juv.mult,rep(nat.mort,max.age-1))
  #browser()
  # Instantaneous fishing mortality
  eef<- -log(1-u)
  ## Set up the fishing mortality vector based on selectivity, max age, age at maturity, etch
  f.vec<-c(rep(0,sel),rep(eef,max.age-sel))      
  # Survival by age including natural mortality and fishing mortality
  si<-exp(-(r.M+f.vec))
  
  # Now calculate survivorship by age
  x<-0:(max.age)
  lx<-0
  lx[1]<-1
  for(i in 2:(length(si)+1)) lx[i]<-lx[i-1]*si[i-1]
  #browser()
  # Multiply survivorship by fecundity to get reproductive contribution of each age
  lxmx<-lx[1:length(si)]*mx

  # and now solve the lotka r minimizing function to get an estimate of reproductive capacity.
  junk<-nlminb(start = 0.03,obj = minimise, lower=-1,upper=1)
  return(junk)
}