# this function uses the parameter estimates to get an r max value for a population which is then used to set an upper ceiling on r (in the RUNME code) , you also get samples ready from this for MCMC forward projections
# The sampling as is is very specific, gets us natural mortality and juvenile mortality (similar to the parms.calc method from Porbeagle)

# age.mat:     Age at maturity from params_calc function
# n.offspring: Number of offspring a female produces on average, for fish we'll need this for each mature age class, so should be a vector
# gest:        The gestation period, for fish well make this 1 (i.e. this will effectively enable an individual to reproduce once a year
# nat.mort:    Natural mortality, probably will be 0.2 for most stocks
# juv.mult:    A multiplier to allow for juvenile mortality if that's of interest.
# max.age:     Maximum age
# repo.lag:    Lag in years between giving birth and conceiving again.  Likely 0 for fish.
# u:           Exploitation rate (annual not instantaneous), currently set up to be 1 value.
# sel:         Selectivity, currently set up to be 1 value
#q.ext           quasi-extinction threshold for forward projections
#f.sigma         lognormal shape parameter for forward projections
#f.auto          autocorrelation parameter for forward projections
#forward.n.years project for 30 years (i.e. possible increase in Curtis)
#n.sims          SET NUMBER OF SIMS

parms.short<-function( age.mat=9,n.offspring =4, gest =2, nat.mort=0.112, juv.mult=1, 
                       max.age=40,lag.f=0, u=0,sel=2,n.sims = 1000,q.ext=20, f.sigma=0.01,f.auto = 0.03,forward.n.year = 30)
{
  # deterministic Rmax
  det.r<-lotka.r(age.mat,n.offspring,gest,nat.mort,juv.mult,max.age,lag.f,u,sel)
  print("determinisitc rmax")
  print(det.r$par)
  # The rest of this code is very white shark specific, I'm not sure how we'd generalize this?
  r.cutoff<-2*det.r$par 
  age.mat.vec<-sample(8:11,n.sims*2,replace=T) #9.5 - note that proportionately, this is about the same level of variability (16% above, 16% below estimate - same as in short scenario)
  litter.vec<-sample(2:6,n.sims*2,replace=T)    #8 total switched to females =4
  gestation.vec<-sample(2:3,n.sims*2,replace=T) #2.5
  nat.mort<-seq(0.062,0.162,by=0.001)
  nat.mort.vec<-sample(nat.mort,n.sims*2,replace = T) #0.112+0.05
  juv.mult.vec<-rep(1,n.sims*2) #0
  lag.vec<-rep(0,n.sims*2) #0
  max.age.vec<-sample(35:45,n.sims*2,replace=T) #30
  N.end<-c(sample(2000:2499,size=n.sims/2,replace=T),sample(2500:3000,size=n.sims/2,replace=T)) # see Burgess et al. 2014
}