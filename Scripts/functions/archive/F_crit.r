# Finds the critical value of F

#age.mat   The age at maturity estimated in parms_calc using life history information
#max.age   The longevity of the species estimated in parms_calc using life history information
#mx        The fecundity ogive, so this integrates age at maturity and weight information to sort out how many, estimated in parms_calc using life history information
#r.M       The natural mortality estimated in parms_calc using life history information
#sel       The selectivity of the stock.  Currently this is set up as
#removals  The average removals over the time series.

F.crit<-function(age.mat,max.age,mx,r.M,sel,removals)
{
  ## note that mx and r.M are outputs of the parms_calc function
  M<-r.M

  # Store in expression frame # DK note -- Heather can explain the advantage of this to me!
  assign(".M",M)
  assign(".mx",mx) #
  assign(".sel",sel)
  assign(".max.age",max.age)
  assign(".age.mat",age.mat)

  # The objective function for nlminb, mx is the fecundity ogive, so integrates size and maturity into one component, that's handy!
  minimize<-function(f)
  {
    # DK note -- Should max.age and age.mat be .max.age and .age.mat, I think so , but also think it doesn't matter
    f.vec<-c(rep(0,.sel+1),rep(f,.age.mat-.sel),rep(0,max.age-age.mat))  ## fishing only on juveniles
    si<-exp(-(.M+f.vec)) #### So this is the survivorship for each age class... Z,
    lx<-0
    lx[1]<-1                     #lx is the proportion alive at the start of each age class
    for(i in 2:(length(si)+1)) lx[i]<-lx[i-1]*si[i-1] # So here we make this the proportion alive in each age class
    
    spr<-sum(lx[1:length(si)]*.mx) # This gets the lifetime contribution of recruits for each spawner thanks to some heavy lifting from the mx vector.
    # sumsq <- (spr-(1/.litter))^2  #litter is the number of female age-0 offspring
    # DK Note: I believe we are going to want to generalize this, (see the fishmethods sbpr function) to look at other possibilities (make this a percentage of maximum spawner potential (MSP)?)
    sumsq <- (spr-1)^2 # So here the SPR is minimized against 1, which I think means the Fcrit to ensure that each mamma on average produces one recruit.
    return(sumsq)
  }

  junk<-nlminb(start = 0.1,obj = minimize, lower=0.0001,upper=4)

  # So this is F critical
  f<-junk$par
  # DK Note: This is the vector of f values for each age class, how this is set up only immature individuals are given an F.crit, will need to discuss this with Heather.
  f.vec<-c(rep(0,.sel+1),rep(f,.age.mat-.sel),rep(0,max.age-age.mat))
  # Survival vector for each age includes M and F
  si.f<-exp(-(.M+f.vec))         #values or vectors labelled .f are at f.crit
  # initialize the lx vector (survviorship)
  lx.f<-0
  lx.f[1]<-1
  # Get the survivorship and calculate spawners per recruit with the above survivorship
  for(i in 2:(length(si.f)+1)) lx.f[i]<-lx.f[i-1]*si.f[i-1]
  spr.f<-sum(lx.f[1:length(si.f)]*.mx) # This should be same value as whatever we minimized spr against in the minimize function.


  # Get the survivorship and calculate spawners per recruit based only on natural mortality.
  si.f0<-exp(-(.M))       #values or vectors labelled .f0 are in the absence of fishing
  lx.f0<-0
  lx.f0[1]<-1
  for(i in 2:(length(si.f0)+1)) lx.f0[i]<-lx.f0[i-1]*si.f0[i-1]
  spr.f0<-sum(lx.f0[1:length(si.f0)]*.mx)
  #N.crit
  #given f.crit and catch and C=N(1-exp(-F))
  #Get an estimate for N based on the average removals for the stock and the value of F critical for the population.  This N is only the N of ages that have a non-zero selectivity
  N<-removals/(1-exp(-f)) # these are within the selectivity range
  # And then you can scale this up to the total population using the survivorship.
  #Take the ratio of sum of all survivorship and divide by the sum of all survivorship in the selective ages and multiply by N to get the number of individuals at F.crit is N.crit... cool beans!
  N.crit<-N*(sum(lx.f)/sum(lx.f[(.sel+1):(.age.mat+1)]))#add in number below and above the age of selectivity, calculated using lx

  return(list(mx=mx,M=M,f.vec=f.vec,si.f=si.f,lx.f=lx.f,spr.f=spr.f,f=f,si.f0=si.f0,lx.f0=lx.f0,spr.f0=spr.f0,N=N,N.crit=N.crit))
}
