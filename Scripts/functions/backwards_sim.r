# The backwards simulation part of the ICM.sim function, this likely will just be me using this for testing purposes.

# All the data you need to make it dance...
# years
# n.years
# n.sims
# forward.n.years
# age.mat         The age at maturity estimated in parms_calc using life history information
# max.age         The longevity of the species estimated in parms_calc using life history information
# nat.mort:       Natural mortality, probably will be 0.2 for most stocks
# juv.mult:       A multiplier to allow for juvenile mortality if that's of interest.
# sel             The selectivity of the stock.  Currently this is set up as a single number and it is tweaked in the code, we'll need to make this more complex.
# removals        Removals from the fishery in a given year.
# N               Population size in a given year
# repo.lag:       Lag in years between giving birth and conceiving again.  Likely 0 for fish.
# gest:           The gestation period, for fish well make this 1 (i.e. this will effectively enable an individual to reproduce once a year
# u:              Exploitation rate (annual not instantaneous), currently set up to be 1 value.
# pop.model:      What method you going to use to get the population growth modeled for the backwards model.  You can use exponential model and logistic model.


icm.sim<-function(years,n.years,n.sims,forward.n.years,
                  age.mat,max.age,nat.mort,juv.mult,sel,removals,N,repo.lag,gest,u,pop.model = "exp")
{
  source("D:/Github/ICM/Scripts/functions/Lotka_r.r")
  source("D:/Github/ICM/Scripts/functions/F_crit.r")
  source("D:/Github/ICM/Scripts/functions/u_calc.r")
  source("D:/Github/ICM/Scripts/functions/backwards_project.r")
  
  #Initialize a bunch of objects
  n.years<-length(years)
  Pop<-matrix(rep(0,n.years*n.sims),n.sims,n.years)    #matrix(data=NA, nrow=<<see below>>, ncol=<<see below>>,byrow=F, dimnames=NULL) 
  Forward.pop<-matrix(rep(0,forward.n.years*n.sims),n.sims,forward.n.years)
  Pop.vec<-rep(NA,n.years)
  rem<-rep(NA,n.years)
  Forward.pop.vec<-rep(0,forward.n.years)
  r.vec<-rep(0,n.sims)
  vars.neg<-matrix(rep(0,10*n.sims*2),n.sims*2,10) ## this is how many arguemnts are in the function; 2*n.sims
  r.fishing.vec<-rep(0,n.sims)
  yr<-seq(1,forward.n.years)
  F.crit.vec<-rep(0,n.sims)
  N.crit.vec<-rep(0,n.sims)
  spr.crit.vec<-rep(0,n.sims)
  spr.f0.vec<-rep(0,n.sims)
  mean.u.vec<-rep(0,n.sims)
  B1.vec<-rep(0,n.sims) ### backwards projections 
  B2.vec<-rep(0,n.sims) ### forward projections
  #removals<-matrix(rep(0,n.years*n.sims),n.sims,n.years)

  
  #create vectors of random pars
  
  for(i in 1:(1*n.sims*2))   #2 is a patch to get n.sims values for r that are <r.cutoff 
  {
    junk<-lotka.r(age.mat.vec[i],litter.vec[i],gestation.vec[i],nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],lag.vec[i],0,0)   #exploitation rate = 0; sel=0  
    vars<-c(junk$par,age.mat.vec[i],litter.vec[i],gestation.vec[i],nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],lag.vec[i],0,0)
    #browser()
    vars.neg[i,]<-vars
    #### fishing only on juveniles OR constant selectivity on all individuals older than age 2
    junk2<-F.crit.calc(age.mat.vec[i],litter.vec[i],gestation.vec[i],nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],lag.vec[i],sel,100)
    #junk2<-F.crit.calc.original(age.mat.vec[i],litter.vec[i],gestation.vec[i],nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],lag.vec[i],sel,1) ## note, here you are taking one animal each year
    
    r.vec[i]<-junk$par 
    F.crit.vec[i]<-junk2$f      
    N.crit.vec[i]<-junk2$N.crit
    spr.crit.vec[i]<-junk2$spr.f
    spr.f0.vec[i]<-junk2$spr.f0
    
  }
  temp.r.vec<-r.vec[r.vec>0 & r.vec<r.cutoff]
  temp.F.crit.vec<-F.crit.vec[r.vec>0 & r.vec<r.cutoff]
  temp.N.crit.vec<-N.crit.vec[r.vec>0 & r.vec<r.cutoff]
  temp.spr.crit.vec<-spr.crit.vec[r.vec>0 & r.vec<r.cutoff]
  temp.spr.f0.vec<-spr.f0.vec[r.vec>0 & r.vec<r.cutoff]
  
  r.vec<-temp.r.vec[1:n.sims]
  F.crit.vec<-temp.F.crit.vec[1:n.sims]
  N.crit.vec<-temp.N.crit.vec[1:n.sims]
  spr.crit.vec<-temp.spr.crit.vec[1:n.sims]
  spr.f0.vec<-temp.spr.f0.vec[1:n.sims]
  print(r.vec)
  
  
  
  
  ############################################################
  
  #then do backward population projections for each population decline scenario
  
  for(i in 1:n.sims)
  {
    # Get the final year estimate of your population
    Pop.vec[n.years]<-N.end[i]
    #Now run your model backwards with the r from the Lotka function and
    # if you use the logistic growth model the K estimated for the population.
    for(y in n.years:2)
    {
      # DK Note: So for our removals time series, we put the removals between t+1 and t 
      # down as year t+1.  We can change this, but that's how this is set up at the moment.
      removals.next <- removals[y]
      # The exponential model
      if(pop.model == 'exp') 
      {
        exp.res <- back.proj(option = "exponential",pop.next = Pop.vec[y],r=r,removals = removals.next)
        pop.next <- exp.res$Pop.current
      }
      # If you are running the logistic model
      if(pop.model == 'logistic')
      {
        log.res <- back.proj(option = "logistic",pop.next = pop.next,K=K,r=r,removals = removals.next)
        pop.next <- min(log.res$Pop.current)
      }
      Pop.vec[y-1] <- pop.next
    } #Loop through all the years.
    
    Pop[i,]<-Pop.vec
    
    ### calc removals (increase due to r plus amount population declined by):   
    # removals[i,]<-Pop.vec-(Pop.vec/exp(r.vec[i]))+(Pop.vec*(decl.rate[i])) WRONG
    removals[i,]<-Pop.vec*(-1+exp(r.vec[i])+decl.rate[i]) ## redone with Jamie.
    
    ## removals in final year of population decline to calculate forward exploitation rate   
    rem<- removals[n.years]
    
    #print(Pop.vec)
    B1.vec[i]<-lm(log(Pop.vec)~years)$coef[2]  ## logic check
    
    
    #calculate exploitation
    junk<-c()
    for(y in 1:n.years) junk[y]<-u.calc(nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],age.mat.vec[i],sel,rem*0.7,Pop.vec[y])
    mean.u.vec[i]<-mean(junk) 
  } #end simulation
  
  
  