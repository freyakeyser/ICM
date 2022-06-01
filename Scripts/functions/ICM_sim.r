# The simulation function

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
#pop.model:       What method you going to use to get the population growth modeled for the backwards model.  You can use exponential model, logistic model, mean removals model.


icm.sim<-function(years,n.years,n.sims,forward.n.years,
                  age.mat,max.age,nat.mort,juv.mult,sel,removals,N,repo.lag,gest,u,pop.model = "Baum")
{
  source("D:/Github/ICM/Scripts/functions/decline_rate.r")
  source("D:/Github/ICM/Scripts/functions/Lotka_r.r")
  source("D:/Github/ICM/Scripts/functions/F_crit.r")
  source("D:/Github/ICM/Scripts/functions/u_calc.r")
  
  #Initialize a bunch of objects
  n.years<-length(years)
  Pop<-matrix(rep(0,n.years*n.sims),n.sims,n.years)    #matrix(data=NA, nrow=<<see below>>, ncol=<<see below>>,byrow=F, dimnames=NULL) 
  Forward.pop<-matrix(rep(0,forward.n.years*n.sims),n.sims,forward.n.years)
  Pop.vec<-rep(0,n.years)
  rem<-rep(0,n.years)
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
  removals<-matrix(rep(0,n.years*n.sims),n.sims,n.years)
  
  
  ## specific to allowable harm assessment - varying levels of fishing mortality.
  #removals<-rep(0,length(years))  ### NO harm
  #removals<-c(1,0,0,2,2,0,1,0,0,0,0,0,2,0,0,1,0,0,0,0,0,1,0,0,0,0,0,2,0,0,0,0,0,1,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0) #last 0 is a place holder, doesn't matter 
  #removals<-rep(3,67)
  #removals<-rep(10,67)
  #removals<-rep(20,67)
  
  #mean.removals<-0
  #mean.removals<-median(removals[removals>0])  #assume low levels of by-catch annually (one animal).
  #mean.removals<-3
  #mean.removals<-10
  #mean.removals<-20
  
  #sel<-2  ### this is an assumption - evaluate (age 2 become suceptible to the fishery)
  #K<-2000   #upper limit on population size used (high level of d.d.) 
  
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
    for(y in 1:(n.years-1))
    {
      pop.next <- Pop.vec[n.years-y+1]
      removals.next <- removals[n.years-y]
      r <- r.vec[i]
      ### for a population with a given decline rate (note decline rate is an increase because projecting backwards:
      
      
    }
    
    Pop[i,]<-Pop.vec
    
    ### calc removals (increase due to r plus amount population declined by):   
    # removals[i,]<-Pop.vec-(Pop.vec/exp(r.vec[i]))+(Pop.vec*(decl.rate[i])) WRONG
    removals[i,]<-Pop.vec*(-1+exp(r.vec[i])+decl.rate[i]) ## redone with Jamie.
    
    ## removals in final year of population decline to calculate forward exploitation rate   
    rem<-Pop.vec[n.years]*(-1+exp(r.vec[i])+decl.rate[i])  ## note that decline rate is a positive value
    (Pop.vec[n.years]*(decl.rate[i]))
    
    #print(Pop.vec)
    B1.vec[i]<-lm(log(Pop.vec)~years)$coef[2]  ## logic check
    
    
    #calculate exploitation
    junk<-c()
    for(y in 1:n.years)
    {
      junk[y]<-u.calc(nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],age.mat.vec[i],sel,rem*0.7,Pop.vec[y])
      #junk[y]<-u.calc.original(nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],sel,removals[y],Pop.vec[y])
      #junk[y]<-u.calc.original(nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],sel,rem*0,Pop.vec[y]) ## removals 10% of historical
    }
    mean.u.vec[i]<-mean(junk) 
  } #end simulation
  
  
  #then do forward projections:
  
  
  for(i in 1:n.sims)
  {
    #first recalc r with fishing mortality      
    junk<-lotka.r(age.mat.vec[i],litter.vec[i],gestation.vec[i],nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],lag.vec[i],mean.u.vec[i],sel)  
    r.fishing.vec[i]<-junk$par 
    
    #random deviate vectors 
    log.dev.vec<-rep(0,forward.n.years) #sigma=0 case
    r.forward<-rep(r.fishing.vec[i],forward.n.years)
    if(f.sigma>0)  #overwrite above if sigma>0 
    {
      w<-rnorm(forward.n.years,0,1) 
      log.dev.vec[1]<-w[1]*f.sigma 
      r.forward[1]<-r.fishing.vec[i]+log.dev.vec[1]#-sigma^2/2
      if(r.forward[1]>r.cutoff){r.forward[1]<-r.cutoff} #sets an upper bound on r in the sims
      
      for (y in 2:forward.n.years)
      {
        log.dev.vec[y]<-log.dev.vec[y-1]*f.auto+w[y]*f.sigma
        r.forward[y]<- r.fishing.vec[i]+log.dev.vec[y]#-sigma^2/2
        if(r.forward[y]>r.cutoff){r.forward[y]<-r.cutoff} #sets an upper bound on r in the sims}
      }    
    }
    #print(summary(r.forward))
    
    #then project forward
    Forward.pop.vec[1]<-N.end[i]
    for(y in 1:(forward.n.years-1))
    {
      Forward.pop.vec[y+1]<-Forward.pop.vec[y]*exp(r.forward[y])
    }
    Forward.pop[i,]<-Forward.pop.vec
    #print(Forward.pop[i,])
    #browser()
    B2.vec[i]<-lm(log(Forward.pop.vec)~yr)$coef[2]
    
    
  } #end simulation
  
  
  end.pop.size<-Forward.pop[,forward.n.years]
  
  
  r.summary<-quantile(r.vec,c(0.1,0.5,0.9))   
  F.crit.summary<-quantile(F.crit.vec,c(0.1,0.5,0.9))   
  N.crit.summary<-quantile(N.crit.vec,c(0.1,0.5,0.9)) 
  spr.crit.summary<-quantile(spr.crit.vec, c(0.1,0.5,0.9))
  spr.f0.summary<-quantile(spr.f0.vec,c(0.1,0.5,0.9))
  proportion.declining<-length(B2.vec[B2.vec<0])/length(B2.vec)
  end.pop.summary<-quantile(end.pop.size,c(0.1,0.5,0.9)) 
  prop.ext<-length(end.pop.size[end.pop.size<q.ext])/length(end.pop.size)
  
  junk<-Forward.pop[,forward.n.years]/Forward.pop[,1]
  prop.curtis<-(length(junk[junk>2]))/n.sims ## proportion of sims that doubled
  prop.10<-(length(junk[junk>1.1]))/n.sims ## proportion that increased by 10%
  
  
  #### NOTE - this is NOT what was done for the manuscript
  #PBR<-r.vec*0.5*Forward.pop[,1]*0.5 ### potential biological removals where f=0.5
  #PBR1<-r.vec*0.5*Forward.pop[,1]*0.1  ## f=0.1
  
  return(list(vars.neg=vars.neg,Pop=Pop,removals=removals,Forward.pop=Forward.pop,r.vec=r.vec,r.fishing.vec=r.fishing.vec,r.forward=r.forward,B1.vec=B1.vec,B2.vec=B2.vec,mean.u.vec=mean.u.vec,r.summary=r.summary,
              F.crit.vec=F.crit.vec,F.crit.summary=F.crit.summary,N.crit.summary=N.crit.summary,proportion.declining=proportion.declining,
              end.pop.summary=end.pop.summary,prop.ext=prop.ext,decl.rate=decl.rate,prop.curtis=prop.curtis,
              prop.10=prop.10,spr.crit.summary=spr.crit.summary,spr.f0.summary=spr.f0.summary,spr.f0.vec=spr.f0.vec))
  
  ### note Jamie's additions - for plotting histograms - also increase # years to 60 to get doubling time
  # kick out Fcrit values not summary
  #kick out spr.f0 values not summary
  #calculate median # years to doubling. see end of next section
  
}