### white shark - evaluation of life history information relative to status.
## started Feb 2018 - HB

################################################

## life history functions used in simulation

###################################################

setwd('C:/mydocs/white')

lotka.r<-function(age.mat,litter,gestation,nat.mort,juv.mult,max.age,lag.f,u,sel) 
{
  
  # age.mat = 9.5
  # litter = 8 (both sexes) - 4 (female only) ### NOTE Euler-Lotka is relative to female pups only.
  # gest. = 2.5 yrs
  # M = 0.112
  # juv.mult = 0
  # max.age = 40
  # lag.f = 0 (lag in years between giving birth and concieving again)
  # u = 0 (no fishery)
  # sel = 0 (no fishery)
  
  gest.lag<-gestation+lag.f
  
  temp<-c(0,c(rep(0,age.mat+gestation-1),litter,rep(c(rep(0,gest.lag-1),litter),100))) 
  mx<-temp[1:(max.age+1)]
  
  r.M<-c(nat.mort*juv.mult,rep(nat.mort,max.age-1))
  
  
  f<--log(1-u)
  f.vec<-c(rep(0,sel),rep(f,max.age-sel))      
  si<-exp(-(r.M+f.vec)) 
  
  x<-0:(max.age)
  lx<-0
  lx[1]<-1
  for(i in 2:(length(si)+1))
  {lx[i]<-lx[i-1]*si[i-1]}
  
  
  lxmx<-lx*mx
  assign(".lxmx",lxmx)     # Store in expression frame
  assign(".x",x)
  #print("expected female output from one newly born female")
  #print(sum(lxmx))
  #print("mx")
  #print(mx)

  ### here is actual function to minimize; minimization done based on the sum of squares.
  minimize<-function(start.r)
  {
    rx<-exp(-start.r*.x)
    lotka<-sum(.lxmx*rx)
    sumsq <- sum((lotka-1)^2) 
    return(sumsq)
  }
  junk<-nlminb(start = 0.03,objective = minimize, lower=-1,upper=1)
  return(junk)
}

F.crit.calc.original<-function(age.mat,litter,gestation,nat.mort,juv.mult,max.age,lag.f,sel,removals)
{
  #this function finds the value of F such that the population growth rate = 1
  #this is based on an SPR calculation
  #the population growth rate =1 where SPR = 1/number of age-0 female offspring
  #fishing selectivity is knife-edge at age=sel
  #N.crit is the size the total population would have to be such that F would be less than F.crit 
  #here, removals are assumed to be 1 animal
  
  gest.lag<-gestation+lag.f
  temp<-c(0,c(rep(0,age.mat+gestation-1),litter,rep(c(rep(0,gest.lag-1),litter),100))) 
  mx<-temp[1:(max.age+1)]
  
  M<-c(nat.mort*juv.mult,rep(nat.mort,max.age-1))
  
  assign(".M",M)     # Store in expression frame
  assign(".mx",mx)
  assign(".sel",sel)
  assign(".max.age",max.age)
  assign(".litter",litter)
  minimize<-function(f)
  {
    f.vec<-c(rep(0,.sel),rep(f,.max.age-.sel))         
    si<-exp(-(.M+f.vec)) ### survivorship
    lx<-0
    lx[1]<-1                     #lx is the proportion alive at the start of each age class
    for(i in 2:(length(si)+1))
    {lx[i]<-lx[i-1]*si[i-1]}
    spr<-sum(lx*.mx)
    # sumsq <- (spr-(1/.litter))^2  #litter is the number of female age-0 offspring
    sumsq <- (spr-1)^2
    return(sumsq)
  }
  
  junk<-nlminb(start = 0.1,obj = minimize, lower=0.0001,upper=4)
  
  
  f<-junk$par
  f.vec<-c(rep(0,.sel),rep(f,.max.age-.sel))          
  si.f<-exp(-(.M+f.vec))         #values or vectors labelled .f are at f.crit
  lx.f<-0
  lx.f[1]<-1
  for(i in 2:(length(si.f)+1))
  {lx.f[i]<-lx.f[i-1]*si.f[i-1]}
  spr.f<-sum(lx.f*.mx)
  
  si.f0<-exp(-(.M))       #values or vectors labelled .f0 are in the absence of fishing 
  lx.f0<-0
  lx.f0[1]<-1
  for(i in 2:(length(si.f0)+1))
  {lx.f0[i]<-lx.f0[i-1]*si.f0[i-1]}
  spr.f0<-sum(lx.f0*.mx)
  
  
  
  #N.crit
  #given f.crit and catch and C=N(1-exp(-F))
  N<-removals/(1-exp(-f)) # these are at or above the age of selectivity
  #add in number below the age of selectivity, calculated using lx
  N.crit<-N*(sum(lx.f)/sum(lx.f[(sel+1):length(lx.f)]))
  
  
  return(list(mx=mx,M=M,f.vec=f.vec,si.f=si.f,lx.f=lx.f,spr.f=spr.f,f=f,si.f0=si.f0,lx.f0=lx.f0,spr.f0=spr.f0,N=N,N.crit=N.crit))
}

F.crit.calc<-function(age.mat,litter,gestation,nat.mort,juv.mult,max.age,lag.f,sel,removals)
{
  #this function finds the value of F such that the population growth rate = 1
  #this is based on an SPR calculation
  #the population growth rate =1 where SPR = 1/number of age-0 female offspring
  #fishing selectivity is knife-edge at age=sel
  #N.crit is the size the total population would have to be such that F would be less than F.crit 
  
  gest.lag<-gestation+lag.f
  temp<-c(0,c(rep(0,age.mat+gestation-1),litter,rep(c(rep(0,gest.lag-1),litter),100))) 
  mx<-temp[1:(max.age+1)]
  
  M<-c(nat.mort*juv.mult,rep(nat.mort,max.age-1))
  
  assign(".M",M)     # Store in expression frame
  assign(".mx",mx)
  assign(".sel",sel)
  assign(".max.age",max.age)
  assign(".age.mat",age.mat)
  assign(".litter",litter)
  minimize<-function(f)
  {
    f.vec<-c(rep(0,.sel),rep(f,.age.mat-.sel),rep(0,max.age-age.mat))      
    si<-exp(-(.M+f.vec)) ### survivorship
    lx<-0
    lx[1]<-1                     #lx is the proportion alive at the start of each age class
    for(i in 2:(length(si)+1))
    {lx[i]<-lx[i-1]*si[i-1]}
    spr<-sum(lx*.mx)
    # sumsq <- (spr-(1/.litter))^2  #litter is the number of female age-0 offspring
    sumsq <- (spr-1)^2
    return(sumsq)
  }
  
  junk<-nlminb(start = 0.1,obj = minimize, lower=0.0001,upper=4)
  
  
  f<-junk$par
  f.vec<-c(rep(0,.sel),rep(f,.age.mat-.sel),rep(0,max.age-age.mat))       
  si.f<-exp(-(.M+f.vec))         #values or vectors labelled .f are at f.crit
  lx.f<-0
  lx.f[1]<-1
  for(i in 2:(length(si.f)+1))
  {lx.f[i]<-lx.f[i-1]*si.f[i-1]}
  spr.f<-sum(lx.f*.mx)
  
  si.f0<-exp(-(.M))       #values or vectors labelled .f0 are in the absence of fishing 
  lx.f0<-0
  lx.f0[1]<-1
  for(i in 2:(length(si.f0)+1))
  {lx.f0[i]<-lx.f0[i-1]*si.f0[i-1]}
  spr.f0<-sum(lx.f0*.mx)
  
  #browser()
  
  #N.crit
  #given f.crit and catch and C=N(1-exp(-F))
  N<-removals/(1-exp(-f)) # these are within the selectivity range
  #add in number below and above the age of selectivity, calculated using lx
  N.crit<-N*(sum(lx.f)/sum(lx.f[(.sel+1):(.age.mat+1)]))
  
  
  return(list(mx=mx,M=M,f.vec=f.vec,si.f=si.f,lx.f=lx.f,spr.f=spr.f,f=f,si.f0=si.f0,lx.f0=lx.f0,spr.f0=spr.f0,N=N,N.crit=N.crit))
}

## note that this includes fishing on adults - would have to restrict if decide to use
# F.crit.calc function above.
u.calc<-function(nat.mort,juv.mult,max.age,age.mat,sel,removals,N)
{
  #this function finds the value of f (of the vulnerable, not total, population taken each year)
  #this is based on an equilibrium age structure assumption
  #selectivity is knife-edged
  
  
  M<-c(nat.mort*juv.mult,rep(nat.mort,max.age-1))
  
  assign(".M",M)     # Store in expression frame
  assign(".sel",sel)
  assign(".max.age",max.age)
  assign(".age.mat",age.mat)
  assign(".N",N)
  assign(".removals",removals)
  
  minimize<-function(f)
  {
    f.vec<-c(rep(0,.sel),rep(f,.age.mat-.sel),rep(0,max.age-age.mat))        
    si<-exp(-(.M+f.vec)) 
    lx.f<-0
    u<-1-exp(-f)
    lx.f[1]<-1                     #lx is the proportion alive at the start of each age class
    for(i in 2:(length(si)+1))
    {lx.f[i]<-lx.f[i-1]*si[i-1]}
    
    #add in number below the age of selectivity, calculated using lx
    vulnerable<-.N*(sum(lx.f[(.sel+1):length(lx.f)])/sum(lx.f))
    #   print(vulnerable)
    pred.C<-vulnerable*u
    #   print(pred.C) 
    sumsq <- (.removals-pred.C)^2
    return(sumsq)
  }
  
  junk<-nlminb(start = 0.1,obj = minimize, lower=0.0001,upper=10)
  
  f.out<-junk$par
  u.out<-1-exp(-f.out)
  
  return(u.out)
}


##############################################################################

#### white parms 

###############################################################################


parms.short<-function()
{
  age.mat<<-9
  litter<<-4
  gestation<<-2
  nat.mort<<-0.112
  juv.mult<<-1
  max.age<<-40
  lag.f<<-0
  u<<-0
  sel<<-2
  
  # deterministic Rmax
  det.r<-lotka.r(age.mat,litter,gestation,nat.mort,juv.mult,max.age,lag.f,u,sel)
  print("determinisitc rmax")
  print(det.r$par)
  
  r.cutoff<<-2*det.r$par 
  
  ## MC parameters
  n.sims<<-1000 ## SET NUMBER OF SIMS
  
  age.mat.vec<<-sample(8:11,n.sims*2,replace=T) #9.5 - note that proportionately, this is about the same level of variability (16% above, 16% below estimate - same as in short scenario)
  litter.vec<<-sample(2:6,n.sims*2,replace=T)    #8 total switched to females =4
  gestation.vec<<-sample(2:3,n.sims*2,replace=T) #2.5

  nat.mort<-seq(0.062,0.162,by=0.001)
  nat.mort.vec<<-sample(nat.mort,n.sims*2,replace = T) #0.112+0.05
  juv.mult.vec<<-rep(1,n.sims*2) #0
  lag.vec<<-rep(0,n.sims*2) #0
  
  max.age.vec<<-sample(35:45,n.sims*2,replace=T) #30
  N.end<<-c(sample(2000:2499,size=n.sims/2,replace=T),sample(2500:3000,size=n.sims/2,replace=T)) # see Burgess et al. 2014

  q.ext<<-20 #quasi-extinction threshold for forward projections
  f.sigma<<-0.01 #lognormal shape parameter for forward projections
  f.auto<<-0.03 #autocorrelation parameter for forward projections
  forward.n.years<<-30 #project for 30 years (i.e. possible increase in Curtis)
  
}

parms.long<-function()
{
  age.mat<<-30
  litter<<-4
  gestation<<-2
  nat.mort<<-0.063
  juv.mult<<-1
  max.age<<-70
  lag.f<<-0
  u<<-0
  sel<<-0
  
  # deterministic Rmax
  det.r<-lotka.r(age.mat,litter,gestation,nat.mort,juv.mult,max.age,lag.f,u,sel)
  #print("determinisitc rmax")
  #print(det.r$par)
  
  r.cutoff<<-2*det.r$par 
  
    # MC parameters
  n.sims<<-1000  ## SET NUMBER OF SIMS
  
  age.mat.vec<<-sample(25:35,n.sims*4,replace=T) #30 - note that proportionately, this is about the same level of variability (16% above, 16% below estimate - same as in short scenario)
  litter.vec<<-sample(2:6,n.sims*4,replace=T)    #8 total switched to females =4
  gestation.vec<<-sample(2:3,n.sims*4,replace=T) #2.5
  
  nat.mort<-seq(0.053,0.073,by=0.0001)
  nat.mort.vec<<-sample(nat.mort,n.sims*4,replace = T) #0.063+0.01
  juv.mult.vec<<-rep(1,n.sims*4) #0
  lag.vec<<-rep(0,n.sims*4) #0
  
  max.age.vec<<-sample(60:80,n.sims*4,replace=T) #70
  N.end<<-c(sample(2000:2499,size=n.sims/2,replace=T),sample(2500:3000,size=n.sims/2,replace=T)) # see Burgess et al. 2014
  
  q.ext<<-20 #quasi-extinction threshold for forward projections
  f.sigma<<-0.01 #lognormal shape parameter for forward projections
  f.auto<<-0.03 #autocorrelation parameter for forward projections
  forward.n.years<<-30 #first year gets n.end and then project for 30 years as in Curtis.
  
}

Baum<-function() 
{years<<-1986:2000 #15 years
decl.junk<-seq(1.058,1.137,by=0.0001)  #  200 values - note that this is the old parameterization of decline
decl.rate<<-sample(decl.junk,n.sims*4,replace=T) 
## decline rate 59% over 15 years to 89% over 15 years -Baum et al. 2003)
#to calculate annual rate from total decline:
#1-.59 = 0.41 to 1-.89 = 0.11  total decline (beginning value = 1; ending value = 0.41) ending/beginning
#0.41^(1/14)-1 = -0.05770781 per year ### NOTE - 14 time steps over 15 years
#0.11^(1/14)-1 = -0.1368369 per year  
}

Curtis<-function()
{years<<-1970:1989 ## 20 years
decl.junk<-seq(0.05098,0.06659,by=0.0001)  #  approx 150 values
decl.rate<<-sample(decl.junk,n.sims*4,replace=T) 
## 63% to 73% decline reported in Curtis et al. (2014)
#to calculate annual rate from total decline:
#1-.63 = 0.37 to 1-.73 = 0.27  total decline
#0.37^(1/19)-1 =  -0.05098347 per year ## 19 timesteps over 20 years.
#0.27^(1/19)-1 =  -0.06659144 per year 
}


#################################################################################################

#### simulations

#####################################################################################################
# to run analysis
#parms.short() or parms.long()
#Baum() or Curtis()
#white.sim.result<-white.sim

## NOTE for MS, only used Curtis decline rate.
set.seed(10)

white.sim<-function()
{

  
  #ajf gibson Apr. 26/08
  ## modified Heather Sept.2016
  #assign output to "white.sim.result" which is an object used by the plotting function below 
  
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
    Pop.vec[n.years]<-N.end[i]
    for(y in 1:(n.years-1))
    {
      #Pop.vec[n.years-y]<-(Pop.vec[n.years-y+1]+removals[n.years-y])/exp(r.vec[i])  ### simple exponential population growth
      #Pop.vec[n.years-y]<-(Pop.vec[n.years-y+1]+removals[n.years-y])/(exp(r.vec[i]))*(1-(Pop.vec[n.years-y+1]/K)) #this is approximate at the moment as the density dependence should work forward
      #Pop.vec[n.years-y]<-(Pop.vec[n.years-y+1]+mean.removals)/exp(r.vec[i]) #average value

      ### for a population with a given decline rate (note decline rate is an increase because projecting backwards:
      Pop.vec[n.years-y]<-Pop.vec[n.years-y+1]/(1-decl.rate[i]) 
      
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

### create object: white.sim.result











######################

## Manuscript info

#######################

y<-data.frame(r.long=white.sim.LC0$r.summary,r.short=white.sim.SC0$r.summary,
    spr.long=white.sim.LC0$spr.f0.summary,spr.short=white.sim.SC0$spr.f0.summary)

## pre-decline abundance
z<-data.frame(short=round(quantile(white.sim.SC0$Pop[,1],c(0.1,0.5,0.9)),0),long=round(quantile(white.sim.LC0$Pop[,1],c(0.1,0.5,0.9)),0))

## response to decreases in mortality
## note - only using the fast decline scenario here
### also note - forward.n.years has to = 30 for this code to work!!!
response<-data.frame(reduction=c('30%','50%','70%','90%','100%'),
                     prop10.S=c(white.sim.SC3$prop.10,white.sim.SC5$prop.10,white.sim.SC7$prop.10,white.sim.SC9$prop.10,white.sim.SC0$prop.10),
                     propC.S=c(white.sim.SC3$prop.curtis,white.sim.SC5$prop.curtis,white.sim.SC7$prop.curtis,white.sim.SC9$prop.curtis,white.sim.SC0$prop.curtis),
                     
                     prop10.L=c(white.sim.LC3$prop.10,white.sim.LC5$prop.10,white.sim.LC7$prop.10,white.sim.LC9$prop.10,white.sim.LC0$prop.10),
                     propC.L=c(white.sim.LC3$prop.curtis,white.sim.LC5$prop.curtis,white.sim.LC7$prop.curtis,white.sim.LC9$prop.curtis,white.sim.LC0$prop.curtis))

## Status relative to reference points

## Nmin calculations 
Nmin.SC0<-c()
for (i in 1:20)
{Nmin.SC0[i]<-quantile(white.sim.SC0$Pop[,i],c(0.2))  }
Nmin.SC0

#Nmin.SB0<-c()
#for (i in 1:15)
#{Nmin.SB0[i]<-quantile(white.sim.SB0$Pop[,i],c(0.2))  }
#Nmin.SB0

Nmin.LC0<-c()
for (i in 1:20)
{Nmin.LC0[i]<-quantile(white.sim.LC0$Pop[,i],c(0.2))  }
Nmin.LC0

#Nmin.LB0<-c()
#for (i in 1:15)
#{Nmin.LB0[i]<-quantile(white.sim.LB0$Pop[,i],c(0.2))  }
#Nmin.LB0

### PBR, f=0.1, 20th percentile population sizes as above.
ref.point<-data.frame(PBR1.SC0=sum(median(white.sim.SC0$r.vec)*0.5*Nmin.SC0*0.1),
                      #PBR1.SB0=sum(median(white.sim.SB0$r.vec)*0.5*Nmin.SB0*0.1),
                      #PBR1.LB0=sum(median(white.sim.LB0$r.vec)*0.5*Nmin.LB0*0.1),
                      PBR1.LC0=sum(median(white.sim.LC0$r.vec)*0.5*Nmin.LC0*0.1))  ## f=0.1
## for plotting
PBR1.SC0=data.frame(years=c(1970:1989),pbr=(median(white.sim.SC0$r.vec)*0.5*Nmin.SC0*0.1))
#PBR1.SB0=data.frame(years=c(1986:2000),pbr=(median(white.sim.SB0$r.vec)*0.5*Nmin.SB0*0.1))
#PBR1.LB0=data.frame(years=c(1986:2000),pbr=(median(white.sim.LB0$r.vec)*0.5*Nmin.LB0*0.1))
PBR1.LC0=data.frame(years=c(1970:1989),pbr=(median(white.sim.LC0$r.vec)*0.5*Nmin.LC0*0.1))

# median removals by year
rem.SC0<-c()
SC0.90<-c()
SC0.10<-c()
for (i in 1:20)
{rem.SC0[i]<-median(white.sim.SC0$removals[,i]) 
SC0.90[i]<-quantile(white.sim.SC0$removals[,i],c(0.9))
SC0.10[i]<-quantile(white.sim.SC0$removals[,i],c(0.1))}
SC0<-data.frame(removals=rem.SC0,q.90=SC0.90,q.10=SC0.10)
SC0$sim<-'Short lifespan; slow decline'
SC0$year<-c(1970:1989)

rem.SB0<-c()
SB0.90<-c()
SB0.10<-c()
for (i in 1:15)
{rem.SB0[i]<-median(white.sim.SB0$removals[,i]) 
SB0.90[i]<-quantile(white.sim.SB0$removals[,i],c(0.9))
SB0.10[i]<-quantile(white.sim.SB0$removals[,i],c(0.1))}
SB0<-data.frame(removals=rem.SB0,q.90=SB0.90,q.10=SB0.10)
SB0$sim<-'Short lifespan; fast decline'
SB0$year<-c(1986:2000)

rem.LC0<-c()
LC0.90<-c()
LC0.10<-c()
for (i in 1:20)
{rem.LC0[i]<-median(white.sim.LC0$removals[,i]) 
LC0.90[i]<-quantile(white.sim.LC0$removals[,i],c(0.9))
LC0.10[i]<-quantile(white.sim.LC0$removals[,i],c(0.1))}
LC0<-data.frame(removals=rem.LC0,q.90=LC0.90,q.10=LC0.10)
LC0$sim<-'Long lifespan; slow decline'
LC0$year<-c(1970:1989)

rem.LB0<-c()
LB0.90<-c()
LB0.10<-c()
for (i in 1:15)
{rem.LB0[i]<-median(white.sim.LB0$removals[,i]) 
LB0.90[i]<-quantile(white.sim.LB0$removals[,i],c(0.9))
LB0.10[i]<-quantile(white.sim.LB0$removals[,i],c(0.1))}
LB0<-data.frame(removals=rem.LB0,q.90=LB0.90,q.10=LB0.10)
LB0$sim<-'Long lifespan; fast decline'
LB0$year<-c(1986:2000)

### removals relative to PBR
library(ggplot2)
library(reshape2)


### make 4 ggplots with inset panels
# combine with grid.arrange

p1<-ggplot(SB0, aes(year,y=removals))+geom_line()
g<-ggplotGrob(ggplot(PBR1.SB0,aes(years,pbr))+geom_line()+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())+scale_y_continuous(limits=c(0,25)))
p1<-p1+annotation_custom(grob=g,xmin=1970,xmax=1985,ymin=1500,ymax=3000)
p1<-p1+geom_ribbon(aes(ymin=SB0$q.10,ymax=SB0.90),alpha=0.2)
p1<-p1+scale_x_continuous(limits=c(1970,2000))+scale_y_continuous(limits=c(0,3000))
p1<-p1+xlab('')+ylab('Removals (# animals)')+ggtitle('Short lifespan; fast decline')
p1<-p1+theme(axis.title.y = element_text(size=rel(1.2)))

p2<-ggplot(SC0, aes(year,y=removals))+geom_line()
g<-ggplotGrob(ggplot(PBR1.SC0,aes(years,pbr))+geom_line()+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())+scale_y_continuous(limits=c(0,26)))
p2<-p2+annotation_custom(grob=g,xmin=1981,xmax=1989,ymin=1200,ymax=2000)
p2<-p2+geom_ribbon(aes(ymin=SC0$q.10,ymax=SC0.90),alpha=0.2)
p2<-p2+scale_x_continuous(limits=c(1970,1990))+scale_y_continuous(limits=c(0,2000))
p2<-p2+xlab('')+ylab('Removals (# animals)') ##+ggtitle('Short lifespan; slow decline')
p2<-p2+theme(axis.title.y = element_text(size=rel(1.5)))

p3<-ggplot(LB0, aes(year,y=removals))+geom_line()
g<-ggplotGrob(ggplot(PBR1.LB0,aes(years,pbr))+geom_line()+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())+scale_y_continuous(limits=c(0,25)))
p3<-p3+annotation_custom(grob=g,xmin=1970,xmax=1985,ymin=1500,ymax=3000)
p3<-p3+geom_ribbon(aes(ymin=LB0$q.10,ymax=LB0.90),alpha=0.2)
p3<-p3+scale_x_continuous(limits=c(1970,2000))+scale_y_continuous(limits=c(0,2000))
p3<-p3+xlab('Year')+ylab('Removals (# animals)')+ggtitle('Long lifespan; fast decline')
p3<-p3+theme(axis.title.x =  element_text(size=rel(1.2)))
p3<-p3+theme(axis.title.y = element_text(size=rel(1.2)))

p4<-ggplot(LC0, aes(year,y=removals))+geom_line()
g<-ggplotGrob(ggplot(PBR1.LC0,aes(years,pbr))+geom_line()+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())+scale_y_continuous(limits=c(0,26)))
p4<-p4+annotation_custom(grob=g,xmin=1981,xmax=1989,ymin=1000,ymax=1700)
p4<-p4+geom_ribbon(aes(ymin=LC0$q.10,ymax=LC0.90),alpha=0.2)
p4<-p4+scale_x_continuous(limits=c(1970,1990))+scale_y_continuous(limits=c(0,1700))
p4<-p4+xlab('Year')+ylab('Removals (# animals)')##+ggtitle('Long lifespan; slow decline')
p4<-p4+theme(axis.title.x =  element_text(size=rel(1.5)))
p4<-p4+theme(axis.title.y = element_text(size=rel(1.5)))

## build the figure
library(gridExtra)
Fig2<-grid.arrange(p2,p4)
Fig2

## high res
png(file="Fig2.png",units='in', width=6, height=8,res=800)
grid.arrange(p2,p4)
dev.off()

#############
junk1<-white.sim.LC0$Forward.pop  ## project forward 100 years
junk2<-white.sim.SC0$Forward.pop

junk11<-ifelse(junk1[,]>2*junk1[,1],0,1)
junk22<-ifelse(junk2[,]>2*junk2[,1],0,1)

white.sim.LC0$doubling<-rowSums(junk11)
white.sim.SC0$doubling<-rowSums(junk22)


# build data
test<-data.frame(lifespan=c(rep('long',4000),rep('short',4000)),
                scenario= rep(c(rep('Growth rate (r)',1000),rep('Instantaneous fishing mortality (F.crit)',1000), rep('Population doubling time',1000),rep('Lifetime reproductive output',1000)),2),
           value=c(white.sim.LC0$r.vec, white.sim.LC0$F.crit.vec,white.sim.LC0$doubling,white.sim.LC0$spr.f0.vec,
                   white.sim.SC0$r.vec, white.sim.SC0$F.crit.vec, white.sim.SC0$doubling,white.sim.SC0$spr.f0.vec))
## get rid of the weird doubling-time estimates
test<-test[!test$value%in%c(4,5,6,100),]

### create medians
library(dplyr) # With dplyr for example
test <- test %>% group_by(scenario,lifespan) %>%
  mutate(med = median(value))

### histograms of parameters.
p1<-ggplot(test,aes(x=value, fill=lifespan,color=lifespan))+
  geom_histogram(position="identity",alpha=0.85)+geom_vline(aes(xintercept=med))+ 
  facet_wrap(~scenario,scales='free')+ scale_fill_grey()+labs(x = 'Value', y = 'Count')

## high res
png(file="Fig1.png",units='in', width=6, height=6,res=800)
p1
dev.off()

       
#############################################

period<-c(1970:2019)
zero.mort<-c(apply(white.sim.LC0$Pop,2,FUN=median),apply(white.sim.LC9$Forward.pop[,1:30],2,FUN=median))
zero.mort1<-c(apply(white.sim.LC0$Pop,2,FUN=quantile, prob=0.1),apply(white.sim.LC9$Forward.pop[,1:30],2,FUN=quantile, prob=0.1))
zero.mort9<-c(apply(white.sim.LC0$Pop,2,FUN=quantile, prob=0.9),apply(white.sim.LC9$Forward.pop[,1:30],2,FUN=quantile, prob=0.9))

mort3<-c(apply(white.sim.LC0$Pop,2,FUN=median),apply(white.sim.LC3$Forward.pop[,1:30],2,FUN=median))
mort31<-c(apply(white.sim.LC0$Pop,2,FUN=quantile, prob=0.1),apply(white.sim.LC3$Forward.pop[,1:30],2,FUN=quantile, prob=0.1))
mort39<-c(apply(white.sim.LC0$Pop,2,FUN=quantile, prob=0.9),apply(white.sim.LC3$Forward.pop[,1:30],2,FUN=quantile, prob=0.9))

mort5<-c(apply(white.sim.LC0$Pop,2,FUN=median),apply(white.sim.LC5$Forward.pop[,1:30],2,FUN=median))
  
mort7<-c(apply(white.sim.LC0$Pop,2,FUN=median),apply(white.sim.LC7$Forward.pop[,1:30],2,FUN=median))

D<-data.frame(period=period,zero.mort=zero.mort,zero.mort1=zero.mort1,zero.mort9=zero.mort9,
              mort3=mort3,mort31=mort31,mort39=mort39)

#plot(period, zero.mort, type='l',ylim=c(0,15000))

fig3<-ggplot(D, aes(period,y=zero.mort))+geom_line(cex=1)
fig3<-fig3+geom_ribbon(aes(ymin=zero.mort1,ymax=zero.mort9),alpha=0.2)
fig3<-fig3+xlab('Year')+ylab('Abundance') 
fig3<-fig3+theme(axis.title.y = element_text(size=rel(1.2)))+theme(axis.title.x = element_text(size=rel(1.2)))
fig3<-fig3+geom_line(aes(period, mort3),lty=5,cex=1)+geom_ribbon(aes(ymin=mort31,ymax=mort39),alpha=0.4)
fig3

Db<-data.frame(period=c(rep(period,4)),median=c(zero.mort,mort3,mort5,mort7),scenario=c(rep('90%',length=length(period)),rep('30%',length=length(period)),rep('50%',length=length(period)),rep('70%',length=length(period))))
  
fig3b<-ggplot(Db, aes(period,y=median,lty=scenario))+geom_line(cex=1)
fig3b<-fig3b+xlab('Year')+ylab('Abundance') 
fig3b<-fig3b+theme(axis.title.y = element_text(size=rel(1.2)))+theme(axis.title.x = element_text(size=rel(1.2)))
fig3b


## high res
png(file="Fig3b.png",units='in', width=6, height=6,res=800)
fig3b
dev.off()










#####################################################

## old summary plots

#######################################################
#### if need to store multiple panels in a list and then combine
p.list = lapply(sort(unique(all.rem$sim)), function(i) {
  ggplot(all.rem[all.rem$sim==i,], aes(Var2,median(value))) +geom_smooth(se=FALSE)
})
p2.list<-lapply(sort(unique(ref.annual$sim)), function(i) {
  ggplotGrob(ggplot(ref.annual[ref.annual$sim==i,], aes(years,PBR)) +geom_line())
})


plot.white.sim<-function()
{
  pop<-white.sim.result$Pop
  remov<-white.sim.result$removals[,2:length(years)]
  
  r<-white.sim.result$r.vec
  B1<-white.sim.result$B1.vec
  ann.change<-white.sim.result$ann.change
  
  #### all simulations
  windows()
  #par(mfcol=c(2,1),mar=c(5,5,2,2))
  #plot(years,pop[1,],type="l",ylim=c(0,max(pop[,1])),xlab="Year",ylab="Abundance",lty=1)
  #for (i in 2:n.sims)
  #{lines (years,pop[i,])}
  
  ## boxplots
  #library(ggplot2)
  #library(reshape2)
  #remov.m <- melt(remov)
  #ggplot(remov.m, aes(as.factor(Var1), value)) + geom_boxplot()+
  #  stat_summary(fun.y=mean, geom="line", aes(group=1),col='red',lwd=2)
  ## can't really see the mean line anyways - too many simulations.
  
  boxplot(t(remov),medcol="red",whisklty = 0, staplelty = 0,xlab="Simulation number",ylab='Removals (# of animals)')
  
  windows()
  par(mfcol=c(2,2),omi=c(0.5,0.7,0.5,0.25),mar=c(4,4,2,2)) 
  
  hist(pop[,length(years)],nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  mtext("N (2000)",1,outer=F,cex=1.2,line=3)
  abline(v=median(pop[,length(years)]),col='red')
  
  hist(pop[,1],nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  mtext("N (1986)",1,outer=F,cex=1.2,line=3)
  abline(v=median(pop[,1]),col='red')
  
  
  hist(r,nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  #mtext("Probability Density",2,outer=T,cex=1.4,line=2)
  mtext("r",1,outer=F,cex=1.2,line=3)
  abline(v=median(r),col='red')
  
  ## instead of the annual rate of change, it would be useful to plot median removals.
  hist(rowMeans(remov),nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  abline(v=median(rowMeans(remov)),col='red')
  mtext("Mean Number of Removals",1,outer=F,cex=1.2,line=3)
  
  #hist(B1,nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  #mtext("Annual Rate of Change",1,outer=F,cex=1.2,line=3)
  #mtext("Slope",1,outer=F,cex=1.2,line=3)
  #abline(v=median(B1),col='red')
  
}

