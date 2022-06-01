
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
