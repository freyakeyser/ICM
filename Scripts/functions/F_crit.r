# Finds the critical value of F I think...



F.crit<-function(age.mat,max.age,mx,r.M,sel,removals)
{
  ## note that mx and r.M are outputs of the parameter function
  M<-r.M
  
  assign(".M",M)     # Store in expression frame
  assign(".mx",mx)
  assign(".sel",sel)
  assign(".max.age",max.age)
  assign(".age.mat",age.mat)
  minimize<-function(f)
  {
    f.vec<-c(rep(0,.sel+1),rep(f,.age.mat-.sel),rep(0,max.age-age.mat))  ## fishing only on juveniles     
    si<-exp(-(.M+f.vec)) ### survivorship
    lx<-0
    lx[1]<-1                     #lx is the proportion alive at the start of each age class
    for(i in 2:(length(si)+1))
    {lx[i]<-lx[i-1]*si[i-1]}
    spr<-sum(lx[1:length(si)]*.mx)
    # sumsq <- (spr-(1/.litter))^2  #litter is the number of female age-0 offspring
    sumsq <- (spr-1)^2
    return(sumsq)
  }
  
  junk<-nlminb(start = 0.1,obj = minimize, lower=0.0001,upper=4)
  
  
  f<-junk$par
  f.vec<-c(rep(0,.sel+1),rep(f,.age.mat-.sel),rep(0,max.age-age.mat))       
  si.f<-exp(-(.M+f.vec))         #values or vectors labelled .f are at f.crit
  lx.f<-0
  lx.f[1]<-1
  for(i in 2:(length(si.f)+1))
  {lx.f[i]<-lx.f[i-1]*si.f[i-1]}
  spr.f<-sum(lx.f[1:length(si.f)]*.mx)
  
  si.f0<-exp(-(.M))       #values or vectors labelled .f0 are in the absence of fishing 
  lx.f0<-0
  lx.f0[1]<-1
  for(i in 2:(length(si.f0)+1))
  {lx.f0[i]<-lx.f0[i-1]*si.f0[i-1]}
  spr.f0<-sum(lx.f0[1:length(si.f0)]*.mx)
  
  #N.crit
  #given f.crit and catch and C=N(1-exp(-F))
  N<-removals/(1-exp(-f)) # these are within the selectivity range
  #add in number below and above the age of selectivity, calculated using lx
  N.crit<-N*(sum(lx.f)/sum(lx.f[(.sel+1):(.age.mat+1)]))
  
  return(list(mx=mx,M=M,f.vec=f.vec,si.f=si.f,lx.f=lx.f,spr.f=spr.f,f=f,si.f0=si.f0,lx.f0=lx.f0,spr.f0=spr.f0,N=N,N.crit=N.crit))
}
