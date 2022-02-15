#this function finds the value of f (of the vulnerable, not total, population taken each year)
# note, because the age of selectivity is so low, this is similar to the total.
#this is based on an equilibrium age structure assumption
#selectivity is knife-edged at sel and 0 after maturation
#this calculates u only for the known removals component


u.calc<-function(age.mat,max.age,r.M,sel,removals,N)
{

  M<-r.M
  
  assign(".M",M)     # Store in expression frame
  assign(".sel",sel)
  assign(".max.age",max.age)
  assign(".age.mat",age.mat)
  assign(".N",N)
  assign(".removals",removals)
  
  
  minimise<-function(f)
  {
    
    f.vec<-c(rep(0,.sel+1),rep(f,.age.mat-.sel),rep(0,.max.age-.age.mat))  		
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
  
  junk<-nlminb(start = 0.1,obj = minimise, lower=0.0001,upper=10)
  
  f.out<-junk$par
  u.out<-1-exp(-f.out)
  
  return(u.out)
}