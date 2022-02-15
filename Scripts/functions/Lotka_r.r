# Function to calculate Lotka R

# Function Arguements
# 1: age.mat: Age at maturity
# 2: max.age: Maximum age
# 3: mx: natural mortality
# 4: r.M: Reproduction
# 5: u: Exploitation rate (annual not instantaneous)
# 6: sel: Selectivity


lotka.r<-function(age.mat,max.age,mx,r.M,u,sel) 
{
  ## ouput of parameter function to get mx and r.M vectors.
  
  f<--log(1-u)
  ## note that the plus 1 is to deal with the parameterization of an age0 in Enric's code
  f.vec<-c(rep(0,sel+1),rep(f,age.mat-sel),rep(0,max.age-age.mat))       
  
  si<-exp(-(r.M+f.vec)) 
  
  x<-0:(max.age)
  lx<-0
  lx[1]<-1
  for(i in 2:(length(si)+1))
  {lx[i]<-lx[i-1]*si[i-1]}
  
  lxmx<-lx[1:length(si)]*mx
  assign(".lxmx",lxmx)     # Store in expression frame
  assign(".x",x)
  # print("expected female output from one newly born female")
  #print(sum(lxmx))
  #print("mx")
  #print(mx)
  minimise<-function(start.r)
  {
    rx<-exp(-start.r*.x)
    lotka<-sum(.lxmx*rx)
    sumsq <- sum((lotka-1)^2) 
    return(sumsq)
  }
  junk<-nlminb(start = 0.03,obj = minimise, lower=-1,upper=0.5)
  return(junk)
}