# Function to calculate Lotka R

# Function Arguements
# age.mat:   Age at maturity from params_calc function
# max.age:   Maximum age
# mx:        Age specific fecundity
# r.M:       Age specific Natural mortality
# u:         Exploitation rate (annual not instantaneous), currently set up to be 1 value.
# sel:       Selectivity, currently set up to be 1 value

lotka.r<-function(age.mat,max.age,mx,r.M,u,sel)
{
  # Function to solve for lotka r.
  minimise<-function(start.r)
  {
    # Get a starting estimate of lotka's r
    rx<-exp(-start.r*.x)
    # Multiply the repro contribution of each age by current value of rx and sum those up, as this approaches the correct value this approaches 1
    lotka<-sum(.lxmx*rx)
    # Then minimize the sum of squares difference of this lotka term from 1 to determine what r is
    sumsq <- sum((lotka-1)^2)
    return(sumsq)
  }

  # Instantaneous fishing mortality
  f<- -log(1-u)
  ## Set up the fishing mortality vector based on selectivity, max age, age at maturity, etch
  #note that the plus 1 is to deal with the parameterization of an age0 in Enric's code
  f.vec<-c(rep(0,sel+1),rep(f,age.mat-sel),rep(0,max.age-age.mat))
  # Survival by age including natural mortality and fishing mortality
  si<-exp(-(r.M+f.vec))

  # Now calculate survivorship by age
  x<-0:(max.age)
  lx<-0
  lx[1]<-1
  for(i in 2:(length(si)+1)) lx[i]<-lx[i-1]*si[i-1]

  # Multiply survivorship by fecundity to get reproductive contribution of each age
  lxmx<-lx[1:length(si)]*mx
  assign(".lxmx",lxmx)     # Store in expression frame
  assign(".x",x)
  # print("expected female output from one newly born female")
  #print(sum(lxmx))
  #print("mx")
  #print(mx)
  # and now solve the lotka r minimizing function to get an estimate of reproductive capacity.
  junk<-nlminb(start = 0.03,obj = minimise, lower=-1,upper=0.5)
  return(junk)
}