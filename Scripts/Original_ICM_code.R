##### Incidental Catch MOdel

### incorporates code from POR.sim.R - NOTE that the u.calc and f.crit.calc functions in the original code were modified
## and MC simulation_matrix_POR_v1.R

### H.Bowlby April 2020

setwd('c:/mydocs/ICCAT/POR2020')

## required libraries for Enric's code
library(MASS)
library(popbio)
library(tidyverse)
library(grid)
library(truncnorm)

library(ggplot2)


#################################################################

## set up functions 

#################################################################
### note that the parms.calc function generates the information required for the main simulation model
## this function is simplified from standard.SHK - runs faster than calling the full standard.SHK function

### the standard.SHK function called by the plotting code is specifically used to do calculations
## relative to reference points. By separating the main simulations from the evaluation components
## we can run one and not the other if need be - speeds up overall evaluation.
## the functions could be better integrated with more time.


#####################################
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


u.calc<-function(age.mat,max.age,r.M,sel,removals,N)
{
  #this function finds the value of f (of the vulnerable, not total, population taken each year)
  # note, because the age of selectivity is so low, this is similar to the total.
  #this is based on an equilibrium age structure assumption
  #selectivity is knife-edged at sel and 0 after maturation
  #this calculates u only for the known removals component
  
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

parms.calc<-function(repro.cycle)
{
  ## use to get max age, age at maturity, mx and M for each iteration of the simulation
  e.cor <- matrix(c(1, -0.99, -0.88,
                    -0.99, 1, 0.94,
                    -0.88, 0.94, 1),
                  ncol = 3, dimnames = list(c("Linf", "k", "t0"))) 
  colnames(e.cor) <- c("Linf", "k", "t0")
  
  # Correlation matrix for population maturity ogive
  mat.cor <- matrix(c(1, -0.99,
                      -0.99,  1),
                    ncol = 2, dimnames = list(c("Intercept", "slope")))
  colnames(mat.cor) <- c("Intercept", "slope")

  #browser()
  # VBGF parameter
  Linf.se <- 13.36
  k.se <- 0.007 
  t0.se <-  0.474 
  
  # covariance matrix that uses the correlation matrix of VBGF parameters
  e.cov <- e.cor*as.matrix(c(Linf.se, k.se, t0.se))%*%t(as.matrix(c(Linf.se, k.se, t0.se))) # The covariance matrix
  VB_pars <- mvrnorm(1, mu = c(309.8, 0.061, -5.90), Sigma = e.cov)
  Linf <- as.numeric(VB_pars["Linf"])
  k <- as.numeric(VB_pars["k"])
  t0 <- as.numeric(VB_pars["t0"])
  
  #Maturity ogive
  Mat_intercept.se  <-  1.6793 
  Mat_slope.se <- 0.1179 
  Mat.cov <- mat.cor*as.matrix(c(Mat_intercept.se, Mat_slope.se))%*%t(as.matrix(c(Mat_intercept.se, Mat_slope.se)))
  Mat_pars <- mvrnorm(1, mu = c(-10.2899, 0.7299), Sigma = Mat.cov)
  Mat_intercept <- as.numeric(Mat_pars["Intercept"])
  Mat_slope <-as.numeric(Mat_pars["slope"])
  
  # litter size which is transformed to number of female pups per female per year.
  litter.size <- rnorm(1,4,0) # average litter size (both sexes).
  gest.period <- 9/12 # gestation period in years
  #repro.cycle <- 2 # reproductive periodicity in years
  n.fem.pups <- litter.size/2# average litter size divided by two to account for just females. 
  n.pups.year <- n.fem.pups/repro.cycle #Number of female pups per year
  max.age.lower.bound <- 25
  tmax <- 7*log(2)/k
  max.age.upper.bound <- tmax
  
  max.age <- round(runif(1,max.age.lower.bound,max.age.upper.bound))
  Age <- 0:max.age # age vector for age specific values
  
  # Length-at-age from multivariate VBGF dist
  Lt <- Linf*(1-exp(-k*(Age-t0)))
  # Weight at age from W-L relationship
  W_a <- 1.482E-5
  W_b <- 2.9641
  Wt <- W_a*Lt^W_b
  Maturity_ogive <- round(1/(1+exp(-(Age*Mat_slope+Mat_intercept))),2)
  Fecundity_ogive <- round(1/(1+exp(-((Age-gest.period)*Mat_slope+Mat_intercept))),2)
  age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
  age_at_first_repro <- round(gest.period+age_at_maturity) # for adult elasticity calculations
  
  Maturity_ogive <- round(1/(1+exp(-(Age*Mat_slope+Mat_intercept))),2)
  age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
  
  Fecundity_ogive <- round(1/(1+exp(-((Age-gest.period)*Mat_slope+Mat_intercept))),2)
  mx <- Fecundity_ogive * n.pups.year
  
  #===========================================================================
  # Mortality estimation
    # Jensen methods  
  Jensen.mat <- 1.65/(age_at_maturity)
  Jensen.k <- 1.5*k
    #Peterson & Wroblewski method
  Peterson <- 1.92*(Wt*1000)^-0.25
    #Then et al. methods that update the old Hoenig and Pauly methods
  Then_hoenig <- 4.899*tmax^-0.916
  Then_pauly <- 4.118*k^0.73*Linf^-0.33
    # List all M estimates in a table
  M.table <-  cbind(Jensen.k, Jensen.mat, Then_hoenig, Then_pauly, Peterson)
    #"sample" above randomly gets 1 value of mortality from those listed in the table
  #"min" for maximum survivorship from all methods - note, taking a random sample from one row
  Mortality <- M.table[,min(x=ncol(M.table),size=1)]
  
  return(list(r.M=Mortality,mx=mx,max.age=max.age,age.mat=age_at_maturity))
  
  }
  

# Standard Leslie matrix models with MC simulations------------------------------------------------------------

# Runs matrix analysis with a post-breeding census
## n is number of MC simulations
## Age at first capture (AAFC) and Fishing mortality (F) can be specified
## repro.cycle refers to periodicity: 1 = annual, 1.5 = annual + biennial, 2 = biennial

standard.SHK <- function(n,AAFC,F,repro.cycle){
  
  
  #setup windows progress bar
  pb <- winProgressBar(title="Demography progress", label="0% done", min=0, max=100, initial=0)
  
  # correlation matrix for VBGF parameters
  #=====================================================
  # Correlation matrices for multivariate normal distributions
  #=====================================================
  
  # correlation matrix for VBGF parameters
  
  e.cor <- matrix(c(1, -0.99, -0.88,
                    -0.99, 1, 0.94,
                    -0.88, 0.94, 1),
                  ncol = 3, dimnames = list(c("Linf", "k", "t0"))) 
  colnames(e.cor) <- c("Linf", "k", "t0")
  
  # Correlation matrix for population maturity ogive
  mat.cor <- matrix(c(1, -0.99,
                      -0.99,  1),
                    ncol = 2, dimnames = list(c("Intercept", "slope")))
  colnames(mat.cor) <- c("Intercept", "slope")
  #######################
  
  
  # Create empty containers for the simulations to fill
  lambda.table <-  vector(length = n)
  r.table <-  vector(length = n)
  R0.table <- vector(length = n)
  G.table <- vector(length = n)
  s0.table <- vector(length = n)
  alpha.hat.table <- vector(length = n)
  steepness.table <- vector(length = n)
  R.table <- vector(length = n)
  SPRmer.table <- vector(length = n)
  stable.age.table <-  matrix(ncol = n, nrow = 81)
  repro.value.table <-  matrix(ncol = n, nrow = 81)
  
  #Elasticity tables
  juv.ratio.table <- vector(length = n)
  adult.ratio.table <- vector(length = n)
  juv.elast.table <- vector(length = n)
  fert.elast.table <- vector(length = n)
  adult.elast.table <- vector(length = n)
  
  
  # main loop that runs the Monte Carlo simulations
  for(i in 1:n){
    # Initiate progress bar
    info <- sprintf("%d%% done", round((i/n)*100))
    setWinProgressBar(pb, i/(n)*100, label=info)
    
    # VBGF parameter
    Linf.se <- 13.36
    k.se <- 0.007 
    t0.se <-  0.474 
    
    # covariance matrix that uses the correlation matrix of VBGF parameters
    e.cov <- e.cor*as.matrix(c(Linf.se, k.se, t0.se))%*%t(as.matrix(c(Linf.se, k.se, t0.se))) # The covariance matrix
    VB_pars <- mvrnorm(1, mu = c(309.8, 0.061, -5.90), Sigma = e.cov)
    ## note these VB parameters are identical to those in Campana et al. 2010
    Linf <- as.numeric(VB_pars["Linf"])
    k <- as.numeric(VB_pars["k"])
    t0 <- as.numeric(VB_pars["t0"])
    
    #Maturity ogive
    Mat_intercept.se  <-  1.6793 
    Mat_slope.se <- 0.1179 
    Mat.cov <- mat.cor*as.matrix(c(Mat_intercept.se, Mat_slope.se))%*%t(as.matrix(c(Mat_intercept.se, Mat_slope.se)))
    Mat_pars <- mvrnorm(1, mu = c(-10.2899, 0.7299), Sigma = Mat.cov)
    Mat_intercept <- as.numeric(Mat_pars["Intercept"])
    Mat_slope <-as.numeric(Mat_pars["slope"])
    
    
    # litter size which is transformed to number of female pups per female per year.
    litter.size <- rnorm(1,4,0) # average litter size (both sexes).
    gest.period <- 9/12 # gestation period in years
   # repro.cycle <- 2 # reproductive periodicity in years
    n.fem.pups <- litter.size/2# average litter size divided by two to account for just females. 
    n.pups.year <- n.fem.pups/repro.cycle #Number of female pups per year
    max.age.lower.bound <- 25
    tmax <- 7*log(2)/k
    max.age.upper.bound <- tmax
    
    max.age <- round(runif(1,max.age.lower.bound,max.age.upper.bound))
    Age <- 0:max.age # age vector for age specific values
    
    # Length-at-age from multivariate VBGF dist
    Lt <- Linf*(1-exp(-k*(Age-t0)))
    # Weight at age from W-L relationship
    W_a <- 1.482E-5
    W_b <- 2.9641
    Wt <- W_a*Lt^W_b
    Maturity_ogive <- round(1/(1+exp(-(Age*Mat_slope+Mat_intercept))),2)
    Fecundity_ogive <- round(1/(1+exp(-((Age-gest.period)*Mat_slope+Mat_intercept))),2)
    age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
    age_at_first_repro <- round(gest.period+age_at_maturity) # for adult elasticity calculations
    
    
    #===========================================================================
    # Mortality estimation
    
    # Jensen methods  
    Jensen.mat <- 1.65/(age_at_maturity)
    Jensen.k <- 1.5*k
    
    #Peterson & Wroblewski method
    Peterson <- 1.92*(Wt*1000)^-0.25
    
    #Then et al. methods that update the old Hoenig and Pauly methods
    Then_hoenig <- 4.899*tmax^-0.916
    Then_pauly <- 4.118*k^0.73*Linf^-0.33
    
    # List all M estimates in a table
    M.table <-  cbind(Jensen.k, Jensen.mat, Then_hoenig, Then_pauly, Peterson)
    
    #"sample" above randomly gets 1 value of mortality from those listed in the table
    #"min" for maximum survivorship from all methods - note, taking a random sample from one row
    Mortality <- M.table[,min(x=ncol(M.table),size=1)]
    
    # Estimate Survivorship from M and F estimates. If an AAFC is set then F is only applied to ages
    # older than the AAFC - note that this will be called in the simulation function
    survivorship <- NULL
    # browser()
    for(j in 1:max.age){
      f.vec<-c(rep(0,AAFC+1),rep(F,age_at_maturity-AAFC),rep(0,max.age-age_at_maturity))  
      if(Age[j]<=AAFC-1){survivorship[j] <- exp(-(Mortality[j]))}
      else{survivorship[j] <- exp(-(Mortality[j]+f.vec[j]))}
    }
    
    # Add a zero at the end of the survivorship vector to identify that no individuals survive
    # the final age class
    sx <- append(survivorship,0,after=length(survivorship))
    
    ##fecundity at age
    mx <- Fecundity_ogive * n.pups.year
    
    # line up vector so sx and mx are offset by 1 so that eventually fx = mx*sx+1...
    a <- mx[-1]
    
    # add zero to the end so both vectors are the same length again.
    MX <- append(a,0,after=length(a))
    
    # Create top row of matrix
    fx <- MX*sx 
    
    # create a matrix entirely of zeros
    Sm <- matrix(rep(0),nrow=length(survivorship),ncol=length(survivorship)) 
    
    #insert survival along the diagonal
    diag(Sm) <- survivorship 
    
    # create extra row to designate 0 values for final age class
    col <- rep(0,length(survivorship)) 
    Sm <- cbind(Sm,col,deparse.level = 0)
    
    # join fecundity vector to survival matrix without row names (deparse.level=0)
    matrix <- rbind(fx,Sm,deparse.level = 0)
    
    ############ Matrix projection ##########
    
    # calculate lambda, convert to r for other calculations and list Lambda in a table for future analysis
    lambda <- eigen.analysis(matrix)$lambda
    lambda.table[i] <-  lambda 
    r <- log(lambda)
    r.table[i] <- r
    
    #convert sx to lx values for R0 and G calculations
    lx <- NULL
    for(j in 3:length(Age)){
      lx[1] <- 1
      lx[2] <- survivorship[1]
      lx[j] <- survivorship[j]^(j-1)
      lx[max.age+1] <- 0
    }
    
    # Net reproductive rate 
    R0 <- sum(mx*lx)
    R0.table[i] <- R0 
    
    #Generation time
    G <- sum(Age*(exp(-r*Age)*(mx*lx)))
    G.table[i] <- G 
    
    #Maximum lifetime reproductive rate, steepness, position of inflection point, and SPRmer
    #First year survivorship (s0)
    s0 <- lx[2]
    s0.table[i] <- s0
    alpha.hat <- s0*R0
    alpha.hat.table[i] <- alpha.hat
    steepness <- alpha.hat/(4+alpha.hat)
    steepness.table[i] <- steepness
    R <- 0.633-(0.187*log(r*G))
    R.table[i] <- R
    SPRmer <- 1/sqrt(alpha.hat)
    SPRmer.table[i] <- SPRmer
    
  } # End of loop for individual simulations
  
  #=========================================================================
  # Analyse Monte Carlo outputs
  
  # results with na's removed if any simulations produced errors
  lambda.table <- na.omit(lambda.table)
  r.table <- na.omit(r.table)
  G.table <- na.omit(G.table)
  R0.table <- na.omit(R0.table)
  s0.table <- na.omit(s0.table)
  
  # median and quantiles of demographic results
  lambda.median <-  median(lambda.table, na.rm = T) 
  lambda.quantile <- quantile(lambda.table, c(.025, .975),na.rm = T)
  lambda.results <- as.vector(append(lambda.median,lambda.quantile))
  r.median <-  median(r.table, na.rm = T) 
  r.quantile <- quantile(r.table, c(.025, .975),na.rm = T)
  r.results <- as.vector(append(r.median,r.quantile))
  R0.median <-  median(R0.table, na.rm = T)
  R0.quantile <- quantile(R0.table, c(.025, .975),na.rm = T)
  R0.results <- as.vector(append(R0.median,R0.quantile))
  G.median <-  median(G.table, na.rm = T)
  G.quantile <- quantile(G.table, c(.025, .975),na.rm = T)
  G.results <- as.vector(append(G.median,G.quantile))
  s0.median <-  median(s0.table, na.rm = T)
  s0.quantile <- quantile(s0.table, c(.025, .975),na.rm = T)
  s0.results <- as.vector(append(s0.median,s0.quantile))
  alpha.hat.median <-  median(alpha.hat.table, na.rm = T)
  alpha.hat.quantile <- quantile(alpha.hat.table, c(.025, .975),na.rm = T)
  alpha.hat.results <- as.vector(append(alpha.hat.median,alpha.hat.quantile))
  steepness.median <-  median(steepness.table, na.rm = T)
  steepness.quantile <- quantile(steepness.table, c(.025, .975),na.rm = T)
  steepness.results <- as.vector(append(steepness.median,steepness.quantile))
  R.median <-  median(R.table, na.rm = T)
  R.quantile <- quantile(R.table, c(.025, .975),na.rm = T)
  R.results <- as.vector(append(R.median,R.quantile))
  SPRmer.median <-  median(SPRmer.table, na.rm = T)
  SPRmer.quantile <- quantile(SPRmer.table, c(.025, .975),na.rm = T)
  SPRmer.results <- as.vector(append(SPRmer.median,SPRmer.quantile))
  
  # Create result tables with titles to be printed by the function
  Results <- rbind(lambda.results,r.results,G.results,R0.results,s0.results,alpha.hat.results,steepness.results,SPRmer.results,R.results)
  colnames(Results) <- c("Median",".025",".975")
  rownames(Results) <- c("Lambda",
                         "rmax",
                         "Generation Time",
                         "Net Repro Rate",
                         "Age-0 Survivorship",
                         "Alpha hat",
                         "Steepness",
                         "SPR_MER",
                         "R")
  close(pb)
  
  return(list(Results=Results,M.table=M.table,mx.cortes=mx,M.cortes=Mortality,R0.table=R0.table,r.table=r.table,sel=AAFC,age.mat=age_at_maturity,alpha.hat.table=alpha.hat.table))
  
}# End of function










########################################################

### start analysis

########################################################
## to actually run the analysis, load in all of the functions, including por.sim() below.
## then call: por.sim.result<-por.sim(n.sims=??, repro.cycle=1 or 2)
## then call: plot.por.sim()

### note that the number of simulations, the age at first selectivity to the fishery and 
# the removals scenarios that shoudl be evaluated need to be set in the por.sim() 
# and standard.SHK functions

### NOTE that I left the number of iterations hardcoded into standard.SHK so that the main 
## simulation model coudl be run with a different number of iterations than the life hsitory evaluation.






## read in removals data from Task 1 - Provided by Carlos JUNE 18 - incorporates modifications to address
## changes to catch series from the meeting - estimation of recent catches .
# see associated details - t1nc_POR-all_20200605_v3_FINALcatchSERIES.xlsx
calc.N<-read.csv('calc.N.NWupdated.csv')
## note - calculation of N in the data file still uses mean size from Canada 39.34

## note the assumption is a median length of 149 cm and median weight of 39.34kg to gets numbers from biomass
## from 147 and W_a and W_b above, get 39.35 rather than 
nafo<-read.csv('calc.N.NAFO.csv')
#### multiple options for removals series.

## option 1: using size frequency all years combined. 
### add size distribution from Rui - Provided June 3, 2020
size<-read.csv('NW.size.csv')
### use growth parameters above to transform this directly into biomass.
W_a <- 1.482E-5
W_b <- 2.9641

size$wt<-W_a*size$Size_cm_FL^W_b

junk<-hist(size$wt,breaks=43,plot=F)  ### every 5cm size class
##proportion of animals in each weight bin
split<-junk$counts/sum(junk$counts)

## biomass * split / weight in each bin
N.ratio<-c()
for (i in 1:length(calc.N$por.kg))
{
  ## weighted average in 5 cm length bins
  N.ratio[i]<-sum(calc.N$por.kg[i]*split/junk$breaks[2:length(junk$breaks)])
}


## option 2: using median from NW size frequency, all years combined
N.med<-calc.N$por.kg/median(size$wt)

## year-specific length-frequency or median size.

N.med.yr<-c()
N.ratio.yr<-c()
year<-calc.N$year

for (i in 1:length(year))
{
  if(i < 35)
  {
    junk<-hist(size$wt[size$Year==1994],breaks=43,plot=F)  ### every 5cm size class
    ##proportion of animals in each weight bin
    split<-junk$counts/sum(junk$counts)
       ## biomass * split / weight in each bin
       N.ratio.yr[i]<-sum(calc.N$por.kg[i]*split/junk$breaks[2:length(junk$breaks)])
       N.med.yr[i]<-calc.N$por.kg[i]/median(size$wt[size$Year==1994])
  }
  else
  {
    junk<-hist(size$wt[size$Year==year[i]],breaks=43,plot=F)  ### every 5cm size class
    ##proportion of animals in each weight bin
    split<-junk$counts/sum(junk$counts)
    ## biomass * split / weight in each bin
    N.ratio.yr[i]<-sum(calc.N$por.kg[i]*split/junk$breaks[2:length(junk$breaks)])
    N.med.yr[i]<-calc.N$por.kg[i]/median(size$wt[size$Year==year[i]])
  }
}

N.options<-data.frame(calc.N,N.ratio,N.ratio.yr,N.med,N.med.yr)

write.csv(N.options,'options.for.biomass.to.number.animals.NW.csv')

# average removals for the most recent 5 years (2014-2018) are 39 mt = 39000kg. 
# mean size of removals = 147cm for Canada; used 149 because listed in Figure 3 of Campana et al. 2012
# 149cm = average weight of 39.34 kg
## gives an average removals of 991 animals (39000/39.34) in 2014-2018: used below to represent recent time period in forward projections
# HOWEVER, average removals thorughout the time series are 1789mt; gives 45,475 animals (using Task 1 series)


  por.sim<-function(n.sims,q.ext=500,sigma=0,u.long=F,repro.cycle)
  {
    #assign output to "por.sim.result" which is an object used by the plotting function below 
    
    set.seed(204) #so that different scenarios are directly comparable

    #setup windows progress bar
    pb <- winProgressBar(title="Simulation progress", label="0% done", min=0, max=100, initial=0)
    
    years<-c(N.options$year,2019)
    n.years<-length(years)
    forward.n.years<-52 #first year gets n.end, second year gets assumed removals and then project for 50 years
    Pop<-matrix(rep(0,n.years*n.sims),n.sims,n.years)     
    Forward.pop<-matrix(rep(0,forward.n.years*n.sims),n.sims,forward.n.years)
    Pop.vec<-rep(0,n.years)
    ext.time<-rep(110,n.sims)   #NOT USED: values will be overwritten with values from 1 to 100 if extinction happens. 110 means>100 years 
    Forward.pop.vec<-rep(0,forward.n.years)
    r.nomort.vec<-rep(0,n.sims)  ## productivity in the absence of fishing
    r.allmort.vec<-rep(0,n.sims) #NOT USED: includes observed and other mortality
    F.crit.vec<-rep(0,n.sims)
    N.crit.vec<-rep(0,n.sims)
    mean.u.long.vec<-rep(0,n.sims) # NOT USED
    mean.u.short.vec<-rep(0,n.sims) # NOT USED

    N.crit.list<-list(rep(NA,length(n.sims)))
    u.future.list<-list(rep(NA,length(n.sims)))
    F.future.list<-list(rep(NA,length(n.sims)))
    r.scenarios<-list(rep(NA,length(n.sims)))
    Forward.pop.scenarios<-list(rep(NA,length(n.sims)))
    
    B1.vec<-rep(0,n.sims)
    
    removals<-c(N.options$N.med,0) 
    #removals<-removals*2
    mean.removals<-mean(N.options$N.med[N.options$year<2010])  ## NOT USED: 2009 was the year that removals of porbeagle REALLY dropped.
    
    sel<-1 ## basically vulnerable to the fishery right after birth. year 1 onwards.
    r.cutoff<-0.2 #NOT USED. note that this is extremely high
    K<-20000000  # carrying capacity of 20 million for forward projections
   
    age.mat.vec<-rep(0,n.sims)
    max.age.vec<-rep(0,n.sims)
    mx.list<-list(rep(NA,length(n.sims)))
    r.M.list<-list(rep(NA,length(n.sims)))
    
    #create parameters -  consistent with estimation method used in standard.SHK below
    ## note that lotka.r function just removes necessity of developing the survival by age matrix
    ## and is used here as a short-cut. Parameter estimates are within three decimal places of rmax
    # from standard.SHK
    
   for (i in 1:n.sims)
   {junk<-parms.calc(repro.cycle)
    age.mat.vec[i]<-junk$age.mat
    max.age.vec[i]<-junk$max.age
    mx.list[[i]]<-junk$mx
    r.M.list[[i]]<-junk$r.M
   }
     
    N.end<-c(sample(100000:200000,size=n.sims/2,replace=T),sample(200000:300000,size=n.sims/2,replace=T)) 
    #N.end<-c(sample(150000:200000,size=n.sims/2,replace=T),sample(200000:250000,size=n.sims/2,replace=T)) 
    #logic here is approx. 200000 sharks currently - from SCA in Campana et al. 2010. 
    # sampling done this way so that half are guaranteed above and below 
 
    for(i in 1:(1*n.sims))   #2 is a patch to get n.sims values for r that are <r.cutoff 
    {
      #print(junk$par)
      junk2<-F.crit(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],sel,mean.removals)
      #print(junk2$spr.f0)
      junk3<-lotka.r(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],0,sel) #no human induced mortality  
      F.crit.vec[i]<-junk2$f      
      N.crit.vec[i]<-junk2$N.crit
      r.nomort.vec[i]<-junk3$par
            }
   
    ############################################################
    
    #then do backward population projections
    
    for(i in 1:n.sims)
    {
      # Initiate progress bar
      info <- sprintf("%d%% done", round((i/n.sims)*100))
      setWinProgressBar(pb, i/(n.sims)*100, label=info)
      
      Pop.vec[n.years]<-N.end[i]
      for(y in 1:(n.years-1))
      {
        Pop.vec[n.years-y]<-(Pop.vec[n.years-y+1]+removals[n.years-y])/exp(r.nomort.vec[i])  ## productivity in the absense of fishing
      }
      
      Pop[i,]<-Pop.vec
      #print(Pop.vec)
      B1.vec[i]<-lm(log(Pop.vec)~years)$coef[2]
      
      #calculate historical exploitation - note: as is, this code is not used right now.
      junk<-c()
      
      for(y in 1:n.years)
      {
        # first calc is relative to observed removals.
        junk[y]<-u.calc(age.mat.vec[i],max.age.vec[i],r.M.list[[i]],sel,removals[y],Pop.vec[y])

      }

      mean.u.long.vec[i]<-mean(junk[(length(junk)-15):length(junk)-11])  # NOT USED: years are 2004 - 2008; 
      mean.u.short.vec[i]<-mean(junk[(length(junk)-10):length(junk)-1])  # NOT USED; years are 2009-2018

      ### this is what is used to evaluate future removals.
      ## calculate exploitation relative to current pop size for hypothetical future removal scenarios      
      junkx<-c()
      rem.future<-c(0,seq(1000,24000,by=1000))  ### added zero removals as a scenario
      scenarios=length(rem.future)
      
      for (k in 1:length(rem.future))
      {
        ## this list is relative to hypothetical levels we are interested in testing
        junkx[k]<-u.calc(age.mat.vec[i],max.age.vec[i],r.M.list[[i]],sel,rem.future[k],Pop.vec[length(Pop.vec)])
      }
      u.future.list[[i]]<-junkx
      F.future.list[[i]]<--log(1-junkx)

           }
      
    ### evaluated removals table    
   #then do forward projections:
    for(i in 1:n.sims)
    {
     
       #first recalc r with fishing mortality at the different levels of removals  
      xxx<-as.vector(u.future.list[[i]])
      
      x1<-c()
      x2<-c()
      temp<-c()
      junkxx<-c()
     for(k in 1:length(xxx))
     {
      x1<-lotka.r(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],xxx[k],sel)
      x2[k]<-x1$par
      temp<-F.crit(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],sel,xxx[k]*Pop.vec[n.years])
      junkxx[k]<-temp$N.crit
     }
     r.scenarios[[i]]<-x2
     N.crit.list[[i]]<-junkxx
     
      ## note - this is not used. 
     ## hold-over code from basking shark - where we wanted to assume future removals based on historical
      if(u.long==T)
      {
        junk1<-lotka.r(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],mean.u.long.vec[i],sel)  
      }
      else
      {
        junk1<-lotka.r(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],mean.u.short.vec[i],sel)  
      }
      
      r.allmort.vec[i]<-junk1$par 
      
      ## note that this also isn't used. Incorporates autocorrelated deviates to annaul values for r
      ## in future projections. I would need to make a matrix that the deviates were applied over
      ## and I just didn't have time.
      #random deviate vectors 
      log.dev.vec<-rep(0,forward.n.years) #sigma=0 case
      r.forward<-rep(r.allmort.vec[i],forward.n.years)
      if(sigma>0)  #overwrite above if sigma>0 
      {
        w<-rnorm(forward.n.years,0,1) 
        log.dev.vec[1]<-w[1]*sigma 
        r.forward[1]<-r.allmort.vec[i]+log.dev.vec[1]#-sigma^2/2
        #if(r.forward[1]>r.cutoff){r.forward[1]<-r.cutoff} #sets an upper bound on r in the sims #removed not logical 
        
        for (y in 2:forward.n.years)
        {
          log.dev.vec[y]<-w[y]*sigma
          r.forward[y]<- r.allmort.vec[i]+log.dev.vec[y]#-sigma^2/2
          if(r.forward[y]>r.cutoff){r.forward[y]<-r.cutoff} #sets an upper bound on r in the sims}
        }    
      }
      #  print(summary(r.forward))
#browser()      
      R.ave<-mean(N.options$N.med[seq(length(year)-2,length(year))])
 
      #then project forward
      ## 2019
      Forward.pop.vec[1]<-N.end[i]*exp(r.nomort.vec[i])-R.ave 
      r.first<-r.nomort.vec[i]
      
      #print(r.first)
      ## I am being lazy here by not incorporating future variability in r
      xx<-as.vector(r.scenarios[[i]])
      
      Forward.pop.mat<-matrix(rep(0,forward.n.years*length(xx)),length(xx),forward.n.years)
      Forward.pop.mat[,1]<-N.end[i]*exp(r.nomort.vec[i])-R.ave 
       
      for(y in 1:(forward.n.years-1))
      {

        if(y==1)
        {
          # 2020
          Forward.pop.vec[y+1]<-Forward.pop.vec[y]*exp(r.first)*(1-(Forward.pop.vec[y]/K))-R.ave
          Forward.pop.mat[,y+1]<-(Forward.pop.mat[,y]*exp(r.first))*(1-(Forward.pop.mat[,y]/K))-R.ave
          
        }
        
        #2021 onwards
        else
        {
          Forward.pop.vec[y+1]<-(Forward.pop.vec[y]*exp(r.forward[y]))*(1-(Forward.pop.vec[y]/K))
          
          ### old code: used to determine an extinction threshold for species at risk questions.
          ## note, I have left q.ext initialized in the rest of the simulation.
          #if(Forward.pop.vec[y]<q.ext)
          #{
          #  Forward.pop.vec[y+1]<-0
          #  if(ext.time[i]>101){ext.time[i]<-y}
          #}
          
          ###project forward - note that these don't have variability in the annual values of future r
          Forward.pop.mat[,y+1]<-(Forward.pop.mat[,y]*exp(xx))*(1-(Forward.pop.mat[,y]/K))
          
        }
       
      }
      Forward.pop[i,]<-Forward.pop.vec
      Forward.pop.scenarios[[i]]<-Forward.pop.mat
     #print(Forward.pop.mat)
    } #end simulation
    close(pb)
    
    #=========================================================================
    # Summarize projections 
    
    ### note that a lot of these summary values aren't used in the current assessment. 
    ## they are hold-overs from previous use of the code.
    end.pop.size<-Forward.pop[,forward.n.years]
    
    r.nomort.summary<-quantile(r.nomort.vec,c(0.1,0.5,0.9))  #no human induced mortality
    r.allmort.summary<-quantile(r.allmort.vec,c(0.1,0.5,0.9))  #includes both "other" and known removals 
    mean.u.long.summary<-quantile(mean.u.long.vec,c(0.1,0.5,0.9))  
    mean.u.short.summary<-quantile(mean.u.short.vec,c(0.1,0.5,0.9))  
    F.crit.summary<-quantile(F.crit.vec,c(0.1,0.5,0.9))   
    N.crit.summary<-quantile(N.crit.vec,c(0.1,0.5,0.9))   
    N.1950.summary<-quantile(Pop[,1],c(0.1,0.5,0.9))   
    N.2019.summary<-quantile(Forward.pop[,2],c(0.1,0.5,0.9))   
    N.2119.summary<-quantile(Forward.pop[,forward.n.years],c(0.1,0.5,0.9))   
    
    proportion.declining<-length(B1.vec[B1.vec<0])/length(B1.vec)
    end.pop.summary<-quantile(end.pop.size,c(0.1,0.5,0.9)) 
    prop.ext<-length(end.pop.size[end.pop.size<q.ext])/length(end.pop.size)
    
    return(list(N.crit.list=N.crit.list,u.future.list=u.future.list,F.future.list=F.future.list,r.scenarios=r.scenarios,rem.future=rem.future,Forward.pop.scenarios=Forward.pop.scenarios,scenarios=scenarios,n.sims=n.sims,years=years, q.ext=q.ext,sigma=sigma,mean.removals=mean.removals,
                Pop=Pop,Forward.pop=Forward.pop,forward.n.years=forward.n.years,r.allmort.vec=r.allmort.vec,r.nomort.vec=r.nomort.vec,
                B1.vec=B1.vec,mean.u.long.vec=mean.u.long.vec,r.nomort.summary=r.nomort.summary,
                r.allmort.summary=r.allmort.summary,F.crit.summary=F.crit.summary,N.crit.summary=N.crit.summary,N.1950.summary=N.1950.summary,
                mean.u.long.summary=mean.u.long.summary,mean.u.short.summary=mean.u.short.summary,N.2019.summary=N.2019.summary,
                N.2119.summary=N.2119.summary,proportion.declining=proportion.declining,end.pop.summary=end.pop.summary,prop.ext=prop.ext,ext.time=ext.time))
    
  } #end function
  ###############################################################
  ###############################################################
  
  plot.por.sim<-function(n.sims,repro.cycle)
  {
      
    #
    pop<-por.sim.result$Pop/1000
    n.sims<-por.sim.result$n.sims
    scenarios<-por.sim.result$scenarios
    forward.years<-seq(1:por.sim.result$forward.n.years)+2018
    rem.future<-por.sim.result$rem.future
    
    r<-por.sim.result$r.nomort.vec
    B1<-por.sim.result$B1.vec
    ann.change<-por.sim.result$ann.change
    #par(mfrow=c(1,1),las=1,omi=c(1,1,0.5,0.25),mar=c(3,3,1,1)) 

    ### plot of the medians of the projections at different removal scenarios    
    xxx<-do.call(rbind.data.frame,por.sim.result$Forward.pop.scenarios)
    xxx$scenario<-rep(1:scenarios,n.sims) ## check

    med.scen<-xxx %>% group_by(scenario) %>% summarize_at(vars(V1:V52),median)
    ## need long data format for plotting
    values<-names(med.scen[2:length(med.scen)])

    
    ### determine median F for input into standard.SHK below.
    xxx<-data.frame(do.call(rbind,por.sim.result$F.future.list))
    #names(xxx)<-as.character(por.sim.result$rem.future)
    
   # med.F<-xxx %>% summarize_all(median) 
    ## note that this is lazy - you inappropriately calculated medians for the majority of columns
   
    ####### get r value at diff levels F
    x0<-standard.SHK(n.sims,1,0,repro.cycle)

    ## note that there are NaN produced as r estimate goes negative. (warnings)
   # x1<-list()
    # for(i in 1:length(med.F))
    #{
    # temp<-standard.SHK(n.sims,1,as.numeric(med.F[i]),repro.cycle)  
    # x1[[i]]<-temp$Results[2,]
    #}
 #x2<-do.call(rbind,x1)
    
## summary table of how R values change with different fishing scenarios.
  #  r.decline<-data.frame(rbind(x0$Results[2,],x2))
   # rownames(r.decline)<-c('none',as.character(por.sim.result$rem.future))
    #r.decline$biomass<-c(0,as.numeric(por.sim.result$rem.future)*39.34/1000)
    #write.csv(r.decline, paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/r.decline.csv'))

write.csv(x0$Results,paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/parameter.vals.csv'))
        
    ## get every 5th year on scale
    every_nth = function(n) {
      return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
    }
    labs<-forward.years[seq(1,length(forward.years),5)]

    #pdf(paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.pdf'),onefile = TRUE)
    
### plots of future projections under different evaluated fishing scenarios.
    scen.plot<-gather(med.scen,year,abundance,all_of(values),factor_key = T)
    labs2<-as.character(rem.future)
 p1<-ggplot(scen.plot,aes(year,abundance/1000,group=scenario,colour=as.factor(scenario)))+geom_line()
 p1<-p1+scale_colour_discrete(name='Median Removals',labels = c(labs2))+labs(y='Total Abundance (1000s)',x='Year')
 p1<-p1+scale_x_discrete(breaks=every_nth(n=5),labels=c(labs)) #+theme(axis.text.x = element_text(angle = 90))
p1

 png(file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.future.png'),units='in', width=10, height=6,res=800)
 print(p1)
dev.off()
#browser()
     ## historical trajectory 
    years<-c(N.options$year,2019)
    abund<-c()
    a.90<-c()
    a.10<-c()
    for (i in 1:length(years))
    {abund[i]<-median(por.sim.result$Pop[,i]) 
    a.90[i]<-quantile(por.sim.result$Pop[,i],c(0.9))
    a.10[i]<-quantile(por.sim.result$Pop[,i],c(0.1))}

    png(file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.abund.png'),units='in', width=12, height=6,res=800)
    par(mfrow=c(2,1),mar=c(4,4,0,4))
    plot(years, abund/1000,type='l',lwd=2,xlab='',ylab='Abundance (thousands)',ylim=c(0,1.4*max(abund/1000))) 
    lines(years,a.90/1000,lty=5)
    lines(years, a.10/1000,lty=5)
    plot(N.options$year,N.options$N.med/1000, type='l',xlab=c('Year'),ylab=c('Removals (thousands)'))
  dev.off()

  print('% decline from historical maximum to present')
  print((max(abund)-abund[length(abund)])/max(abund))
  print('max to min abundance decline rate')
  print((max(abund)-min(abund))/max(abund))
  print(data.frame(years,abund))
  
#browser()  
    ### determine if overfished - allowing variability in SPRmer.
    alpha.hat.x0<-x0$alpha.hat.table
    SPRmer.S0<-(sqrt(alpha.hat.x0)-1)/(alpha.hat.x0-1)
    scaled.decline<-abund[length(abund)]/abund[1]
    ind<-scaled.decline/median(SPRmer.S0)
    upper<-1#-(median(x0$M.cortes[1]))
    hist(ind,breaks=30,main='Variability in derivation of SPRmer',xlab='Critical value',xlim=c(0.2,1))
    abline(v=0.5,col='blue')
    abline(v=upper,col='red')

  values<-c(length(ind[ind<0.5])/length(ind),length(ind[ind<upper])/length(ind))    
  crit.levels<-c('p<0.5','p<1')


    ### determine if overfished - allowing variability in population estimates.
    alpha.hat.x0<-x0$Results[6,1]
    SPRmer.S0<-(sqrt(alpha.hat.x0)-1)/(alpha.hat.x0-1)
    scaled.decline<-por.sim.result$Pop[,ncol(por.sim.result$Pop)]/por.sim.result$Pop[,1]  ### note that this is specific to the N. Atlantic right now.
    ind<-scaled.decline/median(SPRmer.S0)
    upper<-1#-(median(x0$M.cortes[1]))  ## this corresponds to MSY proxy, rather than MSY proxy minus M

  png(file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.overfished.png'),units='in', width=9, height=6,res=800)
    hist(ind,breaks=50, xlim=c(0.2,1),main='',xlab='Critical value')
    abline(v=0.5,col='blue')
    abline(v=upper,col='red')
   dev.off()
    
    values2<-c(length(ind[ind<0.5])/length(ind),length(ind[ind<upper])/length(ind)) ## proportion of sims that ARE overfished (i.e. below p reference point)   
    overfishing.eval<-data.frame(crit.levels,var.in.SPRmer=values,var.in.Popsize=values2)
    write.csv(overfishing.eval,paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/overfishing.eval.csv'))

    ### table of reaching recovery target under different forward fishing scenarios
    ## first, calculate the values for the decline in the index:
    
    trajectories<-do.call(rbind.data.frame,por.sim.result$Forward.pop.scenarios)
    trajectories$scen<-rep(rem.future,length=n.sims)
#browser()    
    ## set up the columns to sequence
    #test<-trajectories%>%select(c(1,2,seq(3,52,5)))
    #name.breaks<-colnames(test)
    ind.decline<-trajectories/por.sim.result$Pop[,1]
    #test$scen<-rep(rem.future,length=n.sims)
## use SPRmer from above - note you can't have variability in SPRmer plus in pop size - use pop size
    ind<-ind.decline/median(SPRmer.S0)
    ind$scen<-rep(rem.future,length=n.sims)
    values3<-ind%>%group_by(scen)%>%summarize_all(~sum(.>0.5)/n.sims) ## proportion of values > 0.5 i.e. not overfished
    values3<-values3%>%select(c(scen,V2,V7,V12,V17,V22,V27,V32,V37,V42,V47,V52)) ## every 5 years
    names(values3)<-c("Scenario",seq(2020,2070,5)) ## NOTE this is specific to the projection timeframe of 50 years + 2019, 2020
    values4<-ind%>%group_by(scen)%>%summarize_all(~sum(.>upper)/n.sims)
    values4<-values4%>%select(c(scen,V2,V7,V12,V17,V22,V27,V32,V37,V42,V47,V52))
    names(values4)<-c("Scenario",seq(2020,2070,5))
 ## output the summarized tables
    write.csv(values3,paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/SPRmer.0.5.csv'))
    write.csv(values4,paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/SPRmer.M.csv'))
    
### biomass at SPRmer.S0 ref point = 353,000  
    test2<-trajectories%>%group_by(scen)%>%summarize_at(vars(V1:V52),median)
    values<-names(test2[2:length(test2)])
    unfish<-mean(por.sim.result$Pop[,1])  ### using historical maximum as unfished abundance
    target<-unfish*median(SPRmer.S0)
    labs<-forward.years[seq(1,length(forward.years),5)]

    write.csv(data.frame(unfish=unfish,target=target,SPRmer.S0=SPRmer.S0,hist.decline=(max(abund)-abund[length(abund)])/max(abund)),
              file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/biomass.ref.points.csv'))
    
    png(file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/b_bsprmer.projection.png'),units='in', width=9, height=6,res=800)
    scen.bmsy<-gather(test2,year,abundance,all_of(values),factor_key = T)
    labs2<-as.character(rem.future)
    pmsy<-ggplot(scen.bmsy,aes(year,abundance/unfish/SPRmer.S0,group=scen,colour=as.factor(scen)))+geom_line()
    pmsy<-pmsy+scale_colour_discrete(name='Median Removals',labels = c(labs2))+labs(y='B/Bsprmer',x='Year')+geom_hline(yintercept=SPRmer.S0/SPRmer.S0)
    pmsy<-pmsy+scale_x_discrete(breaks=every_nth(n=5),labels=c(labs)) #+theme(axis.text.x = element_text(angle = 90))
    print(pmsy)
    dev.off()
    
    
       
  #      years.f<-c(2019:2069)  
    ## future trajectory
   # abund<-c()
    #a.90<-c()
    #a.10<-c()
    #for (i in 1:por.sim.result$forward.n.years)
    #{abund[i]<-median(por.sim.result$Forward.pop[,i]) 
    #a.90[i]<-quantile(por.sim.result$Forward.pop[,i],c(0.9))
    #a.10[i]<-quantile(por.sim.result$Forward.pop[,i],c(0.1))}
    
    #plot(years.f, abund/1000,type='l',lwd=2,xlab='Years',ylab='Abundance (thousands)',ylim=c(0,1.2*max(abund/1000))) 
    #lines(years.f,a.90/1000,lty=5)
    #lines(years.f, a.10/1000,lty=5)

    
    
    #par(mfcol=c(2,2),las=1,omi=c(1,1,0.5,0.25),mar=c(4,4,2,2)) 
    
  #  hist(pop[,length(years)],nclass=20,cex=0.7,xlab=" ",xlim=c(50,350))
   #mtext("N (2019)",1,outer=F,cex=1.2,line=3)
    
  #  hist(pop[,1],nclass=20,cex=0.7,xlab=" ",xlim=c(0,2000))
  #  mtext("N (1950)",1,outer=F,cex=1.2,line=3)
    
   # hist(r,nclass=20,cex=0.7,xlab=" ",plot=T)
  #  mtext("Probability Density",2,outer=T,cex=1.4,line=4)
  #  mtext("r",1,outer=F,cex=1.2,line=3)
    
   # hist(B1,nclass=20,cex=0.7,xlab=" ",plot=T)
    #mtext("Probability Density",2,outer=T,cex=1.4,line=4)
    #mtext("Annual Rate of Change",1,outer=F,cex=1.2,line=3)
    #mtext("Slope",1,outer=F,cex=1.2,line=3)

    
  }
  

  
  