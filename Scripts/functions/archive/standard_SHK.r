# This function is an extension of parms_calc() does everything it does, but also has a simulation
# step that is used to determine the status of the stock and get reference points (see Brooks & friends 2010 - Analytical reference points for age-structured models: application to data-poor fisheries)

# Function Arguments
# First up the VBGF parameters
# linf.mn:            The mean estimate of the Von B asymptote. Defaults to the value for Porbeagle (309.8)
# k.mn:               The mean estimate of the von B rate parameter. Defaults to the value for Porbeagle (0.061)
# t0.mn:              The mean estimate of the von B age at size 0. Defaults to the value for Porbeagle (-5.90)
# linf.se:            The standard error of the the Von B asymptote. Defaults to the value for Porbeagle (13.36)
# k.se:               The standard error of the the Von B  rate parameter. Defaults to the value for Porbeagle (0.007)
# t0.se:              The standard error of the the Von B age at size 0. Defaults to the value for Porbeagle (0.474)

# Maturity parameters
# mat.int.mn          The mean of the maturity model intercept. Defaults to the value for Porbeagle (-10.2899)
# mat.slop.mn         The mean of the maturity model slope Defaults to the value for Porbeagle (0.7299)
# mat.int.se          The standard error of the maturity model intercept. Defaults to the value for Porbeagle (1.6793)
# mat.slope.se        The standard error of the maturity model slope. Defaults to the value for Porbeagle (0.1179)

# Weight at age from W-L relationship. Note there is no uncertainty in the current W-L relationship
#W.a                  The Weight-length relationship intercept. Defaults to the value for Porbeagle (1.482E-5)
#W.b                  The Weight-length relationship exponent. Defaults to the value for Porbeagle (2.9641)


# Reproduction parameters. We'll need to tweak language here for fish.
# n.offspring         Average number of offspring which is transformed to number of female progeny per female per year.  Defaults to the value for Porbeagle (4)
# gest.period         The gestation period in years. Defaults to the value for Porbeagle (9/12)
# repro.cycle         The number of times an individual reproduces in a year

# max.age.lb         The lowest likely 'maximum age'. Defaults to the value for Porbeagle (25)
# repro.cycle        The number of times you reproduce in a year.  Defaults to 1 now
# vb.cor.mat         The correlation matrix for the von B parameters to ensure that we don't get weird parameter combinations.  Default is from PorBeagle, will have to sort out if this is fine for all.
# mat.cor.mat        The correlation matrix for the maturity parameters to ensure that we don't get weird parameter combinations.  Default is from PorBeagle, will have to sort out if this is fine for all.

# Other parameters
# n.mc                The number of Monte Carlo simulations to run
# AAFC                I think Age At First Capture. Seems like this assumes selectivity to be 1
# F                   Fishing mortality from the future removals scenarios, one for each scenario.




standard.SHK <- function(linf.mn = 309.8, k.mn = 0.061, t0.mn = -5.90, linf.se = 13.36, k.se = 0.007, t0.se = 0.474,
                         mat.int.mn = -10.2899, mat.slop.mn = 0.7299,  mat.int.se = 1.6793,mat.slope.se = 0.1179,
                         w.a = 1.482E-5, w.b = 2.9641,
                         max.age.lb = 25,n.offspring = 4,gest.period = 9/12, repro.cycle =1,
                         vb.cor.mat = 'default',mat.cor.mat = 'default', 
                         n,AAFC=0,F=0)
  {
  

  #setup windows progress bar
  pb <- winProgressBar(title="Demography progress", label="0% done", min=0, max=100, initial=0)
  
  # correlation matrix for VBGF parameters
  #=====================================================
  # Correlation matrices for multivariate normal distributions
  #=====================================================
  
  # correlation matrix for VBGF parameters
  if(vb.cor.mat == 'default')
  {
    e.cor <- matrix(c(1, -0.99, -0.88,
                      -0.99, 1, 0.94,
                      -0.88, 0.94, 1),
                    ncol = 3, dimnames = list(c("Linf", "k", "t0"))) 
    colnames(e.cor) <- c("Linf", "k", "t0")
  }
  # covariance matrix that uses the correlation matrix of VBGF parameters
  e.cov <- e.cor*as.matrix(c(linf.se, k.se, t0.se))%*%t(as.matrix(c(linf.se, k.se, t0.se))) # The covariance matrix
  
  # Correlation matrix for population maturity ogive
  if(mat.cor.mat == 'default')
  {
    mat.cor <- matrix(c(1, -0.99,
                        -0.99,  1),
                      ncol = 2, dimnames = list(c("Intercept", "slope")))
    colnames(mat.cor) <- c("Intercept", "slope")
  }
  #Maturity ogive
  Mat.cov <- mat.cor*as.matrix(c(mat.int.se, mat.slope.se))%*%t(as.matrix(c(mat.int.se, mat.slope.se)))
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
    
    # Now get the VB_pars for each run of the simulation.
    VB_pars <- mvrnorm(1, mu = c(309.8, 0.061, -5.90), Sigma = e.cov)
    ## note these VB parameters are identical to those in Campana et al. 2010
    Linf <- as.numeric(VB_pars["Linf"])
    k <- as.numeric(VB_pars["k"])
    t0 <- as.numeric(VB_pars["t0"])
    
    #Maturity ogive
    Mat_pars <- mvrnorm(1, mu = c(-10.2899, 0.7299), Sigma = Mat.cov)
    Mat_intercept <- as.numeric(Mat_pars["Intercept"])
    Mat_slope <-as.numeric(Mat_pars["slope"])
    
    
    # litter size which is transformed to number of female pups per female per year.
    litter.size <- rnorm(1,n.offspring,0) # average litter size (both sexes).
    # repro.cycle <- 2 # reproductive periodicity in years
    n.fem.pups <- litter.size/2# average litter size divided by two to account for just females. 
    n.pups.year <- n.fem.pups/repro.cycle #Number of female pups per year
    tmax <- 7*log(2)/k
    max.age.upper.bound <- tmax
    max.age <- round(runif(1,max.age.lb,max.age.upper.bound))
    Age <- 0:max.age # age vector for age specific values
    
    # Length-at-age from multivariate VBGF dist
    Lt <- Linf*(1-exp(-k*(Age-t0)))
    # Weight at age from W-L relationship

    Wt <- w.a*Lt^w.b
    Maturity_ogive <- round(1/(1+exp(-(Age*Mat_slope+Mat_intercept))),2)
    Fecundity_ogive <- round(1/(1+exp(-((Age-gest.period)*Mat_slope+Mat_intercept))),2)
    age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
    age_at_first_repro <- round(gest.period+age_at_maturity) # for adult elasticity calculations
    
    
    #===========================================================================
    # Mortality estimation
    
    # Jensen methods  
    Jensen.mat <- 1.65/(age_at_maturity)
    Jensen.k <- 1.5*k
    
    #Then et al. methods that update the old Hoenig and Pauly methods
    Then_hoenig <- 4.899*tmax^-0.916
    Then_pauly <- 4.118*k^0.73*Linf^-0.33
    
    #Peterson & Wroblewski method, this method may give time varying m ###
    Peterson <- 1.92*(Wt*1000)^-0.25

    
    # List all M estimates in a table
    M.table <-  cbind(Jensen.k, Jensen.mat, Then_hoenig, Then_pauly, Peterson)
    
    #"sample" above randomly gets 1 value of mortality from those listed in the table
    #"min" for maximum survivorship from all methods - note, taking a random sample from one row
    #DK note: There is no sample in this as far as I can tell, so it is just taking the first input from this table
    Mortality <- M.table[,min(x=ncol(M.table),size=1)]
    
    # Estimate Survivorship from M and F estimates. If an AAFC is set then F is only applied to ages
    # older than the AAFC - note that this will be called in the simulation function
    survivorship <- NULL
    # browser()
    for(j in 1:max.age)
    {
      # So these gets F as a vector, assumes F is same for all ages that can be harvested, 
      # DK note: But those ages appear only to be ages that haven't matured yet, which must be a very Porbeagle thing that we'll need to change
      # I think this can come outside this loop as it won't change with each j as far as I can tell.
      f.vec<-c(rep(0,AAFC+1),rep(F,age_at_maturity-AAFC),rep(0,max.age-age_at_maturity))  
      # Below AAFC survivorship is just natural mortality
      if(Age[j]<=AAFC-1) survivorship[j] <- exp(-(Mortality[j]))
      # Above AAFC it includes the f term.  DK note: I don't think this if-else is really needed because in the f.vec you've set F = 0 below AAFC.
      else{survivorship[j] <- exp(-(Mortality[j]+f.vec[j]))}
    } # end for(j in 1:max.age)
    
    # Add a zero at the end of the survivorship vector to identify that no individuals survive
    # the final age class
    sx <- append(survivorship,0,after=length(survivorship))
    mx <- Fecundity_ogive * n.pups.year  ##fecundity at age
    a <- mx[-1]# line up vector so sx and mx are offset by 1 so that eventually fx = mx*sx+1...
    MX <- append(a,0,after=length(a)) # add zero to the end so both vectors are the same length again.
    fx <- MX*sx # Create top row of matrix
    Sm <- matrix(rep(0),nrow=length(survivorship),ncol=length(survivorship)) # create a matrix entirely of zeros
    diag(Sm) <- survivorship #insert survival along the diagonal
    col <- rep(0,length(survivorship)) # create extra row to designate 0 values for final age class
    Sm <- cbind(Sm,col,deparse.level = 0) # survival matrix with 0's for final age class
    matrix <- rbind(fx,Sm,deparse.level = 0)# join fecundity vector to survival matrix without row names (deparse.level=0)
    
    ############ Matrix projection ##########
    # calculate lambda, convert to r for other calculations and list Lambda in a table for future analysis
    lambda <- eigen.analysis(matrix)$lambda
    lambda.table[i] <-  lambda 
    r <- log(lambda) # Convert lambda to r value
    r.table[i] <- r
    
    #convert sx to lx values for R0 and G calculations
    lx <- NULL
    lx[1] <- 1
    for(j in 2:max.age) 
    {
      if(j == 2) lx[j] <- survivorship[1]
      if(j > 2) lx[j] <- survivorship[j]^(j-1)
    } # end for(j in 2:max.age) 
    # Make the last age a 0
    lx[max.age+1] <- 0
    
    # Net reproductive rate, fecundity times survivorship
    R0 <- sum(mx*lx)
    R0.table[i] <- R0 
    
    #Generation time
    G <- sum(Age*(exp(-r*Age)*(mx*lx)))
    G.table[i] <- G 
    
    #Maximum lifetime reproductive rate, steepness, position of inflection point, and SPRmer
    #First year survivorship (s0)
    s0 <- lx[2]
    s0.table[i] <- s0 # Table of the survivorship values for age 1
    # Alpha hat is the maximum lifetime reproductive rate, e.g the maximum number of female spawners that can be produced by a female during her life
    alpha.hat <- s0*R0 # So this is reproductive rate times survivorship of age 1. 
    # DK Note, so I think the alpha.hat estimate relying on assuming a constant M across age classes, is it ok to have M vary (time or by age) here?  There is a SPR method in fishMethods package that can
    # handle variable M, wondering if we should explore that?  I probably should read Brooks 2010 - Analytical reference points for age-structured models: application to data-poor fisheries
    alpha.hat.table[i] <- alpha.hat # Put it into a table
    steepness <- alpha.hat/(4+alpha.hat) # DK note: Where does the 4 come from?
    steepness.table[i] <- steepness
    # So this is some sort of reproductive metric accounting for generation time
    R <- 0.633-(0.187*log(r*G)) # DK note: Where does this come from
    R.table[i] <- R
    SPRmer <- 1/sqrt(alpha.hat) # Spawning potential ratio at maximum excess recruitment, Reference point from Brooks et al. 2010 paper.
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
