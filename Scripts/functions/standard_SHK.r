# Big important function



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
