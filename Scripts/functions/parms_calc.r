# Calculate your parameters

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