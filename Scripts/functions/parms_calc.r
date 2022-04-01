# This function returns some of the key life history parameters for your stock. This does lots of the heavy lifting and will need to be generalized to work for any old species.

#Arguments for the function.  The original defaults were set to to Porbeagle example default.

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
#w.a                  The Weight-length relationship intercept. Defaults to the value for Porbeagle (1.482E-5)
#w.b                  The Weight-length relationship exponent. Defaults to the value for Porbeagle (2.9641)


# Reproduction parameters. We'll need to tweak language here for fish.
# n.offspring         Average number of offspring which is transformed to number of female progeny per female per year.  Defaults to the value for Porbeagle (4)
# gest.period         The gestation period in years. Defaults to the value for Porbeagle (9/12)


# max.age.lb         The lowest likely 'maximum age'. Defaults to the value for Porbeagle (25)
# repro.cycle        The number of times you reproduce in a year.  Defaults to 1 now
# vb.cor.mat         The correlation matrix for the von B parameters to ensure that we don't get weird parameter combinations.  Default is from PorBeagle, will have to sort out if this is fine for all.
# mat.cor.mat        The correlation matrix for the maturity parameters to ensure that we don't get weird parameter combinations.  Default is from PorBeagle, will have to sort out if this is fine for all.

parms.calc<-function(linf.mn = 309.8, k.mn = 0.061, t0.mn = -5.90, linf.se = 13.36, k.se = 0.007, t0.se = 0.474,
                     mat.int.mn = -10.2899, mat.slop.mn = 0.7299,  mat.int.se = 1.6793,mat.slope.se = 0.1179,
                     w.a = 1.482E-5, w.b = 2.9641,
                     max.age.lb = 25,n.offspring = rnorm(1,4,0),gest.period = 9/12, repro.cycle =1,
                     vb.cor.mat = 'default',mat.cor.mat = 'default')
{

  ## use to get max age, age at maturity, mx and M for each iteration of the simulation
  ############################ Parameters that we'll need to turn into function arguments  ############################



  ############################ Parameters that we'll need to turn into function arguments  ############################

  # Von B calculations.  We have error covariance matrix this is used to get draws of possible reasonable parameter estimates
  # of the Von B parameters for the simulations.  Here we need to sort out if this correlation matrix is sufficient in all situations or if we need to update it
  # for different life histories.  Best idea if needs changed would be to have a few 'general' ones rather than having to make the user specify this matrix every time.
  if(vb.cor.mat == 'default')
  {
    e.cor <- matrix(c(1, -0.99, -0.88,
                    -0.99, 1, 0.94,
                    -0.88, 0.94, 1),
                    ncol = 3, dimnames = list(c("Linf", "k", "t0")))
    colnames(e.cor) <- c("Linf", "k", "t0")
  }
  # covariance matrix that uses the correlation matrix of VBGF parameters to simulate possible von B parameter combos.
  e.cov <- e.cor*as.matrix(c(linf.se, k.se, t0.se))%*%t(as.matrix(c(linf.se, k.se, t0.se))) # The covariance matrix
  VB_pars <- mvrnorm(1, mu = c(linf.mn,k.mn,t0.mn), Sigma = e.cov)
  Linf <- as.numeric(VB_pars["Linf"])
  k <- as.numeric(VB_pars["k"])
  t0 <- as.numeric(VB_pars["t0"])


  # Simulate possible Maturity Ogives  slopes and intercept combos
  # Correlation matrix for population maturity ogive. Should check to make sure this correlation structure is generalizable.
  if(mat.cor.mat == 'default')
  {
    mat.cor <- matrix(c(1, -0.99,
                      -0.99,  1),
                      ncol = 2, dimnames = list(c("Intercept", "slope")))
    colnames(mat.cor) <- c("Intercept", "slope")
  }
  #Maturity ogive
  Mat.cov <- mat.cor*as.matrix(c(mat.int.se, mat.slope.se))%*%t(as.matrix(c(mat.int.se, mat.slope.se)))
  Mat_pars <- mvrnorm(1, mu = c(mat.int.mn,mat.slop.mn), Sigma = Mat.cov)
  Mat_intercept <- as.numeric(Mat_pars["Intercept"])
  Mat_slope <-as.numeric(Mat_pars["slope"])

  # Here we are estimating the fecundity of the population
  #repro.cycle <- 2 # reproductive periodicity in years
  n.fem.pups <- n.offspring/2# average litter size divided by two to account for just females.
  n.pups.year <- n.fem.pups/repro.cycle #Number of female pups per year

  # Getting the bounds on the maximum age
  tmax <- 7*log(2)/k # Where does this relationship come from?
  max.age.upper.bound <- tmax
  # Grab one of the possible max ages for the simulation
  max.age <- round(runif(1,max.age.lb,max.age.upper.bound))
  Age <- 0:max.age # age vector for age specific values

  # Length-at-age from multivariate VBGF dist
  Lt <- Linf*(1-exp(-k*(Age-t0)))
  # Calculate the weight using the length-at-age
  Wt <- w.a*Lt^w.b

  # Maturity
  Maturity_ogive <- round(1/(1+exp(-(Age*Mat_slope+Mat_intercept))),2)
  age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive


  # Fecundity
  Fecundity_ogive <- round(1/(1+exp(-((Age-gest.period)*Mat_slope+Mat_intercept))),2)
  age_at_first_repro <- round(gest.period+age_at_maturity) # for adult elasticity calculations
  mx <- Fecundity_ogive * n.pups.year


  #===========================================================================
  # Natural mortality estimation, note that there are 5 methods used here to estimate natural mortality and for each simulation draw only one of these
  # estimates are choose (at random)

  ### Currently these methods give a fixed natural mortality ###
  # Jensen methods
  Jensen.mat <- 1.65/(age_at_maturity)
  Jensen.k <- 1.5*k
  #Then et al. methods that update the old Hoenig and Pauly methods
  Then_hoenig <- 4.899*tmax^-0.916
  Then_pauly <- 4.118*k^0.73*Linf^-0.33
  ###

  ### This method may give time varying m ###
  #Peterson & Wroblewski method
  Peterson <- 1.92*(Wt*1000)^-0.25
  ###

  # List all M estimates in a table
  M.table <-  cbind(Jensen.k, Jensen.mat, Then_hoenig, Then_pauly, Peterson)
  #"sample" above randomly gets 1 value of mortality from those listed in the table
  #"min" for maximum survivorship from all methods - note, taking a random sample from one row
  ## -- DK Note, I'm pretty sure this just takes the first column, so always is picking Jensen.k. From above comment I don't think this is expected behavior of this -- ##
  Mortality <- M.table[,min(x=ncol(M.table),size=1)]

  return(list(r.M=Mortality,mx=mx,max.age=max.age,age.mat=age_at_maturity))

}