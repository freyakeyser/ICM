# Function to calculate Lotka R, DK revised this with several options of how we want to calculate

# Function Arguments
# yrs:       Do you want to run this for multiple years, if so put the year range here, default = 1 which will just do the calculation once
# age.mat:     Age at maturity if one value it is the age at 50% maturity, if a vector/matrix it is the age at maturities.
# nat.mort:    Natural mortality, this wants the instantaneous natural mortality, not the proportional. There are several options here....
#              You can enter a vector of instantaneous natural moralities, just put a single number of all ages, or put in a character
#              string that will tell the function to calculate the natural mortality based on life history data.
# ages:        What are the ages you are using 
# fecund:      How are you estimating 
# wt.at.age:   The weight of individuals for each age
#mat.ogive.K   Only used if you provide just the age at 50% maturity.  Default is -50 which is effectively a knife edge maturation, as you approach 0 this becomes more smooth.
# sim          Are we doing a retrospective analysis or a forward projection.  default = 'retro', set sim = 'project' for forward projections
# proj.sim     How do you want to get your samples if you are doing the forward projection. You can 'sample' from the past directly 
#              or use past data and add uncertainty through a lognormal distribution based on the historic data using proj.sim = 'dist'

lotka.r<-function(yrs = 1,age.mat=4,nat.mort = NULL,ages = ages,wt.at.age = NULL,fecund = NULL,
                  mat.ogive.K = -50, sim = 'project',proj.sim = 'dist',
                  L.inf = NULL,K = NULL,t0 = NULL, a.len.wgt = NULL, b.len.wgt = NULL, a.fec.len = NULL, b.fec.len = NULL,
                  sd.mat = 0,sd.nm = 0,sd.wt = 0,sd.fecund = 0
                  # The above are needed if 1) no weight at age data (L.inf,K, t0, a.len.wgt, b.len.wgt) 
                  #                         2) Using the "Then_pauly" nat mort (L.inf and K only) 
                  #                         3) using the "Jensen K" nat mort(K only needed)
                  #                         4) if using the 'eggs' method to estimate fecundity. (L.inf,K, t0, a.fec.len,b.fec.len)
) 
{
  if(is.character(nat.mort)[1]) nm.c <- nat.mort else nm.c <- "Using data"
  if(is.character(fecund)[1]) fecund.c <- fecund else fecund.c <- "Using data"
  # Some warning messages if you didn't include something you needed...
  #if(nat.mort == "Jensen Maturity")  stop("maybe)
  if(nm.c == "Jensen K" & is.null(K))   stop("You need to specify K from von B to use Jensen's K method")
  #if(nat.mort == "Then_hoenig")     stop("maybe)
  if(nm.c == "Then_pauly")
  {
    if(is.null(L.inf) | is.null(K)) stop("You need to specify K and L infinity from von B to use the Then-Pauly method")
  }
  
  if(is.null(wt.at.age))
  {
    if(is.null(L.inf) | is.null(K) | is.null(t0)) stop("You need to specify K, L infinity, and t0 if weight at age vector isn't specified")
  }
  
  if(fecund.c == 'eggs')
  {
    if(is.null(L.inf) | is.null(K) | is.null(t0) | is.null(a.fec.len) | is.null(b.fec.len)) stop("You need to specify K, L infinity, t0, a.fec.len, and b.fec.len if you are using the 'eggs' method")
    if(min(ages) > 0) stop(paste("If using the eggs fecundity option your ages must range from 0 to", max(ages),"you also need a natural mortality estimate for age 0"))
  }
  
  # The Euler-Lotka Function to optimize on to find the intrinsic rate of growth for your population
  eulerlotka <- function(r) (sum(lx * mx * exp(-r * ages)) - 1)^2
  # Set up a dataframe
  res <- data.frame(yrs = NA, r = NA)
  if(length(yrs) > 1) n.yrs <- length(yrs) else n.yrs <- 2 # Needs to be 2 if just doing it for one year to work how I have loop set up
  
  
  if(sim == 'retro')
  {
  
  # Get the maximum age
  max.age <- max(ages)
  # If you are using the fecund = ages, then we need to start from age 0
  #if(fecund == 'eggs') ages <- 0:max.age
  
  # Get the maturity ogive and age at 50% maturity (if needed) if you just gave me the age at 50% maturity
  if(length(age.mat) == 1) 
  { 
    # Add some uncertainty to the age at 50% maturity then fit the maturity ogive
    age.at.50.mat <- rlnorm(1,log(age.at.50.mat),sd.mat)
    if(mat.K <= -50) print(paste("Note that you are assuming knife edge selectivity here with 50% maturity at age", age.mat))
    if(nm.c == "Jensen Maturity") age.at.50.mat <- age.mat
    # This gets us a hypothetical maturity ogive, the default (mat.K = -50) is essentially knife edge maturity
    mat.ogive <- 1/(1+exp(-(mat.K)*(age.at.50.mat-ages)))
    
  } # end if(length(age.mat) == 1) 
  
  # If you have given me the maturity ogive than rename it
  if(length(age.mat) > 1) 
  {
    if(is.vector(age.mat)) mat.ogive <- rlnorm(length(age.mat),age.mat,sd.nm)
    if(!is.vector(age.mat)) 
    {
      mat.ogive <- as.data.frame(age.mat) # Get it set up, then overwrite it with the uncertainty
      for(i in 1:nrow(age.mat)) mat.ogive[i,] <- rlnorm(length(age.mat[i,]),as.numeric(log(age.mat[i,])),sd.nm)
    }
    
    
    # If we are doing Jensen Natural mortality method we need age at 50% maturity, which we can estimate from the maturity ogive using a gam
    if(nm.c == "Jensen Maturity") 
    {
      library(mgcv)
      mod <- gam(mat.ogive~s(ages)) # Fit a simple little GAM to estimate age at 50% maturity
      pred.dat <- data.frame(ages = seq(0,max(ages),by=0.01)) # Dataframe to predict on
      pred.dat$a.50 <- predict(mod,pred.dat) # The prediction
      age.at.50.mat <- pred.dat$ages[which.min(abs(0.5-pred.dat$a.50))] # Estimated age @ 50% matuity
    } # end if(nat.mort == "Jensen Maturity")
  } # end if if(length(age.mat) > 1) 
  
  # If we don't have weight at age data we can get that from life history data...
  if(!is.null(wt.at.age)) 
  {  
    W.age <- as.data.frame(wt.at.age)
    # Add in uncertainty
    if(is.vector(wt.at.age)) W.age <- rlnorm(length(wt.at.age),wt.at.age,sd.nm)
    if(!is.vector(wt.at.age)) for(i in 1:nrow(wt.at.age)) W.age[i,] <- rlnorm(length(wt.at.age[i,]),as.numeric(log(wt.at.age[i,])),sd.wt)
  }
  
  if(is.null(wt.at.age) | fecund.c == 'eggs')
  {
    L.age <- L.inf*(1-exp(-K*(ages-t0))) # Need this for both eggs case and missing weight case.
    # We can get weight at age now... this is in grams
    if(is.null(wt.at.age)) 
    {
      W.age <- a.len.wgt*L.age^b.len.wgt
      # Add in uncertainty....
      W.age <- rlnorm(length(W.age),as.numeric(log(W.age)),sd.wt)
    }
    
  } # if(is.null(wt.at.age) | fecund =='eggs')
  
  # Add any uncertainty to the data we want to add.
  if(!is.vector(nat.mort)) 
  {
    nat.mort <- as.data.frame(nat.mort)
    for(i in 1:nrow(nat.mort)) nat.mort[i,] <- rlnorm(length(nat.mort[i,]),as.numeric(log(nat.mort[i,])),sd.nm)
  }
  
  
  # So if we don't have a matrix of natural mortality data...
  if(length(nat.mort) == 1)  
  {
    # If we have one number, turn that into a vector
    if(is.numeric(nat.mort)) nat.mort <- rep(nat.mort,max.age)
    # if we have a character string then we pick our option of interest...
    if(nm.c == "Jensen Maturity") nat.mort <- rep(1.65/(age.at.50.mat),max.age)
    if(nm.c == "Jensen K")        nat.mort <- rep(1.5*K,max.age)
    if(nm.c == "Then_hoenig")     nat.mort <- rep(4.899*max.age^-0.916,max.age)
    if(nm.c == "Then_pauly")      nat.mort <- rep(4.118*K^0.73*L.inf^-0.33,max.age)
    if(nm.c == "Peterson")        nat.mort <- 1.92*(W.age)^-0.25
    
    if(is.vector(nat.mort)) nat.mort <- rlnorm(length(nat.mort),nat.mort,sd.nm)
  }
  
  
  
  # Now get the mx vector, 
  if(is.character(fecund))  
  {
    # if you said it's based on egg production using fecundity-length relationship
    # This seems to go a bit high...
    if(fecund.c == 'eggs') 
    {
      mx.t <- a.fec.len*L.age^b.fec.len
      mx.t <- rlnorm(length(mx.t),as.numeric(log(mx.t)),sd.fecund)
    }
    
    
    # Alternatively we can get it from SPR formulation
    if(fecund.c == "SPR")
    {
      spr <- sbpr(age=ages,
                  ssbwgt=W.age,
                  partial=rep(1,length(ages)),
                  pmat=mat.ogive,
                  M=nat.mort, # This can just be one value if m assumed fixed.
                  plus=T,oldest =max(ages),pF=0,pM=0.5,incrF = 0.1,MSP=0)	
      #So each individual in an age class will contribute this many recruits...
      # Which is what mx should be, I think??
      mx.t <- nat.mort*W.age*mat.ogive/unlist(spr$Reference_Point[2])/2 # Divided by 2 because Euler below in terms of number of females produced
      # note that give all the inputs already have variation attributed to them, there is no need in the SPR to add in more variability.
    } # end if(fecund == "SPR")
    
  } # end if(is.character(fecund))  
 
  # Or if we directly had an estimate of fecundity (number of recruits per individual)
  if(!is.character(fecund)) 
  { 

    mx.t <- fecund
    spr <- NULL
    if(is.null(nrow(fecund))) mx.t <- rlnorm(length(mx.t),mx.t,sd.fecund)
    if(nrow(mx.t) > 1) for(m in 1:nrow(mx.t)) mx.t[m,] <- rlnorm(length(mx.t[m,]),as.numeric(log(mx.t[m,])),sd.fecund)
  }
  # remove columns associated with NAs in mx, and remove NAs from mx
  if(any(is.na(colSums(mx)))) {
    ages <- ages[-which(is.na(colSums(mx)))]
    nat.mort <- nat.mort[-which(is.na(colSums(mx)))]
    mx <- mx[-which(is.na(colSums(mx)))]
    mx.t <- mx.t[-which(is.na(colSums(mx.t)))]
  }
  } # end the section to get mx and lx for retrospective analyses
  ################################### END RETRO SECTION ################################### END RETRO SECTION################################### END RETRO SECTION
  
 
  
  
  #################################### PROJECTION SECTION#################################### PROJECTION SECTION#################################### PROJECTION SECTION
  # Now build the projection script for lotka.r
  if(sim == 'project')
  {
    # Things we aren't using in projections (yet at least)
    spr = NA
    mat.ogive = NA
    W.age = NA
    # going to start simple here and just have it set up to work for the stocks we have the necessary data for
    mx.t <- matrix(nrow=n.yrs,ncol= ncol(fecund))
    nm.t <- nat.mort
    # So grab a sample from the past as our basis
    f.sample <- round(runif(n.yrs,min=1,max = nrow(fecund)))
    m.sample <- round(runif(n.yrs,min=1,max = nrow(nat.mort)))
    # then you get your sample for the projections
    mx.t <- fecund[f.sample,]
    nat.mort <- nat.mort[m.sample,]
    
    # Now if you want to add uncertainty to this run through a loop and add uncertainty.
    if(proj.sim == 'dist')
    {
      for(y in 1:n.yrs)
      {
        # Now if you want to use the fecund/nat.mort data and add noise to that...
        # But need to be a bit careful here, I don't want them all jumbled around across the ages
        # as they all would tend to respond similarly (i.e. a year with high M probably won't just be high m for one ages)
        # So going to take mean of whole vector and adjust them all accordingly...
        log.mean.nm <- log(mean(as.numeric(nat.mort[y,])))
        log.mean.mx <- log(mean(as.numeric(mx.t[y,])))
        nm.deva <- rlnorm(1,log.mean.nm,sd.nm)/exp(log.mean.nm )
        mx.deva <- rlnorm(1,log.mean.mx,sd.fecund)/exp(log.mean.mx )
        mx.t[y,] <- mx.t[y,]*mx.deva
        nat.mort[y,] <- nat.mort[y,]*nm.deva  
      } # end for(y in 1:n.yrs) loop
     } # end if(proj.sim == 'dist')

  }
  
  ################################### END PROJECTION SECTION ################################### END PROJECTION SECTION################################### END PROJECTION SECTION
  
  
  
  #################################### Calculate Lotka SECTION #################################### Calculate Lotka SECTION #################################### Calculate Lotka SECTION
  mx.tmp <- NULL
  lx.tmp <- NULL
  
  # Now run this through all the yrs 
  for(j in 1:(n.yrs))
  { 
    # Get the r estimate for each year if we are running through a bunch of years
    # Just use the first year if we don't have same length nat.mort as years
    if(!is.null(nrow(nat.mort))) lx <- 1-exp(-nat.mort[j,]) else lx <- 1-exp(-nat.mort)
    # Convert to survivorship vector
    si <- 1-lx 
    # Set the first age class to be 1
    lx<-1
    # And get cumulative survivorship for the stock
    for(s in 2:(length(si))) lx[s]<-lx[s-1]*si[s-1]
    
    if(!is.null(nrow(mx.t))) mx <- mx.t[j,] else mx <- mx.t
    # Now we are cooking!
    junk<-nlminb(start = 0.1, obj = eulerlotka)
    res[j,] <- c(yrs[j],junk$par)
    mx.tmp[[j]] <- mx
    lx.tmp[[j]] <- lx
    #browser()
  } # end for(j in 1:(n.yrs-1))
  
  mx.mat <- do.call("rbind",mx.tmp)
  lx.mat <- do.call('rbind',lx.tmp)
  
  
  return(list(res=res,mx=mx.mat,lx=lx.mat,spr = spr,mat.ogive = mat.ogive,wt.at.age = W.age,nat.mort = nat.mort))
}