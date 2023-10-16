# Here we take the ideas from the backwards projections but 'tune' them to the 'known' biomass by adjusting either the fecundity or survivorship variables iteratively until
# the biomass in the current year is withing X% of the actual biomass.

# All the data you need to make it dance...
# years:          The years you are running the backwards calculation for
# abund.ts:       The abundance time series from the assessment

# tuner:          What are we tuning, the fecundity "f", or natural mortality 'm', or both 'b'
# rems            Removals from the fishery in a given year.
# step.size:      How big is the difference (proportional) in the simulation step.  Default = 0.05 which is about a 5% change for each step

# direction       Allow for forwards and backwards simulations
# What we need for the Lotka.r function...
# fecund:      How are you estimating fecundity. If a vector it is the number of recruits produced by the average individual in the age classes. 
#              If a matrix the rows are different simulations. Also have options to set fecund = 'eggs' which we can use if we can get an estimate of egg mortality
#              fecund = 'SPR' uses the spawner per recruit metric, to get this to work we need to know how many recruits there are in a given year
#              We could also let this vary by year, but we'll need to go 3-D array or something for that
# nm:             Natural mortality, this wants the instantaneous natural mortality, not the proportional. There are several options here....
#                 You can enter a vector of instantaneous natural moralities, just put a single number of all ages, or put in a character
#                 string that will tell the function to calculate the natural mortality based on life history data. If a matrix the rows are different simulations
#                 we could also let this vary by year, but we'll need to go 3-D array or something for that
# ages:       Ages...


fast.tunes<-function(years,step.size=0.05,tuner="m",  direction = 'backwards',ages = NULL,
                   nm=NULL,fecund = NULL,N.init = NULL,abund.ts = NULL, rems)
{
  # Download the function to go from inla to sf
  funs <- c(#"https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/simple.lotka.r",
            "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/backwards_project.r",
            "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_project.r"
            
  )
  # Now run through a quick loop to load each one, just be sure that your working directory is read/write!
  for(fun in funs) 
  {
    download.file(fun,destfile = basename(fun))
    source(paste0(getwd(),"/",basename(fun)))
    file.remove(paste0(getwd(),"/",basename(fun)))
  }
  
  source("D:/Github/ICM/Scripts/functions/simple_Lotka_r.r")
  
  
  require(optimx)  || stop("Please load the 'optimx' package which you'll need for the optimations to run")
  
 
  #Initialize a bunch of objects
  n.years<-length(years)
  res <- data.frame(year = years,vpa.abund = N.init, est.abund = NA, diff= NA,per.diff = NA,r=NA,r.vpa.init = NA,
                    removals = NA,mean.fec=NA,mean.vpa.fec = NA,mean.nm = NA,mean.vpa.nm = NA)
  fecund.ts <- data.frame(fecund)
  nm.ts <- data.frame(nm)
  
  # Here I try to make a fast fitting version of the previous model...
  
  if(direction == 'backwards') y.index <- (n.years-1):1 # Starting at the second last year of data makes the indexing below way more consistent between directions
  if(direction == 'forwards')  y.index <- 1:(n.years-1) # Stop at the second last year of data
  
  N.up <- N.init

  for(y in y.index)
  {
   #if(y %in% c(27)) browser()
    print(paste0("Year = ",years[y]))
    #if(y %in% c(43,39,38)) browser()
    pop.est <- N.up
    # The first step is to fit the default lotka.r values
    # When going forwards we use everything in the 'current' year to project to next year
    # When going backwards we use everything in the 'previous' year to project from the next year
    
    nm.tmp <- nm.init <- nm[y,]
    fecund.tmp <- fecund.init <-  fecund[y,]
    
      
    # DK Note: So for our removals time series, we put the removals between y+1 and y 
    # down as year y+1.  We can change this, but that's how this is set up at the moment.
    # So when going forwards we want to use the current year index, but when going backwards
    # we want to use the previous years. This is all handled by the y.index happily.
    removals.next <- rems[y]
    #browser()
    # Going to line up the r with the 'previous' year, so the r in column 2020 is the r used to move between 2020 and 2021 (in either direction)
    lotka.est <- simple.lotka.r(nat.mort = nm.tmp,fecund=fecund.tmp,ages=ages)
    
    r.est = lotka.est$res
    
    # For extremes, get away from that bad estimate area fast as it is numerical instability territory!
    if(r.est < -0.995) 
    {
      fecund.tmp <- 1.5*fecund.tmp
      nm.tmp <- nm.tmp/1.5
    }
    # Then we run the model 1 year (either backwards or forwards) and see how well it fits.
    
    if(direction == 'backwards')
    {
        exp.res <- back.proj(option = "exponential",pop.next = pop.est,r=r.est,removals = removals.next)
        pop.step <- exp.res$Pop.current
        diff.org <- abund.ts[y] - pop.step
        per.diff.org <- 100*((pop.step - abund.ts[y])/abund.ts[y])
    } # end if(direction == 'backwards')
    
    # And forwards 
    if(direction == 'forwards')
    {
        exp.res <- for.proj(option = "exponential",pop.last = pop.est ,r=r.est,removals = removals.next,fishery.timing = 'beginning')
        if(exp.res$Pop.current < 0) exp.res$Pop.current =0 # don't let it drop below 0
        pop.step <- exp.res$Pop.current
        diff.org <- abund.ts[y+1] - pop.step
        per.diff.org <- 100*((pop.step - abund.ts[y+1])/abund.ts[y+1])
    } # end if(direction == 'forwards')
    diff <- diff.org
    per.diff <- per.diff.org
    # So now we compare that value to the 'real' value from the vpa
    #if(y %in% 3:7) browser()
    
    #Initial conditions
      count =0
      last.per.diff <- 1e6 # Make it ridiculous high to make sure we get started...
      while(abs(per.diff) < abs(last.per.diff))
      {
        #browser()
        last.per.diff <- per.diff
        # indexing...
        count = count + 1
        if(count == 1)
        {
        fecund.inc <- fecund.tmp * (1+step.size)
        fecund.dec <- fecund.tmp / (1+step.size)
        # Calculate the lotka both direction and figure out which way is making the fit better...
        lotka.inc <- simple.lotka.r(nat.mort = nm.tmp,fecund=fecund.inc,ages=ages)
        lotka.dec <- simple.lotka.r(nat.mort = nm.tmp,fecund=fecund.dec,ages=ages)
        r.inc <- lotka.inc$res
        r.dec <- lotka.dec$res
        
        # Backwards
        if(direction == 'backwards')
        {
            exp.inc <- back.proj(option = "exponential",pop.next = pop.est,r=r.inc,removals = removals.next)
            pop.est.inc <- exp.inc$Pop.current
            exp.dec <- back.proj(option = "exponential",pop.next = pop.est,r=r.dec,removals = removals.next)
            pop.est.dec <- exp.dec$Pop.current
            # Now we test the two methods...
            diff.inc <- abund.ts[y] - pop.est.inc
            per.diff.inc <- 100*((pop.est.inc - abund.ts[y])/abund.ts[y])
            diff.dec <- abund.ts[y] - pop.est.dec
            per.diff.dec <- 100*((pop.est.dec - abund.ts[y])/abund.ts[y])
        } # end if(direction == 'backwards')
        
        # And forwards
        if(direction == 'forwards')
        {
            exp.inc <- for.proj(option = "exponential",pop.last = pop.est ,r=r.inc,removals = removals.next,fishery.timing = 'beginning')
            if(exp.inc$Pop.current < 0) exp.inc$Pop.current =0 # don't let it drop below 0
            pop.est.inc <- exp.inc$Pop.current
            # And the decrease one
            exp.dec <- for.proj(option = "exponential",pop.last = pop.est ,r=r.dec,removals = removals.next,fishery.timing = 'beginning')
            if(exp.dec$Pop.current < 0) exp.dec$Pop.current =0 # don't let it drop below 0
            pop.est.dec <- exp.dec$Pop.current
            # Now we test the two methods...
            diff.inc <- abund.ts[y+1] - pop.est.inc
            per.diff.inc <- 100*((pop.est.inc - abund.ts[y+1])/abund.ts[y+1])
            diff.dec <- abund.ts[y+1] - pop.est.dec
            per.diff.dec <- 100*((pop.est.dec - abund.ts[y+1])/abund.ts[y+1])
        } # end if(direction == 'forwards')
        
       
        
        if(abs(per.diff.inc) < abs(per.diff.dec))
        {
          going <- 'up'
          fecund.tmp <- fecund.inc
          per.diff <- per.diff.inc
          #pop.tunes <- pop.est.inc
          
        } else 
          {
            going <- 'down'
            fecund.tmp <- fecund.dec
            per.diff <- per.diff.dec
            #pop.tunes <- pop.est.dec
          } # end if else
        } # end count ==1
        
        if(count > 1)
        {
          # Increase the values...
          if(going == 'up') 
          {
            fecund.tmp <- fecund.tmp*(1+step.size)
            # Once r is big, we also turn down nm!!
            if(count > 2)
            {
              #if(r.tunes > 0.5) nm.tmp <- nm.tmp/(1+step.size)
              if(abs(per.diff) > 30) 
              {
                nm.tmp <- nm.tmp/1.05
                fecund.tmp <- fecund.tmp*1.05
              }
            }
            
          }
          if(going == 'down') 
          {
            fecund.tmp <- fecund.tmp/(1+step.size)
            # Once r is small, we also turn up nm!!
            if(count > 2) 
            {
              #if(r.tunes < -0.3) nm.tmp <- nm.tmp*(1+step.size)
              if(abs(per.diff) > 30) 
              {
                nm.tmp <- nm.tmp*1.05
                fecund.tmp <- fecund.tmp/1.05
              }
            }
          }
          
          
          lotka.tunes <- simple.lotka.r(nat.mort = nm.tmp,fecund=fecund.tmp,ages=ages)
          r.tunes <- lotka.tunes$res
        
          
        if(direction == 'backwards')
        {
          exp.tunes <- back.proj(option = "exponential",pop.next = pop.est,r=r.tunes,removals = removals.next)
          pop.tunes <- exp.tunes$Pop.current
          diff.tunes <- abund.ts[y] - pop.tunes
          per.diff.tunes <- 100*((pop.tunes - abund.ts[y])/abund.ts[y])
          per.diff <- per.diff.tunes
        } # end if(direction == 'backwards')
        
        # And forwards
        if(direction == 'forwards')
        {
          exp.tunes <- for.proj(option = "exponential",pop.last = pop.est ,r=r.tunes,removals = removals.next,fishery.timing = 'beginning')
          if(exp.tunes$Pop.current < 0) exp.tunes$Pop.current =0 # don't let it drop below 0
          pop.tunes <- exp.tunes$Pop.current
          diff.tunes <- abund.ts[y+1] - pop.tunes
          per.diff.tunes <- 100*((pop.tunes - abund.ts[y+1])/abund.ts[y+1])
          per.diff <- per.diff.tunes
        } # end if(direction == 'forwards')
        
     
        } # end if count >1
        #print(paste0("Percent difference = ",round(per.diff))) 
        
        #if(y %in% c(43,39,38)) browser()
        #browser()
        #if(y %in% 3:7) browser()
        
        #if(y %in% c(27)) browser()
      } # end while(abs(per.diff) > 1)
      
    #browser()
      # If count = 1 then the original is the best fit!
      if(max(count) == 1) 
      {
          r.tunes <- r.est
          pop.tunes <- pop.step
          diff.tunes <- diff.org
          per.diff <- per.diff.org
          fecund.tmp <- fecund.init
          nm.tmp <- nm.init
      }
      
    if(direction == 'forwards')  
    {
      index <- which(years == years[y+1])
      #browser() 
      res$r[index] <- r.tunes
      #browser()
      res$r.vpa.init[index] <- r.est
      res$removals[index] <- removals.next
      res$mean.fec[index] <- mean(as.numeric(fecund.tmp))
      res$mean.vpa.fec[index-1] <- mean(as.numeric(fecund.init))
      res$mean.nm[index] <- mean(as.numeric(nm.tmp))
      res$mean.vpa.nm[index-1] <- mean(as.numeric(nm.init))
      fecund.ts[index,] <- fecund.tmp
      nm.ts[index,] <- nm.tmp
    } # end forwards
    if(direction == 'backwards')
    {
      index <- which(years == years[y])
      res$r[index+1] <- r.tunes
      res$r.vpa.init[index+1] <- r.est
      res$removals[index+1] <- removals.next
      res$mean.fec[index+1] <- mean(as.numeric(fecund.tmp))
      res$mean.vpa.fec[index] <- mean(as.numeric(fecund.init))
      res$mean.nm[index+1] <- mean(as.numeric(nm.tmp))
      res$mean.vpa.nm[index] <- mean(as.numeric(nm.init))
      fecund.ts[index+1,] <- fecund.tmp
      nm.ts[index+1,] <- nm.tmp
    }
 
    # Summary stuff
    res$vpa.abund[index] <- abund.ts[index]
    res$est.abund[index] <- pop.tunes
    res$diff[index] <- diff.tunes
    res$per.diff[index] <- per.diff
    #browser()
    
    
    # 
    # Get the next N ready...
    N.up <- pop.tunes
    
    # We then save that result and re-run everything with a new value of fecundity or natural mortality, if that value improves our
    # estimate we keep doing that.  Once it no longer improves the value we stop, save the value, and use that population value to start 
    # the next round.
    # I think what we can have the 'both' scenario do is to change fecundity in the first iteration
    # then the second iteration we change mortality, the third iteration changes fecundity... etc etc.
    # Then we need to break out of the iterative loop when the result gets worse (in terms of absolute % difference)
    # Save that result then move to the next year.... 
    #browser() 
  } # end the years loop
  
  return(list(res=res,fecund.opt = fecund.ts,nm.opt = nm.ts,fecund.org = fecund,nm.org=nm))
} # end the function  
  
  
  
  
  
  
  
  
  
  