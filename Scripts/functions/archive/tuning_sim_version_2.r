# Here we take the ideas from the backwards projections but 'tune' them to the 'known' biomass by adjusting either the fecundity or survivorship variables iteratively until
# the biomass in the current year is withing X% of the actual biomass.

# All the data you need to make it dance...
# years:          The years you are running the backwards calculation for
# abund.ts:       The abundance time series from the assessment
# tune.par:       How close to the assessment biomass do we want to get.  default is 0.01 (1%)
# tuner:          What are we tuning, the fecundity "f", or natural mortality 'm', or both 'b'
# rems            Removals from the fishery in a given year.
# n.steps         The number of different m or f options to explore.  Make sure this is an even number!

# mat.age         Age at maturity if one value it is the age at 50% maturity, if a vector it is the age at maturities for each age(or age/year), if a matrix we want this to be unique for each simulation
#                 we could also let this vary by year, but we'll need to go 3-D array or something for that# ages:           What are the ages you are using this is used to calculate max age, which may not be ideal with plus group stocks.

# nm:             Natural mortality, this wants the instantaneous natural mortality, not the proportional. There are several options here....
#                 You can enter a vector of instantaneous natural moralities, just put a single number of all ages, or put in a character
#                 string that will tell the function to calculate the natural mortality based on life history data. If a matrix the rows are different simulations
#                 we could also let this vary by year, but we'll need to go 3-D array or something for that

# sel             The selectivity of the stock.  Currently this is set up as a single number and it is tweaked in the code, we'll need to make this more complex.


# N               Population size in a given year

# u:              Exploitation rate (annual not instantaneous), currently set up to be 1 value.

# pop.model:      What method you going to use to get the population growth modeled for the backwards model.  You can use exponential model, logistic model, 
                  #or a dec.rate (decline rate) model.

# What we need for the Lotka.r function...
# fecund:      How are you estimating fecundity. If a vector it is the number of recruits produced by the average individual in the age classes. 
#              If a matrix the rows are different simulations. Also have options to set fecund = 'eggs' which we can use if we can get an estimate of egg mortality
#              fecund = 'SPR' uses the spawner per recruit metric, to get this to work we need to know how many recruits there are in a given year
#              We could also let this vary by year, but we'll need to go 3-D array or something for that
# mat.ogive.K 
# L.inf 
# K 
# t0 
# a.len.wgt 
# b.len.wgt 
# a.fec.len 

# w.age           The weight of individuals for each age, if a matrix the rows are different simulations
#                 we could also let this vary by year, but we'll need to go 3-D array or something for that
# K
# dec.rate  


tune.sim<-function(years,tune.par=0.01,tuner="fec", n.steps = 30,
                   mat.age = NULL,nm=NULL,w.age = NULL,ages =NULL,fecund = NULL,N.end = NULL,abund.ts = NULL,
                   sel,rems,N,u,pop.model = "exp",
                   dec.rate = NULL,
                   L.inf = NULL,K = NULL,t0 = NULL, a.len.wgt = NULL, b.len.wgt = NULL, a.fec.len = NULL, b.fec.len = NULL,
                   sd.mat = 0,sd.nm = 0,sd.wt = 0,sd.fecund = 0,CC =1e20)
{
  # Download the function to go from inla to sf
  funs <- c("https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/Lotka_r.r",
            "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/backwards_project.r"
  )
  # Now run through a quick loop to load each one, just be sure that your working directory is read/write!
  for(fun in funs) 
  {
    download.file(fun,destfile = basename(fun))
    source(paste0(getwd(),"/",basename(fun)))
    file.remove(paste0(getwd(),"/",basename(fun)))
  }
  
  require(optimx)  || stop("Please load the 'optimx' package which you'll need for the optimations to run")
  
  #Initialize a bunch of objects
  n.years<-length(years)
  #Pop<-data.frame(abund = rep(NA,n.years*n.sims),sim = rep(1:n.sims,n.years),years = sort(rep(years,n.sims)))    #matrix(data=NA, nrow=<<see below>>, ncol=<<see below>>,byrow=F, dimnames=NULL) 
  Pop.vec<-rep(NA,n.years)
  rem<-rep(NA,n.years)
  r.vec<-NULL
  #B1.vec<-rep(0,n.sims) ### backwards projections 
  #rems<-matrix(rep(0,n.years*n.sims),n.sims,n.years)

  
  
  # So I think what I want to create here is a big list of natural mortality and fecundity 'scenarios' that we select from until we find one that 'works'
  
  if(tuner == 'nm')
  {
    #nm.lst <- NULL
    #fecund.lst <- NULL
    # So here what I'm going to do is make a list of fecundity scenarios which we'll select from to run the model
    # Then pick the scenario that minimizes the differences from the VPA.
    # let's just make X steps around the mean either way, changing m by 5% every time....
    nm.vary <- NULL
    count <- 0
    # THis is all temp for testing...
    
    for(s in (n.steps/2):1) 
    {
      count <- count + 1
      if(s == n.steps/2) nm.vary[[s]] <- nm
      if(s < n.steps/2) nm.vary[[s]] <- nm.vary[[s+1]]/(1+tune.par)
    }
    #for(s in ((n.steps/2)+1):n.steps) nm.vary[[s]] <- nm.vary[[s-1]]*(1+ (tune.par*(s-n.steps/2)))
    
    for(s in ((n.steps/2)+1):n.steps) nm.vary[[s]] <- nm.vary[[s-1]]*(1+tune.par)
    #nm.vary <- do.call('rbind',nm.vary)
    # Then bundle this up into a new nm.lst object we'll use later
    #nm.lst[[as.character(years[y])]] <- nm.vary
  } # end   if(tuner == 'nm')
  #browser()
  if(tuner == 'fec')
  {
    #nm.lst <- NULL
    #fecund.lst <- NULL
    # So here what I'm going to do is make a list of fecundity scenarios which we'll select from to run the model
    # Then pick the scenario that minimizes the differences from the VPA.
      # let's just make X steps around the mean either way, changing m by 5% every time....
      fec.vary <- NULL
      count <- 0
      # THis is all temp for testing...

      for(s in (n.steps/2):1) 
      {
       count <- count + 1
       if(s == n.steps/2) fec.vary[[s]] <- fecund
       if(s < n.steps/2) fec.vary[[s]] <- fec.vary[[s+1]]/(1+tune.par)
      }
      for(s in ((n.steps/2)+1):n.steps) fec.vary[[s]] <- fec.vary[[s-1]]*(1+tune.par)
      #fec.vary <- do.call('rbind',fec.vary)
      # Then bundle this up into a new nm.lst object we'll use later
      #nm.lst[[as.character(years[y])]] <- nm.vary
    } # end   if(tuner == 'fec')

  #Calculate r for your examples... going to be a lot of lotka.r's all of a sudden isn't there...
  r.vec <- NULL
  for(ss in 1:n.steps)
  { 
    tmp.y <- NULL
    for(y in 1:n.years)
    {
      
      mat.tmp <- mat.age[y,]
      waa.tmp <- w.age[y,]
      if(tuner == 'fec') { fecund.tmp <- fec.vary[[ss]][y,] } else {fecund.tmp <- fecund[y,]}
      if(tuner == 'nm') { nm.tmp <- nm.vary[[ss]][y,] } else {nm.tmp <- nm[y,]}
      #tmp.s <- NULL
      #browser()
      junk<-lotka.r(yrs = years[y],age.mat = mat.tmp,nat.mort = nm.tmp,ages=ages,wt.at.age=waa.tmp,fecund=fecund.tmp,
                    L.inf = L.inf,K = K,t0 = t0, 
                    a.len.wgt = a.len.wgt, b.len.wgt = b.len.wgt, 
                    a.fec.len = a.fec.len, b.fec.len = b.fec.len,
                    sd.mat = 0,sd.nm = 0,sd.wt = 0,sd.fecund = 0)
      tmp.y[[as.character(years[y])]] <- data.frame(year = years[y],r = junk$res[1,2],mn.fec = mean(as.numeric(fecund.tmp),na.rm=T),s=ss)      
    } # end for(s in 1:30)

    r.vec[[ss]] <- do.call('rbind',tmp.y)
  }

  #unwrap your r vector
  #r.vec <- do.call('rbind',tmp.s)
  #browser()
   #temp.r.vec<-r.vec[r.vec>0 & r.vec<r.cutoff]
  #r.vec<-temp.r.vec[1:n.sims]
  ############################################################
  

  #r.tmp <- r.vec[r.vec$n.sims == i,]
  #Now run your model backwards with the r from the Lotka function and
  # if you use the logistic growth model the K estimated for the population.
  res <- NULL
  ss.diff <- NULL
  for(ss in 1:n.steps)
  {
    r.tmp <- as.data.frame(r.vec[ss])
    pop.next <- c(rep(NA,n.years-1),N.end) # set/reset the final year population abundance.
    for(y in n.years:2)
    {
      # DK Note: So for our removals time series, we put the removals between t+1 and t 
      # down as year t+1.  We can change this, but that's how this is set up at the moment.
      removals.next <- rems[y-1]
      #browser()
      # So what I want to build here is a simulation that tests the r scenarios sequentially, then at the 
      # end I'll pull out the one that fits the best of the simulation runs.
  
        # Need to think here
        r.up <- r.tmp %>% dplyr::filter(year == years[y]) %>% dplyr::pull(r) # Grab the correct value of r
        # So this one is the exponential, but when above whatever bound you have set the population growth rate averages 0 with a little uncertainty
        if(pop.model == 'bounded_exp') 
        {
          #browser()
          if(pop.next[y] > CC & r.up > 0) 
          {
            r.up <- rlnorm(1,0,0.02)-1
          }
          exp.res <- back.proj(option = "exponential",pop.last = pop.next[y] ,r=r.up,removals = removals.next,direction,fishery.timing = 'beginning')
          #browser()
          if(exp.res$Pop.current < 0) exp.res$Pop.current =0 # don't let it drop below 0
          pop.next[y-1]  <- exp.res$Pop.current
        } 
        # The exponential model
        if(pop.model == 'exponential') 
        {
          exp.res <- back.proj(option = "exponential",pop.next = pop.next[y],r=r.up,removals = removals.next)
          if(exp.res$Pop.current < 0) exp.res$Pop.current =0 # don't let it drop below 0
          pop.next[y-1] <- exp.res$Pop.current
        } # end exponential
        # If you are running the logistic model
        if(pop.model == 'logistic')
        {
          #browser()
          log.res <- for.proj(option = "logistic",pop.last = pop.next[y],K=CC,r=r.up,removals = removals.next)
          if(log.res$Pop.current < 0) log.res$Pop.current =0 # don't let it drop below 0
          pop.next[y-1]  <- log.res$Pop.current
        }
        if(pop.model == 'dec.rate')
        {
          # This needs way more thought than it's been given by DK!
          #removals[i,]<-(-1+exp(r.vec[i])+dec.rate[i]) ## redone with Jamie.
          pop.next[y] <- pop.next[y] + pop.next[y]*(-1+exp(r.up)+dec.rate[i])
        } # end dec.rate
        
      } # end the years loop.
  res[[ss]] <- data.frame(abund = pop.next,mx = r.tmp,years = years,rem = rems,vpa.abund = abund.ts,scenario= ss)
  # The differences between the ts.
  #browser()  
  ss.diff[[ss]] <- data.frame(sq.diff = (res[[ss]]$abund - res[[ss]]$vpa.abund)^2,
                              diff = res[[ss]]$abund - res[[ss]]$vpa.abund,
                              prop.diff = ((res[[ss]]$abund - res[[ss]]$vpa.abund)/res[[ss]]$vpa.abund),
                              prop.dist = ((res[[ss]]$abund - res[[ss]]$vpa.abund)/res[[ss]]$vpa.abund),
                              scenario = ss)
  # The negative prop distances need replaced with the crazy legs metric....
  ss.diff[[ss]]$prop.dist[ss.diff[[ss]]$prop.dist < 0] <- (1/(1+ss.diff[[ss]]$prop.dist[ss.diff[[ss]]$prop.dist < 0]))-1
  # When population goes to 0 make it a huge penalty...
  ss.diff[[ss]]$prop.dist[is.infinite(ss.diff[[ss]]$prop.dist)] <- 1000
  }# End the SS 
        
   diffs <- do.call('rbind',ss.diff)

   # So can I invent a symmetric metric around 0 for a proportion? So typing out loud, a percentage increase of 100% 
   # the same distance as a decline of 50%, so a 1 = -0.5, so we just want the inverse -1?
   # So for the negatives does this transform do what we want? (1/(1+prop))-1
   squalid <- diffs %>% dplyr::group_by(scenario) %>% dplyr::summarise(ssd = sum(sq.diff),
                                                                       per.diff = 100*mean(prop.diff),
                                                                       prop.dist = sum(prop.dist))
   
   # Now select the scenario that works the best
   #browser() 
   pickem <- which(squalid$prop.dist == min(squalid$prop.dist))
   res.final <- res[[pickem]]
   ss.diff.final <- ss.diff[[pickem]]
   if(tuner == 'fec') tuning.final <- fec.vary[[pickem]]
   if(tuner == 'nm') tuning.final <- nm.vary[[pickem]]

    return(list(res=res.final,ss.diff = ss.diff.final,tuning = tuning.final))
}
  
  