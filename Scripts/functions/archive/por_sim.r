# This is the primary simulation function where everything needed is set up and feed into the other functions. Any parameters required for the underlying functions are going
# to need to get loaded here....


## Arguments
#1: n.sims          The number of simulations run
#2:  q.ext          Population extinction threshold, basically a value after which the population is probably screwed
#3:  sigma          Used if you are assuming autocorrelated values of r for the population.  Set to 0 means r is random each year
#4:  u.long         This will select more (but still not all how it is currently coded) of the exploitation data for the analysis
#5:  repro.cycle    The number of times an individual reproduces in a year
#6:  N.options      Time series of removals.  The only things currently used here is the N.med and year columns but I believe you can tweak code to use other columns.


por.sim<-function(n.sims,q.ext=500,sigma=0,u.long=F,repro.cycle,N.options = N.options)
{
  #assign output to "por.sim.result" which is an object used by the plotting function below

  set.seed(204) #so that different scenarios are directly comparable

  #setup windows progress bar
  pb <- winProgressBar(title="Simulation progress", label="0% done", min=0, max=100, initial=0)

  # Initialize a shit tonne of different objects
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
  age.mat.vec<-rep(0,n.sims)
  max.age.vec<-rep(0,n.sims)
  mx.list<-list(rep(NA,length(n.sims)))
  r.M.list<-list(rep(NA,length(n.sims)))


##################################################### Section Function arguments - The stuff in here should become arguments to a function ################################################


  # Now we start to set up our input data how we want
  # Here we use N.med as our removals
  removals<-c(N.options$N.med,0)
  #removals<-removals*2
  # Then we take mean removals, DK Note: we will want to uncomment this line and comment the Porbeagle example once this is generalized and up and running
  #mean.removals<-mean(N.options$N.med,na.rm=T)
  #For the Porbeagle example everything after 2009 is tossed because of new regulations.
  mean.removals<-mean(N.options$N.med[N.options$year<2010])  ## NOT USED: 2009 was the year that removals of porbeagle REALLY dropped.

  sel <- 1 ## DK Note: basically vulnerable to the fishery right after birth, in the below code it's hard coded to be age 2 on wards we are going to need to generalize this and make a vector.
  r.cutoff<-0.2 #This sets a maximum value of the Lotka.r (maximum rate of population growth, 0.2 would be incredibly high growth rate
  K<-20e6  # The carrying capacity.  20e6 is super high for the Porbeagle (I believe this effectively makes the model an exponential growth model) this is used for forward projections.  
  #DK note, carrying capacity always hard to estimate but do we care if we are just doing backwards simulation anyways...

  # So this is the number of individuals at the end of the time series.
  #logic here is approx. 200000 sharks currently - from SCA in Campana et al. 2010.
  # sampling done this way so that half are above and below this estimate.
  # DK note -- How would we specify this in the absence of an assessment, if we use the last value of the assessment, then our result isn't independent of the assessment
  # which I think we want it to be. hmm... Also this current code breaks if you use n.sims that isn't an even number!
  N.end<-c(sample(100000:200000,size=n.sims/2,replace=T),sample(200000:300000,size=n.sims/2,replace=T))
  #N.end<-c(sample(150000:200000,size=n.sims/2,replace=T),sample(200000:250000,size=n.sims/2,replace=T))
  # Removals scenarios. 
  # DK note: We'll need to make this an argument for the function so we can alter these scenarios, we'll probably want an option to use these, or use the above 'junk' code 
  # to base future removals on historical removals.
  rem.future<-c(0,seq(1000,24000,by=1000))  ### added zero removals as a scenario
  scenarios=length(rem.future)
  # So this gets the mean number of individuals removed in the final 3 years of data
  # DK note: We may want to use something other than the removals in the final 3 years, so might make this an option..
  R.ave<-mean(N.options$N.med[seq(nrow(N.options)-2,nrow(N.options))])

  # Need to include parameters from the sub-functions in the call to this function.  Or to simplify and have a file we load with all the parameters we need.
  ##################################################### End Section Function arguments                                                        ################################################


  #####################################################  Starting to run the model     #####################################################

  #create parameters -  consistent with estimation method used in standard.SHK (which isn't used in this currently)
  ## note that lotka.r function just removes necessity of developing the survival by age matrix
  ## and is used here as a short-cut. Parameter estimates are within three decimal places of rmax
  # from standard.SHK
  for (i in 1:n.sims)
  {
    junk<-parms.calc()
    age.mat.vec[i]<-junk$age.mat
    max.age.vec[i]<-junk$max.age
    mx.list[[i]]<-junk$mx
    r.M.list[[i]]<-junk$r.M
  }

  # Here we jump into the F crit calculations and the lotka.r calcs.
  # 2 is a patch to get n.sims values for r that are <r.cutoff  # DK note -- Don't understand this comment yet....
  for(i in 1:(1*n.sims))
  {
    #print(junk$par)
    # Function to get f critical and N critical
    junk2<-F.crit(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],sel,mean.removals)
    F.crit.vec[i]<-junk2$f
    N.crit.vec[i]<-junk2$N.crit
    #print(junk2$spr.f0)
    # Function to get lotka R value when there are no removals.
    junk3<-lotka.r(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],0,sel) #no human induced mortality
    # Now we have lotka r value.
    r.nomort.vec[i]<-junk3$par
  }

  ############################################################

  #then do backward population projections
  for(i in 1:n.sims)
  {
    # Initiate progress bar
    info <- sprintf("%d%% done", round((i/n.sims)*100))
    setWinProgressBar(pb, i/(n.sims)*100, label=info)
    # So this is important, we need some estimate of the numbers at the end of the time series
    Pop.vec[n.years]<-N.end[i]
    # If there was no fishery what would the population be (note how we add removals back in here)
    for(y in 1:(n.years-1)) Pop.vec[n.years-y]<-(Pop.vec[n.years-y+1]+removals[n.years-y])/exp(r.nomort.vec[i])  ## productivity in the absense of fishing

    # Make a matrix for each simulation run
    Pop[i,]<-Pop.vec
    #print(Pop.vec)
    
    # Now we grab the slope of the change in abundance over time on log (good!) scale. DK Note: would we ever care about the uncertainty here?
    B1.vec[i]<-lm(log(Pop.vec)~years)$coef[2]

    #calculate historical exploitation - note: as is, this junk and the mean.u.long/short.vec are not used
    # DK note: They would be used if you wanted to use the historical exploitation to predict future exploitation (which we many want to do)
    junk<-c()

    # first calc is relative to observed removals. So we get the estimated fishing mortality for each year in current simulation
    for(y in 1:n.years) junk[y]<-u.calc(age.mat.vec[i],max.age.vec[i],r.M.list[[i]],sel,removals[y],Pop.vec[y])

    # Take the mean of this, DK note: We'll need to tweak these as they remove rather specific years.
    mean.u.long.vec[i]<-mean(junk[(length(junk)-15):length(junk)-11])  # NOT USED: years are 2004 - 2008;
    mean.u.short.vec[i]<-mean(junk[(length(junk)-10):length(junk)-1])  # NOT USED; years are 2009-2018

    ### this is what is used to evaluate future removals.
    ## calculate exploitation relative to current pop size for hypothetical future removal scenarios
    junkx<-c()

    ## this list is relative to hypothetical levels we are interested in testing
    for (k in 1:scenarios) junkx[k]<-u.calc(age.mat.vec[i],max.age.vec[i],r.M.list[[i]],sel,rem.future[k],Pop.vec[length(Pop.vec)])

    # Here we get the exploitation rate under different removals scenarios
    u.future.list[[i]]<-junkx
    F.future.list[[i]]<- -log(1-junkx)
  } # end the n.sims


  ### evaluated removals table
  #then do forward projections:
  for(i in 1:n.sims)
  {

    #first recalc r with fishing mortality for each of your removal sceanrios from each of the simulations, using the proportion removed by fishery.
    xxx<-as.vector(u.future.list[[i]])

    x1<-c()
    x2<-c()
    temp<-c()
    junkxx<-c()
    # Get the estimates of r and N.crit for each removals scenario
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
    #DK note: We may want to use this if we want to use past removals to predict future removals, see code in simulation loop.
    if(u.long==T)
    {
      junk1<-lotka.r(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],mean.u.long.vec[i],sel)
    }
    else { junk1<-lotka.r(age.mat.vec[i],max.age.vec[i],mx.list[[i]],r.M.list[[i]],mean.u.short.vec[i],sel)}

    # The maximum growth rate from the Lotka calculations.
    r.allmort.vec[i]<-junk1$par

    ## note that this also isn't used. Incorporates autocorrelated deviates to annual values for r in future projections. I would need to make a matrix that the deviates were applied over
    ## and I just didn't have time.
    #random deviate vectors
    # DK note, question for Heather if she thinks getting this up and running is worth adding
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
        if(r.forward[y]>r.cutoff){r.forward[y]<-r.cutoff} #sets an upper bound on r in the sims
      } # end for (y in 2:forward.n.years)
    } # end if(sigma>0)
    #  print(summary(r.forward))

    #then project forward
    ## 2019
    # So to get the first year abundance we take the number at the end from our simulation, grow it at 'r' and then remove R.ave from that.
    # DK note: So assumption here is population grows and then the removals come after growth (so effectively at the end of the year)
    Forward.pop.vec[1]<-N.end[i]*exp(r.nomort.vec[i])-R.ave
    r.first<-r.nomort.vec[i]

    #print(r.first)
    ## I am being lazy here by not incorporating future variability in r
    ## DK note: Do we not want to be lazy and put some uncertainty around r?
    xx<-as.vector(r.scenarios[[i]])
    # Set up an empty matrix, each row will be the N in each year for a given removal sceanrio
    Forward.pop.mat<-matrix(rep(0,forward.n.years*length(xx)),length(xx),forward.n.years)
    # Get the starting population, just using the same number for each removal scenario for simplicity.
    Forward.pop.mat[,1]<-N.end[i]*exp(r.nomort.vec[i])-R.ave
    # And run through the scenarios
    for(y in 1:(forward.n.years-1))
    {
      if(y==1)
      {
        # Year 2
        Forward.pop.vec[y+1]<-Forward.pop.vec[y]*exp(r.first)*(1-(Forward.pop.vec[y]/K))-R.ave
        Forward.pop.mat[,y+1]<-(Forward.pop.mat[,y]*exp(r.first))*(1-(Forward.pop.mat[,y]/K))-R.ave
      } # end if(y==1)
      #and beyond
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
      } # end else
    } # end for(y in 1:....)
    Forward.pop[i,]<-Forward.pop.vec
    Forward.pop.scenarios[[i]]<-Forward.pop.mat
    #print(Forward.pop.mat)
  } #end simulation
  close(pb)

  #=========================================================================
  # Summarize projections

  ### note that a lot of these summary values aren't used in the current assessment.
  ## they are hold-overs from previous use of the code.
  
  end.pop.size<-Forward.pop[,forward.n.years] # Population size at the end of the simulations

  r.nomort.summary<-quantile(r.nomort.vec,c(0.1,0.5,0.9),na.rm=T)  #no human induced mortality, summarized from each simulation
  r.allmort.summary<-quantile(r.allmort.vec,c(0.1,0.5,0.9),na.rm=T)  #This currently isn't output that is used, but this includes both "other" and known removals
  mean.u.long.summary<-quantile(mean.u.long.vec,c(0.1,0.5,0.9),na.rm=T) #This currently isn't output that is used, but would feed the r.allmort calculation when u.long=T
  mean.u.short.summary<-quantile(mean.u.short.vec,c(0.1,0.5,0.9),na.rm=T) #This currently isn't output that is used, but would feed the r.allmort calculation when u.long=F
  F.crit.summary<-quantile(F.crit.vec,c(0.1,0.5,0.9),na.rm=T) # Summary of the F.crit for each simulation
  N.crit.summary<-quantile(N.crit.vec,c(0.1,0.5,0.9),na.rm=T) # Summary of the N.crit for each simulation
  N.back.init.summary<-quantile(Pop[,1],c(0.1,0.5,0.9),na.rm=T)    # Initial population estimates from the backward simulation estimates
  N.for.init.summary<-quantile(Forward.pop[,2],c(0.1,0.5,0.9),na.rm=T) # Initial population estimates from the forward simulation projections
  N.for.end.summary<-quantile(Forward.pop[,forward.n.years],c(0.1,0.5,0.9),na.rm=T) # Final population estimates from the forward simulation projections

  proportion.declining<-length(B1.vec[B1.vec<0])/length(B1.vec) # what proportion of the simulations project declines
  end.pop.summary<-quantile(end.pop.size,c(0.1,0.5,0.9),na.rm=T) # Summary of the population size at the end of each simulation
  prop.ext<-length(end.pop.size[end.pop.size<q.ext])/length(end.pop.size) # The number of simulations in which the population is below the extinction threshold

  return(list(N.crit.list=N.crit.list,
              u.future.list=u.future.list,
              F.future.list=F.future.list,
              r.scenarios=r.scenarios,
              rem.future=rem.future,
              Forward.pop.scenarios=Forward.pop.scenarios,
              scenarios=scenarios,
              n.sims=n.sims,
              years=years, 
              q.ext=q.ext,
              sigma=sigma,
              mean.removals=mean.removals,
              Pop=Pop,
              Forward.pop=Forward.pop,
              forward.n.years=forward.n.years,
              r.allmort.vec=r.allmort.vec,
              r.nomort.vec=r.nomort.vec,
              B1.vec=B1.vec,
              mean.u.long.vec=mean.u.long.vec,
              r.nomort.summary=r.nomort.summary,
              r.allmort.summary=r.allmort.summary,
              F.crit.summary=F.crit.summary,
              N.crit.summary=N.crit.summary,
              N.back.init.summary=N.back.init.summary,
              mean.u.long.summary=mean.u.long.summary,
              mean.u.short.summary=mean.u.short.summary,
              N.for.init.summary=N.for.init.summary,
              N.for.end.summary=N.for.end.summary,
              proportion.declining=proportion.declining,
              end.pop.summary=end.pop.summary,
              prop.ext=prop.ext,
              ext.time=ext.time))

} #end function