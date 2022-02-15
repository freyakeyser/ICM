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