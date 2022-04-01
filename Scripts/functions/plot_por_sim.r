# plot the results from your simulation.

# Add N.options as arguement here.
# Function arguments
# dat           The output from your por.sim() simulation
# repro.cycle   The number of times an individual reproduces in a year
# year          The initial year for your forward projections, we drop this one from the output, so first year of projections = year+1
# n.mc          The number of simulation draws to run in the standard_SHK() function

#
#


plot.por.sim<-function(dat = por.sim.result, repro.cycle=1, year = 2018,n.mc = 1e4)
{
  pop<-dat$Pop/1000 # I believe this would make it population in thousands, but need to check
  n.sims<-dat$n.sims # Number of simulations from the por.sim() simulation step
  scenarios<-dat$scenarios # Number of removal scenarios to test
  forward.years<-seq(1:dat$forward.n.years)+year # The years of your forward projection (simulation)
  rem.future<-dat$rem.future # The removals for each scenario
  r<-dat$r.nomort.vec # r estimate with no removals
  B1<-dat$B1.vec # realized population growth rate from linear regression
  #ann.change<-dat$ann.change # DK note: commented out as this isn't used 
  #par(mfrow=c(1,1),las=1,omi=c(1,1,0.5,0.25),mar=c(3,3,1,1)) 
  
  ### plot of the medians of the forward projections at different removal scenarios    
  xxx<-do.call(rbind.data.frame,dat$Forward.pop.scenarios) # just unpacking into a dataframe
  names(xxx) <- as.character(forward.years)
  xxx$scenario<-rep(1:scenarios,n.sims) ## The scenarios for each simulation
  med.scen<-xxx %>% group_by(scenario) %>% summarize_at(vars(as.character(forward.years)),median) # The median population size for each year across the simulations
  ## need long data format for plotting
  
  
  ### DK Note: Not used as far as I can tell... 
  #determine median F for input into standard.SHK below.
  xxx<-data.frame(do.call(rbind,dat$F.future.list)) # The F for each scenario and simulation unpacked into a dataframe 
  #names(xxx)<-as.character(dat$rem.future)
  
  # med.F<-xxx %>% summarize_all(median) 
  ## note that this is lazy - you inappropriately calculated medians for the majority of columns
  
  ####### Now we use the standard.SHK function to get r value at the different levels of F, n.mc should be a pretty big number cause it's a monte carlo simulation
  x0<-standard.SHK(n=n.mc,AAFC=1,F=0)
  
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
  # rownames(r.decline)<-c('none',as.character(dat$rem.future))
  #r.decline$biomass<-c(0,as.numeric(dat$rem.future)*39.34/1000)
  #write.csv(r.decline, paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/r.decline.csv'))
  
  write.csv(x0$Results,paste0('D:/Github/ICM/Results/parameter.vals.csv'))
  
  ## get every 5th year on scale
  sim.years.5 <- as.character(rev(seq(max(forward.years),min(forward.years),by=-5)))
  #pdf(paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.pdf'),onefile = TRUE)
  
  ### plots of future projections under different evaluated fishing scenarios.
  scen.plot<-gather(med.scen,year,abundance,all_of(as.character(forward.years)),factor_key = T)
  
  p1<-ggplot(scen.plot) + geom_line(aes(year,abundance/1000,group=scenario,colour=as.factor(scenario))) + 
                          scale_colour_discrete(name='Median Removals',labels = c(rem.future))+labs(y='Total Abundance (1000s)',x='Year')+
                          scale_x_discrete(breaks=sim.years.5,labels=sim.years.5) #+theme(axis.text.x = element_text(angle = 90))
  
  save_plot(plot = p1,filename = 'D:/Github/ICM/Results/Figures/future_projections.png',base_height = 8,base_width = 11)

  #
  ## historical trajectory  # Get the 10 and 90% quantiles and the abundance, put this all into one tidy object
  abund.dat <- data.frame(years = c(N.options$year,forward.years[1]),
                          abund = apply(dat$Pop,2,median),
                          a.90=apply(dat$Pop,2,function(x) quantile(x,0.9)),
                          a.10=apply(dat$Pop,2,function(x) quantile(x,0.1)),
                          removals = c(N.options$N.med,NA)/1000)
  # Backwards simulation abundance time series
  p.abund <- ggplot(abund.dat) + geom_ribbon(aes(ymin=a.10,ymax=a.90,x=years),alpha = 0.1,fill='blue') + 
                               geom_line(aes(x=years,y=abund)) + scale_x_continuous(breaks=seq(1960,2020,by=5),labels = seq(1960,2020,by=5)) + 
                               xlab("") + ylab("Abundance (thousands)") + theme_bw()
  # Backwards simulation removals time series
  p.rem <- ggplot(abund.dat) + geom_line(aes(x=years,y=removals)) + scale_x_continuous(breaks=seq(1960,2020,by=5),labels = seq(1960,2020,by=5)) + 
    xlab("") + ylab("Removals (thousands)") + theme_bw()
  
  p.abund.rem <- plot_grid(p.abund,p.rem,ncol = 1)
  save_plot(plot=p.abund.rem, filename = "D:/Github/ICM/Results/Figures/Abundance_removals_ts.png",base_height = 11,base_width = 8.5)
  
  print('% decline from historical maximum to present')
  print((max(abund.dat$abund)-abund.dat$abund[nrow(abund.dat)])/max(abund.dat$abund))
  print('max to min abundance decline rate')
  print((max(abund.dat$abund)-min(abund.dat$abund))/max(abund.dat$abund))
  print(data.frame(abund.dat$years,abund.dat$abund))
  
  #browser()  
  ### DK Note: Start the why do this Section ###
  ### determine if overfished - allowing variability in SPRmer.
  alpha.hat.x0<-x0$alpha.hat.table # So these are the alpha.hat estimates from the MC simulations
  SPRmer.S0<-(sqrt(alpha.hat.x0)-1)/(alpha.hat.x0-1) # The threshold value representing the depletion of spawners and recruits at Maximum Excess Recruitment (Beverton-Holt assumption)
  scaled.decline<-abund.dat$abund[nrow(abund.dat)]/abund.dat$abund[1] # The scaled decline estimate on average (well median) from the por.sim() simulations
  # DK Note: There is something weird here as this gives just one value, why make a historgram for 1 value? Should the median() go away?
  hist.dat<-data.frame(ind = scaled.decline/median(SPRmer.S0)) # This gives the 'average' view fro the por.sim() simulations, what on average is happening to the stock, only 1 value here...
  upper<-1#-(median(x0$M.cortes[1]))
  # DK Note: I don't see point of this figure or this bit of cod?
  #ggplot(hist.dat) + geom_histogram(aes(x=ind)) + xlim(c(0.2,1)) + ggtitle('Variability in derivation of SPRmer') + geom_vline(xintercept = 0.5,color='blue') + geom_vline(xintercept = upper,color='red')
  # Again because only one value here, what's the point?
  values<-c(length(hist.dat$ind[hist.dat$ind<0.5])/length(hist.dat$ind),length(hist.dat$ind[hist.dat$ind<upper])/length(hist.dat$ind))    
  crit.levels<-c('p<0.5','p<1')
  ### End the why do this section! ###
  
  
  ### determine if overfished - allowing variability in population estimates.
  alpha.hat.x0<-as.data.frame(x0$Results)$Median[which(rownames(x0$Results)=="Alpha hat")]
  SPRmer.S0<-(sqrt(alpha.hat.x0)-1)/(alpha.hat.x0-1) # The threshold value representing the depletion of spawners and recruits at Maximum Excess Recruitment (Beverton-Holt assumption)
  #So this gets the decline for each of the simulations from the por.sim() function, think of this is Current/Unfished, this is numerator in equation 9 from Heather's paper
  #DK Note: So do we just assume first year of our data is 'unfished', that seems likely problematic for a bunch of stocks?  Maybe it don't matter for our purposes?
  sim.dat <- data.frame(scaled.decline=dat$Pop[,ncol(dat$Pop)]/dat$Pop[,1])  
  # So this is Equation 9 from Heathers paper, when this 'ratio' is < p we have a stock that is overfished. According to MSY alone, if this is < 1 we are overfished I believe (see 'upper' below)
  sim.dat$ind <- sim.dat$scaled.decline/SPRmer.S0
  upper<-1#-(median(x0$M.cortes[1]))  ## this corresponds to MSY proxy, rather than MSY proxy minus M
  # So this plots all the realizations from the por.sim() simulation with the median threshold value of SPR/S0 from the Monte Carlo simulations.
  plt.sim <- ggplot(sim.dat) + geom_histogram(aes(x = ind,y=stat(count)/sum(count))) + 
                               xlab("Critical value") + ylab("Frequency")+
                               geom_vline(xintercept = 0.5,color='blue') + geom_vline(xintercept = upper,color='red')
  
  save_plot(plot= plt.sim,filename = "D:/Github/ICM/Results/Figures/high_productivity_overfished.png", base_height = 8.5,base_width = 11)
  
  # So here the value of 0.5 is used as the 'overfishing' p, DK note: I'm not sure why 0.5 is used, just 50% of MSY p?
  values2<-c(length(sim.dat$ind[sim.dat$ind<0.5])/length(sim.dat$ind),length(sim.dat$ind[sim.dat$ind<upper])/length(sim.dat$ind)) ## proportion of sims that ARE overfished (i.e. below p reference point)   
  # Put values and values 2 into a little data.frame and output it to a csv.
  overfishing.eval<-data.frame(crit.levels,var.in.SPRmer=values,var.in.Popsize=values2)
  write.csv(overfishing.eval,'D:/Github/ICM/Results/overfishing.eval.csv')
  
  ## Now we start looking at recovery when projecting forward
  ### table of reaching recovery target under different forward fishing scenarios
  ## first, calculate the values for the decline in the index:
  trajectories<-do.call(rbind.data.frame,dat$Forward.pop.scenarios) # First we grab all the time series for the forward projection
  trajectories$Scenario<-rep(rem.future,length=n.sims) # put a label on them so we know what removal scenario they came from
  names(trajectories) <- c(forward.years,'Scenario')
  #browser()    
  ## set up the columns to sequence
  #test<-trajectories%>%select(c(1,2,seq(3,52,5)))
  #name.breaks<-colnames(test)
  ind.decline<-cbind(trajectories[,which(names(trajectories) != "Scenario")]/dat$Pop[,1],Scenario = trajectories$Scenario) # Now what is the decline (proportional) from the 'unfished' (or at least initial) population abundance
  
  #test$scen<-rep(rem.future,length=n.sims)
  ## use SPRmer from above - note you can't have variability in SPRmer plus in pop size - use pop size
  ind<-cbind(ind.decline[,which(names(ind.decline) != "Scenario")]/SPRmer.S0,Scenario = ind.decline$Scenario)
  # So here we get the proportion that aren't overfished using the 0.5 criteria.
  values3 <- ind %>% dplyr::group_by(Scenario) %>% 
                     dplyr::summarize_all(~sum(.>0.5)/n.sims) 
  # And here we grab every fifth column to look at it.
  
  # Go backwards from final year and take a look at every 5th year.
  values3<- values3 %>% dplyr::select(c(Scenario,all_of(sim.years.5))) ## every 5 years
  # So here we get the proportion that aren't overfished using the MSY criteria (i.e. p = 1). 
  values4 <- ind %>% group_by(Scenario) %>% 
                     summarize_all(~sum(.>upper)/n.sims)
  values4 <- values4 %>% dplyr::select(c(Scenario,all_of(sim.years.5)))

  
  ## output the summarized tables
  write.csv(values3,'D:/Github/ICM/Results/SPRmer.0.5.csv')
  write.csv(values4,'D:/Github/ICM/Results/SPRmer.M.csv')
  
  ### biomass at SPRmer.S0 ref point = 353,000  
  # Here we summarize the results of the por.sim() forward simultion step, taking the annual median abundance
  scen.bmsy.wide <- trajectories %>% dplyr::group_by(Scenario) %>% 
                            dplyr::summarize_at(vars(as.character(forward.years)),median)
  # Initial abundance  
  unfish<-mean(dat$Pop[,1])  ### using historical maximum as unfished abundance
  # Threshold abundance at 'Maximum Excess Recruitment'
  target<-unfish*SPRmer.S0
  # Output results. 
  write.csv(data.frame(unfish=unfish,
                       target=target,
                       SPRmer.S0=SPRmer.S0,
                       hist.decline=(max(abund.dat$abund)-abund.dat$abund[nrow(abund.dat)])/max(abund.dat$abund)),
                       file='D:/Github/ICM/Results/biomass.ref.points.csv')
  
  # Pivot this to be long format for ggplot.
  scen.bmsy.long<-gather(scen.bmsy.wide,year,abundance,all_of(as.character(forward.years)),factor_key = T)
  # So this plot is the Biomass divided by the Biomass at Maximum Excess Recruitment for the different Removals Scenarios in the future projections.
  # DK note: So the geom_hline has SPRmer.S0/SPRmer.S0... why not just put 1??
  pmsy<-ggplot(scen.bmsy.long,aes(year,abundance/unfish/SPRmer.S0,group=Scenario,colour=as.factor(Scenario))) +
                    geom_line() + geom_hline(yintercept=SPRmer.S0/SPRmer.S0) + 
                    scale_colour_discrete(name='Median Removals',labels =rem.future)+labs(y='B/Bsprmer',x='Year') + 
                    scale_x_discrete(breaks=sim.years.5) 
  # And save the plot.
  save_plot(file='D:/Github/ICM/Results/Figures/b_bsprmer.projection.png',plot=pmsy,base_width = 11,base_height = 8.5)
  
}
