# plot the results 


plot.por.sim<-function(n.sims,repro.cycle)
{
  
  #
  pop<-por.sim.result$Pop/1000
  n.sims<-por.sim.result$n.sims
  scenarios<-por.sim.result$scenarios
  forward.years<-seq(1:por.sim.result$forward.n.years)+2018
  rem.future<-por.sim.result$rem.future
  
  r<-por.sim.result$r.nomort.vec
  B1<-por.sim.result$B1.vec
  ann.change<-por.sim.result$ann.change
  #par(mfrow=c(1,1),las=1,omi=c(1,1,0.5,0.25),mar=c(3,3,1,1)) 
  
  ### plot of the medians of the projections at different removal scenarios    
  xxx<-do.call(rbind.data.frame,por.sim.result$Forward.pop.scenarios)
  xxx$scenario<-rep(1:scenarios,n.sims) ## check
  
  med.scen<-xxx %>% group_by(scenario) %>% summarize_at(vars(V1:V52),median)
  ## need long data format for plotting
  values<-names(med.scen[2:length(med.scen)])
  
  
  ### determine median F for input into standard.SHK below.
  xxx<-data.frame(do.call(rbind,por.sim.result$F.future.list))
  #names(xxx)<-as.character(por.sim.result$rem.future)
  
  # med.F<-xxx %>% summarize_all(median) 
  ## note that this is lazy - you inappropriately calculated medians for the majority of columns
  
  ####### get r value at diff levels F
  x0<-standard.SHK(n.sims,1,0,repro.cycle)
  
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
  # rownames(r.decline)<-c('none',as.character(por.sim.result$rem.future))
  #r.decline$biomass<-c(0,as.numeric(por.sim.result$rem.future)*39.34/1000)
  #write.csv(r.decline, paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/r.decline.csv'))
  
  write.csv(x0$Results,paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/parameter.vals.csv'))
  
  ## get every 5th year on scale
  every_nth = function(n) {
    return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
  }
  labs<-forward.years[seq(1,length(forward.years),5)]
  
  #pdf(paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.pdf'),onefile = TRUE)
  
  ### plots of future projections under different evaluated fishing scenarios.
  scen.plot<-gather(med.scen,year,abundance,all_of(values),factor_key = T)
  labs2<-as.character(rem.future)
  p1<-ggplot(scen.plot,aes(year,abundance/1000,group=scenario,colour=as.factor(scenario)))+geom_line()
  p1<-p1+scale_colour_discrete(name='Median Removals',labels = c(labs2))+labs(y='Total Abundance (1000s)',x='Year')
  p1<-p1+scale_x_discrete(breaks=every_nth(n=5),labels=c(labs)) #+theme(axis.text.x = element_text(angle = 90))
  p1
  
  png(file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.future.png'),units='in', width=10, height=6,res=800)
  print(p1)
  dev.off()
  #browser()
  ## historical trajectory 
  years<-c(N.options$year,2019)
  abund<-c()
  a.90<-c()
  a.10<-c()
  for (i in 1:length(years))
  {abund[i]<-median(por.sim.result$Pop[,i]) 
  a.90[i]<-quantile(por.sim.result$Pop[,i],c(0.9))
  a.10[i]<-quantile(por.sim.result$Pop[,i],c(0.1))}
  
  png(file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.abund.png'),units='in', width=12, height=6,res=800)
  par(mfrow=c(2,1),mar=c(4,4,0,4))
  plot(years, abund/1000,type='l',lwd=2,xlab='',ylab='Abundance (thousands)',ylim=c(0,1.4*max(abund/1000))) 
  lines(years,a.90/1000,lty=5)
  lines(years, a.10/1000,lty=5)
  plot(N.options$year,N.options$N.med/1000, type='l',xlab=c('Year'),ylab=c('Removals (thousands)'))
  dev.off()
  
  print('% decline from historical maximum to present')
  print((max(abund)-abund[length(abund)])/max(abund))
  print('max to min abundance decline rate')
  print((max(abund)-min(abund))/max(abund))
  print(data.frame(years,abund))
  
  #browser()  
  ### determine if overfished - allowing variability in SPRmer.
  alpha.hat.x0<-x0$alpha.hat.table
  SPRmer.S0<-(sqrt(alpha.hat.x0)-1)/(alpha.hat.x0-1)
  scaled.decline<-abund[length(abund)]/abund[1]
  ind<-scaled.decline/median(SPRmer.S0)
  upper<-1#-(median(x0$M.cortes[1]))
  hist(ind,breaks=30,main='Variability in derivation of SPRmer',xlab='Critical value',xlim=c(0.2,1))
  abline(v=0.5,col='blue')
  abline(v=upper,col='red')
  
  values<-c(length(ind[ind<0.5])/length(ind),length(ind[ind<upper])/length(ind))    
  crit.levels<-c('p<0.5','p<1')
  
  
  ### determine if overfished - allowing variability in population estimates.
  alpha.hat.x0<-x0$Results[6,1]
  SPRmer.S0<-(sqrt(alpha.hat.x0)-1)/(alpha.hat.x0-1)
  scaled.decline<-por.sim.result$Pop[,ncol(por.sim.result$Pop)]/por.sim.result$Pop[,1]  ### note that this is specific to the N. Atlantic right now.
  ind<-scaled.decline/median(SPRmer.S0)
  upper<-1#-(median(x0$M.cortes[1]))  ## this corresponds to MSY proxy, rather than MSY proxy minus M
  
  png(file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/high.productivity.overfished.png'),units='in', width=9, height=6,res=800)
  hist(ind,breaks=50, xlim=c(0.2,1),main='',xlab='Critical value')
  abline(v=0.5,col='blue')
  abline(v=upper,col='red')
  dev.off()
  
  values2<-c(length(ind[ind<0.5])/length(ind),length(ind[ind<upper])/length(ind)) ## proportion of sims that ARE overfished (i.e. below p reference point)   
  overfishing.eval<-data.frame(crit.levels,var.in.SPRmer=values,var.in.Popsize=values2)
  write.csv(overfishing.eval,paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/overfishing.eval.csv'))
  
  ### table of reaching recovery target under different forward fishing scenarios
  ## first, calculate the values for the decline in the index:
  
  trajectories<-do.call(rbind.data.frame,por.sim.result$Forward.pop.scenarios)
  trajectories$scen<-rep(rem.future,length=n.sims)
  #browser()    
  ## set up the columns to sequence
  #test<-trajectories%>%select(c(1,2,seq(3,52,5)))
  #name.breaks<-colnames(test)
  ind.decline<-trajectories/por.sim.result$Pop[,1]
  #test$scen<-rep(rem.future,length=n.sims)
  ## use SPRmer from above - note you can't have variability in SPRmer plus in pop size - use pop size
  ind<-ind.decline/median(SPRmer.S0)
  ind$scen<-rep(rem.future,length=n.sims)
  values3<-ind%>%group_by(scen)%>%summarize_all(~sum(.>0.5)/n.sims) ## proportion of values > 0.5 i.e. not overfished
  values3<-values3%>%select(c(scen,V2,V7,V12,V17,V22,V27,V32,V37,V42,V47,V52)) ## every 5 years
  names(values3)<-c("Scenario",seq(2020,2070,5)) ## NOTE this is specific to the projection timeframe of 50 years + 2019, 2020
  values4<-ind%>%group_by(scen)%>%summarize_all(~sum(.>upper)/n.sims)
  values4<-values4%>%select(c(scen,V2,V7,V12,V17,V22,V27,V32,V37,V42,V47,V52))
  names(values4)<-c("Scenario",seq(2020,2070,5))
  ## output the summarized tables
  write.csv(values3,paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/SPRmer.0.5.csv'))
  write.csv(values4,paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/SPRmer.M.csv'))
  
  ### biomass at SPRmer.S0 ref point = 353,000  
  test2<-trajectories%>%group_by(scen)%>%summarize_at(vars(V1:V52),median)
  values<-names(test2[2:length(test2)])
  unfish<-mean(por.sim.result$Pop[,1])  ### using historical maximum as unfished abundance
  target<-unfish*median(SPRmer.S0)
  labs<-forward.years[seq(1,length(forward.years),5)]
  
  write.csv(data.frame(unfish=unfish,target=target,SPRmer.S0=SPRmer.S0,hist.decline=(max(abund)-abund[length(abund)])/max(abund)),
            file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/biomass.ref.points.csv'))
  
  png(file=paste0('C:/mydocs/ICCAT/POR2020/NW/_repro',repro.cycle,'/b_bsprmer.projection.png'),units='in', width=9, height=6,res=800)
  scen.bmsy<-gather(test2,year,abundance,all_of(values),factor_key = T)
  labs2<-as.character(rem.future)
  pmsy<-ggplot(scen.bmsy,aes(year,abundance/unfish/SPRmer.S0,group=scen,colour=as.factor(scen)))+geom_line()
  pmsy<-pmsy+scale_colour_discrete(name='Median Removals',labels = c(labs2))+labs(y='B/Bsprmer',x='Year')+geom_hline(yintercept=SPRmer.S0/SPRmer.S0)
  pmsy<-pmsy+scale_x_discrete(breaks=every_nth(n=5),labels=c(labs)) #+theme(axis.text.x = element_text(angle = 90))
  print(pmsy)
  dev.off()
  
  
  
  #      years.f<-c(2019:2069)  
  ## future trajectory
  # abund<-c()
  #a.90<-c()
  #a.10<-c()
  #for (i in 1:por.sim.result$forward.n.years)
  #{abund[i]<-median(por.sim.result$Forward.pop[,i]) 
  #a.90[i]<-quantile(por.sim.result$Forward.pop[,i],c(0.9))
  #a.10[i]<-quantile(por.sim.result$Forward.pop[,i],c(0.1))}
  
  #plot(years.f, abund/1000,type='l',lwd=2,xlab='Years',ylab='Abundance (thousands)',ylim=c(0,1.2*max(abund/1000))) 
  #lines(years.f,a.90/1000,lty=5)
  #lines(years.f, a.10/1000,lty=5)
  
  
  
  #par(mfcol=c(2,2),las=1,omi=c(1,1,0.5,0.25),mar=c(4,4,2,2)) 
  
  #  hist(pop[,length(years)],nclass=20,cex=0.7,xlab=" ",xlim=c(50,350))
  #mtext("N (2019)",1,outer=F,cex=1.2,line=3)
  
  #  hist(pop[,1],nclass=20,cex=0.7,xlab=" ",xlim=c(0,2000))
  #  mtext("N (1950)",1,outer=F,cex=1.2,line=3)
  
  # hist(r,nclass=20,cex=0.7,xlab=" ",plot=T)
  #  mtext("Probability Density",2,outer=T,cex=1.4,line=4)
  #  mtext("r",1,outer=F,cex=1.2,line=3)
  
  # hist(B1,nclass=20,cex=0.7,xlab=" ",plot=T)
  #mtext("Probability Density",2,outer=T,cex=1.4,line=4)
  #mtext("Annual Rate of Change",1,outer=F,cex=1.2,line=3)
  #mtext("Slope",1,outer=F,cex=1.2,line=3)
  
  
}
