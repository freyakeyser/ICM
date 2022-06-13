### white shark - evaluation of life history information relative to status.
## started Feb 2018 - HB
# Bring in the fuctions

source("D:/Github/ICM/Scripts/functions/decline_rate.r")

#### simulations

dec.fun()
#####################################################################################################
# to run analysis
#parms.short() or parms.long()
#Baum() or Curtis()
#white.sim.result<-white.sim

## NOTE for MS, only used Curtis decline rate.
set.seed(10)



### create object: white.sim.result











######################

## Manuscript info

#######################

y<-data.frame(r.long=white.sim.LC0$r.summary,r.short=white.sim.SC0$r.summary,
    spr.long=white.sim.LC0$spr.f0.summary,spr.short=white.sim.SC0$spr.f0.summary)

## pre-decline abundance
z<-data.frame(short=round(quantile(white.sim.SC0$Pop[,1],c(0.1,0.5,0.9)),0),long=round(quantile(white.sim.LC0$Pop[,1],c(0.1,0.5,0.9)),0))

## response to decreases in mortality
## note - only using the fast decline scenario here
### also note - forward.n.years has to = 30 for this code to work!!!
response<-data.frame(reduction=c('30%','50%','70%','90%','100%'),
                     prop10.S=c(white.sim.SC3$prop.10,white.sim.SC5$prop.10,white.sim.SC7$prop.10,white.sim.SC9$prop.10,white.sim.SC0$prop.10),
                     propC.S=c(white.sim.SC3$prop.curtis,white.sim.SC5$prop.curtis,white.sim.SC7$prop.curtis,white.sim.SC9$prop.curtis,white.sim.SC0$prop.curtis),
                     
                     prop10.L=c(white.sim.LC3$prop.10,white.sim.LC5$prop.10,white.sim.LC7$prop.10,white.sim.LC9$prop.10,white.sim.LC0$prop.10),
                     propC.L=c(white.sim.LC3$prop.curtis,white.sim.LC5$prop.curtis,white.sim.LC7$prop.curtis,white.sim.LC9$prop.curtis,white.sim.LC0$prop.curtis))

## Status relative to reference points

## Nmin calculations 
Nmin.SC0<-c()
for (i in 1:20)
{Nmin.SC0[i]<-quantile(white.sim.SC0$Pop[,i],c(0.2))  }
Nmin.SC0

#Nmin.SB0<-c()
#for (i in 1:15)
#{Nmin.SB0[i]<-quantile(white.sim.SB0$Pop[,i],c(0.2))  }
#Nmin.SB0

Nmin.LC0<-c()
for (i in 1:20)
{Nmin.LC0[i]<-quantile(white.sim.LC0$Pop[,i],c(0.2))  }
Nmin.LC0

#Nmin.LB0<-c()
#for (i in 1:15)
#{Nmin.LB0[i]<-quantile(white.sim.LB0$Pop[,i],c(0.2))  }
#Nmin.LB0

### PBR, f=0.1, 20th percentile population sizes as above.
ref.point<-data.frame(PBR1.SC0=sum(median(white.sim.SC0$r.vec)*0.5*Nmin.SC0*0.1),
                      #PBR1.SB0=sum(median(white.sim.SB0$r.vec)*0.5*Nmin.SB0*0.1),
                      #PBR1.LB0=sum(median(white.sim.LB0$r.vec)*0.5*Nmin.LB0*0.1),
                      PBR1.LC0=sum(median(white.sim.LC0$r.vec)*0.5*Nmin.LC0*0.1))  ## f=0.1
## for plotting
PBR1.SC0=data.frame(years=c(1970:1989),pbr=(median(white.sim.SC0$r.vec)*0.5*Nmin.SC0*0.1))
#PBR1.SB0=data.frame(years=c(1986:2000),pbr=(median(white.sim.SB0$r.vec)*0.5*Nmin.SB0*0.1))
#PBR1.LB0=data.frame(years=c(1986:2000),pbr=(median(white.sim.LB0$r.vec)*0.5*Nmin.LB0*0.1))
PBR1.LC0=data.frame(years=c(1970:1989),pbr=(median(white.sim.LC0$r.vec)*0.5*Nmin.LC0*0.1))

# median removals by year
rem.SC0<-c()
SC0.90<-c()
SC0.10<-c()
for (i in 1:20)
{rem.SC0[i]<-median(white.sim.SC0$removals[,i]) 
SC0.90[i]<-quantile(white.sim.SC0$removals[,i],c(0.9))
SC0.10[i]<-quantile(white.sim.SC0$removals[,i],c(0.1))}
SC0<-data.frame(removals=rem.SC0,q.90=SC0.90,q.10=SC0.10)
SC0$sim<-'Short lifespan; slow decline'
SC0$year<-c(1970:1989)

rem.SB0<-c()
SB0.90<-c()
SB0.10<-c()
for (i in 1:15)
{rem.SB0[i]<-median(white.sim.SB0$removals[,i]) 
SB0.90[i]<-quantile(white.sim.SB0$removals[,i],c(0.9))
SB0.10[i]<-quantile(white.sim.SB0$removals[,i],c(0.1))}
SB0<-data.frame(removals=rem.SB0,q.90=SB0.90,q.10=SB0.10)
SB0$sim<-'Short lifespan; fast decline'
SB0$year<-c(1986:2000)

rem.LC0<-c()
LC0.90<-c()
LC0.10<-c()
for (i in 1:20)
{rem.LC0[i]<-median(white.sim.LC0$removals[,i]) 
LC0.90[i]<-quantile(white.sim.LC0$removals[,i],c(0.9))
LC0.10[i]<-quantile(white.sim.LC0$removals[,i],c(0.1))}
LC0<-data.frame(removals=rem.LC0,q.90=LC0.90,q.10=LC0.10)
LC0$sim<-'Long lifespan; slow decline'
LC0$year<-c(1970:1989)

rem.LB0<-c()
LB0.90<-c()
LB0.10<-c()
for (i in 1:15)
{rem.LB0[i]<-median(white.sim.LB0$removals[,i]) 
LB0.90[i]<-quantile(white.sim.LB0$removals[,i],c(0.9))
LB0.10[i]<-quantile(white.sim.LB0$removals[,i],c(0.1))}
LB0<-data.frame(removals=rem.LB0,q.90=LB0.90,q.10=LB0.10)
LB0$sim<-'Long lifespan; fast decline'
LB0$year<-c(1986:2000)

### removals relative to PBR
library(ggplot2)
library(reshape2)


### make 4 ggplots with inset panels
# combine with grid.arrange

p1<-ggplot(SB0, aes(year,y=removals))+geom_line()
g<-ggplotGrob(ggplot(PBR1.SB0,aes(years,pbr))+geom_line()+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())+scale_y_continuous(limits=c(0,25)))
p1<-p1+annotation_custom(grob=g,xmin=1970,xmax=1985,ymin=1500,ymax=3000)
p1<-p1+geom_ribbon(aes(ymin=SB0$q.10,ymax=SB0.90),alpha=0.2)
p1<-p1+scale_x_continuous(limits=c(1970,2000))+scale_y_continuous(limits=c(0,3000))
p1<-p1+xlab('')+ylab('Removals (# animals)')+ggtitle('Short lifespan; fast decline')
p1<-p1+theme(axis.title.y = element_text(size=rel(1.2)))

p2<-ggplot(SC0, aes(year,y=removals))+geom_line()
g<-ggplotGrob(ggplot(PBR1.SC0,aes(years,pbr))+geom_line()+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())+scale_y_continuous(limits=c(0,26)))
p2<-p2+annotation_custom(grob=g,xmin=1981,xmax=1989,ymin=1200,ymax=2000)
p2<-p2+geom_ribbon(aes(ymin=SC0$q.10,ymax=SC0.90),alpha=0.2)
p2<-p2+scale_x_continuous(limits=c(1970,1990))+scale_y_continuous(limits=c(0,2000))
p2<-p2+xlab('')+ylab('Removals (# animals)') ##+ggtitle('Short lifespan; slow decline')
p2<-p2+theme(axis.title.y = element_text(size=rel(1.5)))

p3<-ggplot(LB0, aes(year,y=removals))+geom_line()
g<-ggplotGrob(ggplot(PBR1.LB0,aes(years,pbr))+geom_line()+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())+scale_y_continuous(limits=c(0,25)))
p3<-p3+annotation_custom(grob=g,xmin=1970,xmax=1985,ymin=1500,ymax=3000)
p3<-p3+geom_ribbon(aes(ymin=LB0$q.10,ymax=LB0.90),alpha=0.2)
p3<-p3+scale_x_continuous(limits=c(1970,2000))+scale_y_continuous(limits=c(0,2000))
p3<-p3+xlab('Year')+ylab('Removals (# animals)')+ggtitle('Long lifespan; fast decline')
p3<-p3+theme(axis.title.x =  element_text(size=rel(1.2)))
p3<-p3+theme(axis.title.y = element_text(size=rel(1.2)))

p4<-ggplot(LC0, aes(year,y=removals))+geom_line()
g<-ggplotGrob(ggplot(PBR1.LC0,aes(years,pbr))+geom_line()+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())+scale_y_continuous(limits=c(0,26)))
p4<-p4+annotation_custom(grob=g,xmin=1981,xmax=1989,ymin=1000,ymax=1700)
p4<-p4+geom_ribbon(aes(ymin=LC0$q.10,ymax=LC0.90),alpha=0.2)
p4<-p4+scale_x_continuous(limits=c(1970,1990))+scale_y_continuous(limits=c(0,1700))
p4<-p4+xlab('Year')+ylab('Removals (# animals)')##+ggtitle('Long lifespan; slow decline')
p4<-p4+theme(axis.title.x =  element_text(size=rel(1.5)))
p4<-p4+theme(axis.title.y = element_text(size=rel(1.5)))

## build the figure
library(gridExtra)
Fig2<-grid.arrange(p2,p4)
Fig2

## high res
png(file="Fig2.png",units='in', width=6, height=8,res=800)
grid.arrange(p2,p4)
dev.off()

#############
junk1<-white.sim.LC0$Forward.pop  ## project forward 100 years
junk2<-white.sim.SC0$Forward.pop

junk11<-ifelse(junk1[,]>2*junk1[,1],0,1)
junk22<-ifelse(junk2[,]>2*junk2[,1],0,1)

white.sim.LC0$doubling<-rowSums(junk11)
white.sim.SC0$doubling<-rowSums(junk22)


# build data
test<-data.frame(lifespan=c(rep('long',4000),rep('short',4000)),
                scenario= rep(c(rep('Growth rate (r)',1000),rep('Instantaneous fishing mortality (F.crit)',1000), rep('Population doubling time',1000),rep('Lifetime reproductive output',1000)),2),
           value=c(white.sim.LC0$r.vec, white.sim.LC0$F.crit.vec,white.sim.LC0$doubling,white.sim.LC0$spr.f0.vec,
                   white.sim.SC0$r.vec, white.sim.SC0$F.crit.vec, white.sim.SC0$doubling,white.sim.SC0$spr.f0.vec))
## get rid of the weird doubling-time estimates
test<-test[!test$value%in%c(4,5,6,100),]

### create medians
library(dplyr) # With dplyr for example
test <- test %>% group_by(scenario,lifespan) %>%
  mutate(med = median(value))

### histograms of parameters.
p1<-ggplot(test,aes(x=value, fill=lifespan,color=lifespan))+
  geom_histogram(position="identity",alpha=0.85)+geom_vline(aes(xintercept=med))+ 
  facet_wrap(~scenario,scales='free')+ scale_fill_grey()+labs(x = 'Value', y = 'Count')

## high res
png(file="Fig1.png",units='in', width=6, height=6,res=800)
p1
dev.off()

       
#############################################

period<-c(1970:2019)
zero.mort<-c(apply(white.sim.LC0$Pop,2,FUN=median),apply(white.sim.LC9$Forward.pop[,1:30],2,FUN=median))
zero.mort1<-c(apply(white.sim.LC0$Pop,2,FUN=quantile, prob=0.1),apply(white.sim.LC9$Forward.pop[,1:30],2,FUN=quantile, prob=0.1))
zero.mort9<-c(apply(white.sim.LC0$Pop,2,FUN=quantile, prob=0.9),apply(white.sim.LC9$Forward.pop[,1:30],2,FUN=quantile, prob=0.9))

mort3<-c(apply(white.sim.LC0$Pop,2,FUN=median),apply(white.sim.LC3$Forward.pop[,1:30],2,FUN=median))
mort31<-c(apply(white.sim.LC0$Pop,2,FUN=quantile, prob=0.1),apply(white.sim.LC3$Forward.pop[,1:30],2,FUN=quantile, prob=0.1))
mort39<-c(apply(white.sim.LC0$Pop,2,FUN=quantile, prob=0.9),apply(white.sim.LC3$Forward.pop[,1:30],2,FUN=quantile, prob=0.9))

mort5<-c(apply(white.sim.LC0$Pop,2,FUN=median),apply(white.sim.LC5$Forward.pop[,1:30],2,FUN=median))
  
mort7<-c(apply(white.sim.LC0$Pop,2,FUN=median),apply(white.sim.LC7$Forward.pop[,1:30],2,FUN=median))

D<-data.frame(period=period,zero.mort=zero.mort,zero.mort1=zero.mort1,zero.mort9=zero.mort9,
              mort3=mort3,mort31=mort31,mort39=mort39)

#plot(period, zero.mort, type='l',ylim=c(0,15000))

fig3<-ggplot(D, aes(period,y=zero.mort))+geom_line(cex=1)
fig3<-fig3+geom_ribbon(aes(ymin=zero.mort1,ymax=zero.mort9),alpha=0.2)
fig3<-fig3+xlab('Year')+ylab('Abundance') 
fig3<-fig3+theme(axis.title.y = element_text(size=rel(1.2)))+theme(axis.title.x = element_text(size=rel(1.2)))
fig3<-fig3+geom_line(aes(period, mort3),lty=5,cex=1)+geom_ribbon(aes(ymin=mort31,ymax=mort39),alpha=0.4)
fig3

Db<-data.frame(period=c(rep(period,4)),median=c(zero.mort,mort3,mort5,mort7),scenario=c(rep('90%',length=length(period)),rep('30%',length=length(period)),rep('50%',length=length(period)),rep('70%',length=length(period))))
  
fig3b<-ggplot(Db, aes(period,y=median,lty=scenario))+geom_line(cex=1)
fig3b<-fig3b+xlab('Year')+ylab('Abundance') 
fig3b<-fig3b+theme(axis.title.y = element_text(size=rel(1.2)))+theme(axis.title.x = element_text(size=rel(1.2)))
fig3b


## high res
png(file="Fig3b.png",units='in', width=6, height=6,res=800)
fig3b
dev.off()










#####################################################

## old summary plots

#######################################################
#### if need to store multiple panels in a list and then combine
p.list = lapply(sort(unique(all.rem$sim)), function(i) {
  ggplot(all.rem[all.rem$sim==i,], aes(Var2,median(value))) +geom_smooth(se=FALSE)
})
p2.list<-lapply(sort(unique(ref.annual$sim)), function(i) {
  ggplotGrob(ggplot(ref.annual[ref.annual$sim==i,], aes(years,PBR)) +geom_line())
})


plot.white.sim<-function()
{
  pop<-white.sim.result$Pop
  remov<-white.sim.result$removals[,2:length(years)]
  
  r<-white.sim.result$r.vec
  B1<-white.sim.result$B1.vec
  ann.change<-white.sim.result$ann.change
  
  #### all simulations
  windows()
  #par(mfcol=c(2,1),mar=c(5,5,2,2))
  #plot(years,pop[1,],type="l",ylim=c(0,max(pop[,1])),xlab="Year",ylab="Abundance",lty=1)
  #for (i in 2:n.sims)
  #{lines (years,pop[i,])}
  
  ## boxplots
  #library(ggplot2)
  #library(reshape2)
  #remov.m <- melt(remov)
  #ggplot(remov.m, aes(as.factor(Var1), value)) + geom_boxplot()+
  #  stat_summary(fun.y=mean, geom="line", aes(group=1),col='red',lwd=2)
  ## can't really see the mean line anyways - too many simulations.
  
  boxplot(t(remov),medcol="red",whisklty = 0, staplelty = 0,xlab="Simulation number",ylab='Removals (# of animals)')
  
  windows()
  par(mfcol=c(2,2),omi=c(0.5,0.7,0.5,0.25),mar=c(4,4,2,2)) 
  
  hist(pop[,length(years)],nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  mtext("N (2000)",1,outer=F,cex=1.2,line=3)
  abline(v=median(pop[,length(years)]),col='red')
  
  hist(pop[,1],nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  mtext("N (1986)",1,outer=F,cex=1.2,line=3)
  abline(v=median(pop[,1]),col='red')
  
  
  hist(r,nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  #mtext("Probability Density",2,outer=T,cex=1.4,line=2)
  mtext("r",1,outer=F,cex=1.2,line=3)
  abline(v=median(r),col='red')
  
  ## instead of the annual rate of change, it would be useful to plot median removals.
  hist(rowMeans(remov),nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  abline(v=median(rowMeans(remov)),col='red')
  mtext("Mean Number of Removals",1,outer=F,cex=1.2,line=3)
  
  #hist(B1,nclass=20,probability=T,cex=0.7,xlab=" ",plot=T,main="")
  #mtext("Annual Rate of Change",1,outer=F,cex=1.2,line=3)
  #mtext("Slope",1,outer=F,cex=1.2,line=3)
  #abline(v=median(B1),col='red')
  
}

