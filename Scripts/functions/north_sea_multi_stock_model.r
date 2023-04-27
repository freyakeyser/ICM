# OK, here's I'm going to try and develop a multi=species projection model for the North Sea, we have 10 stocks with good data, the retro fits aren't perfect
# But they mostly aren't a disaster either.  All have a reasonable length of time series as well.
library(tidyverse)
library(GGally)
library(cowplot)

# Download the function to go from inla to sf
funs <- c("https://raw.githubusercontent.com/freyakeyser/ICM/main/Scripts/functions/Lotka_r.r",
          "https://raw.githubusercontent.com/freyakeyser/ICM/main/Scripts/functions/forwrd_sim.r",
          "https://raw.githubusercontent.com/freyakeyser/ICM/main/Scripts/functions/forwrd_project.r"
)
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}

#source("D:/Github/ICM/Scripts/functions/Lotka_r.r") # For testing purposes, delete when done





load(file = "D:/Github/ICM/Results/model_inputs.Rdata")

# Subset to the north sea stocks.  I need years, end abundance, natural mortality, and mx

Stocks <- names(years.tmp)[grep("WGNSSK",names(years.tmp))]
# I'm going to reorder the stocks here so the correlation bit below flows better. Kinda a-posteriorying this :-P 
Stocks <- c("ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus", "ICES-WGNSSK_NS 4-3aN_Trisopterus_esmarkii",
            "ICES-WGNSSK_NS 4-6- 3a_Pollachius_virens", "ICES-WGNSSK_NS 4-7d,20_Gadus_morhua", "ICES-WGNSSK_NS4 _Scopthalmus_maximus",
            "ICES-WGNSSK_NS 7d_Pleuronectes_platessa", "ICES-WGNSSK_NS 4,20_Pleuronectes_platessa","ICES-WGNSSK_NS 7d._Solea_solea","ICES-WGNSSK_NS4 _Solea_solea",
            "ICES-WGNSSK_NS 4-7d_Merlangius_merlangus")
yrs.ns <- years.tmp[Stocks]
nm.ns <- pnm.tmp[Stocks] # Proportional natural mortality here.
mx.ns <- mx.tmp[Stocks]
vpa.ns <- vpa.tmp[Stocks]
ages.ns <- ages.tmp[Stocks]
waa.ns <- waa.tmp[Stocks]
mat.ns <- am.tmp[Stocks]


# Now I'm going to need to align these so they all start and end in the same year. Looks like we have data from 1990:2016 for all stocks
# Which is shorter than I hoped, but oh well.  Off to do some fancy stick handling
yrs.all <- 1990:2016 # Could automate this...
i <- Stocks[1]
yrs.sub <- NULL
nm.sub <- NULL
mx.sub <- NULL
vpa.sub <- NULL
ages.sub <- NULL
waa.sub <- NULL
mat.sub <- NULL
mx.sum <- NULL
for(i in Stocks)
{
  # Subset data
  yrs.t <- yrs.ns[[i]]
  nm.t <- nm.ns[[i]]
  mx.t <- mx.ns[[i]]
  vpa.t <- vpa.ns[[i]]
  ages.t <- ages.ns[[i]]
  waa.t <- waa.ns[[i]]
  mat.t <- mat.ns[[i]]
  # Pick the right data to subset
  sel <- which(yrs.t %in% yrs.all)
  yrs.t <- yrs.t[sel]
  nm.t <- nm.t[sel,]
  mx.t <- mx.t[sel,]
  vpa.t <- vpa.t[sel]
  waa.t <- waa.t[sel,]
  mat.t <- mat.t[sel,]
  # Now pop back into an object
  yrs.sub[[i]] <- yrs.t
  nm.sub [[i]] <- nm.t
  mx.sub[[i]] <- mx.t
  vpa.sub[[i]] <- vpa.t
  ages.sub[[i]] <- ages.t
  waa.sub[[i]] <- waa.t
  mat.sub[[i]] <- mat.t
  #I might want to make this a weighted average at some point, but this is probably fine for the moment
  mx.sum[[i]] <- as.numeric(rowSums(mx.t,na.rm=T))
}  
  
# Now we can look for correlations between fecundity for the different stocks. Going to leave it at that for the moment because Nat.M is irregularly variable over time.
tst <- do.call('rbind',mx.sum)

tst <- t(tst)
colnames(tst) <- c("Melanogrammus aeglefinus","Trisopterus esmarkii",
                   "Pollachius virens", "Gadus morhua", "Scopthalmus maximus",
                   "Pleuronectes platessa 7d", "Pleuronectes platessa 4,20", "Solea solea 7d", "Solea solea 4",
                   "Merlangius merlangus")
tst <- data.frame(tst)
# Now with this we can set up a matrix to make all the stocks mx's depend on each other!
cor.mat <- cor(tst)

p <- ggpairs(data = tst)

save_plot("D:/Github/ICM/Figures/NS_fecundity_correlation_plot.png",p,base_width = 20,base_height = 20)

# So then what do we do? We could take the mx values and start running simulations on the populations
# So from correlation matrix Cod is positively correlated with virens and maximus (these 3 lump)
# platessa 7d is positively correlated with platessa 4,20, solea 4 and 7d ( these 4 lump)
# leaves a weak + correlation with aeglinfuns and esmarkii
# While merlengus isn't correlated with anyone.
# 

n.yrs.proj <- 100
n.sims <- 1000
sd.nm <- 0.1
sd.mx <- 0.1
# Now we can make a super simple or complex climate effect.  As proof of concept this will change either the fecundity
# or natural mortality terms by 1% or 0.5% a decade, based on how I code things below.  Lots of scope to
# do whatever the heck we wanted here.
c.effect <- round(0.25*n.yrs.proj) # The year in which the climate effect starts, 
climate.mx.effect <- 0.01
climate.nm.effect <- 0.01

# These 4 stocks are the ones we use to get the correlation for each 'functional group' (which I define as the ones that are correlated with each other)
focal.stocks <- c("ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus", 
                  "ICES-WGNSSK_NS 4-6- 3a_Pollachius_virens",
                  "ICES-WGNSSK_NS 7d_Pleuronectes_platessa",
                  "ICES-WGNSSK_NS 4-7d_Merlangius_merlangus")
#sd.cor <- 0.01
tmp.mx <- NULL
tmp.nm <- NULL
tmp.mat <- NULL
tmp.age <- NULL
tmp.waa <- NULL
mx.dev <- NULL
nm.dev <- NULL
samp.mx <- NULL
samp.nm <- NULL
ts.unpack <- NULL
r.unpack <- NULL


###################
## OK NOW TO SET THIS UP SO WE LOOP THROUGH MULTIPLE SIMULATIONS!!!
#################

# So everything will need to get wrapped up in a simulation loop
for(j in 1:n.sims)
{
st.time <- Sys.time()
# This picks the years of the mx and nm samples, note all stocks sample from the same year, tho it is separate for mx and nm
# They are all the same length so it doesn't matter what I pick for the nrow() bit.
samp.mx <- round(runif(n.yrs.proj,1,nrow(mx.sub[[1]])))
samp.nm <- round(runif(n.yrs.proj,1,nrow(mx.sub[[1]])))

for(i in Stocks)
{
  mx.tp <- NULL
  nm.tp <- NULL
  age.tp <- NULL
  waa.tp <- NULL
  mat.tp <- NULL
  # OK, so how I am building in the correlation structure for now is to sample the mx and nm structure for all the stocks in a given year
  # So if the above picks year 3 to sample for the mx data we sample year 3 for all mx. This inherently maintains the correlation structure for
  # each stock.  But I do want to put some variability around this so I'm going to make the deviations for the correlated stocks more simlilar (or dissimilar)
  # So this gets us the uncertainty for the year, usually change in +- 7%, which seems reasonable with sd @ 0.1
  
  # Only calculate this when we need to...
  if(i %in% focal.stocks) mx.dev[[i]] <- rlnorm(n.yrs.proj,0,sd.mx) 
  # Calculate this for all the natural mortality bits as we don't build any correlation into this (yet)...
  nm.dev[[i]] <-  rlnorm(n.yrs.proj,0,sd.mx)
  
  # This is bit hacky, but given I know the correlations ahead of doing this and I know the order of the stocks feeding into this, I'll be hacky.
  if(i == "ICES-WGNSSK_NS 4-3aN_Trisopterus_esmarkii")
  {
    row.sel <- which(row.names(cor.mat) == "Trisopterus.esmarkii")
    col.sel <- which(row.names(cor.mat) =="Melanogrammus.aeglefinus")
    corr <- cor.mat[row.sel,col.sel]
    mx.d.tmp <- mx.dev[["ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus"]]-1
    # Toss in a little variability in the correlation, this weakens the correlation a little bit (approx 10%), but does a pretty solid job.
    mx.dev[[i]] <- rnorm(n.yrs.proj,mx.d.tmp,corr)+1
  }# end Functional Group 1 if statement
  
  # More hacky
  if(i %in% c("ICES-WGNSSK_NS 4-7d,20_Gadus_morhua","ICES-WGNSSK_NS4 _Scopthalmus_maximus"))
  {
    if(i == "ICES-WGNSSK_NS 4-7d,20_Gadus_morhua") row.sel <- which(row.names(cor.mat) == "Gadus.morhua")
    if(i == "ICES-WGNSSK_NS4 _Scopthalmus_maximus") row.sel <- which(row.names(cor.mat) == "Scopthalmus.maximus")
    col.sel <- which(row.names(cor.mat) =="Pollachius.virens")
    corr <- cor.mat[row.sel,col.sel]
    mx.d.tmp <- mx.dev[["ICES-WGNSSK_NS 4-6- 3a_Pollachius_virens"]]-1
    # Toss in a little variability in the correlation, this weakens the correlation a little bit (approx 10%), but does a pretty solid job.
    mx.dev[[i]] <- rnorm(n.yrs.proj,mx.d.tmp,corr)+1
  } # end Functional Group 2 if statement
  
  # Final hacky section of this hacky section
  if(i %in% c("ICES-WGNSSK_NS 4,20_Pleuronectes_platessa","ICES-WGNSSK_NS 7d._Solea_solea","ICES-WGNSSK_NS4 _Solea_solea"))
  {
    if(i == "ICES-WGNSSK_NS 4,20_Pleuronectes_platessa") row.sel <- which(row.names(cor.mat) == "Pleuronectes.platessa.4.20")
    if(i == "ICES-WGNSSK_NS 7d._Solea_solea") row.sel <- which(row.names(cor.mat) == "Solea.solea.7d")
    if(i == "ICES-WGNSSK_NS4 _Solea_solea") row.sel <- which(row.names(cor.mat) == "Solea.solea.4")
    col.sel <- which(row.names(cor.mat) =="Pleuronectes.platessa.7d")
    corr <- cor.mat[row.sel,col.sel]
    mx.d.tmp <- mx.dev[["ICES-WGNSSK_NS 7d_Pleuronectes_platessa"]]-1
    # Toss in a little variability in the correlation, this weakens the correlation a little bit (approx 10%), but does a pretty solid job.
    mx.dev[[i]] <- rnorm(n.yrs.proj,mx.d.tmp,corr)+1
  } # end Functional Group 3 if statement

  # Easy enough now to incorporate the above uncertainty to our data...
  tmp.nm[[i]] <- nm.sub[[i]][samp.nm,] * nm.dev[[i]]
  # So we pull out the sample for the year and then add a bit of noise to it.
  tmp.mx[[i]]  <- mx.sub[[i]][samp.mx,]*mx.dev[[i]]
  # the other pieces in case I need them... Just going to take the samples from the fecundity for now, might want to add some uncertainty
  # to the waa someday
  tmp.age[[i]] <- ages.sub[[i]]
  tmp.waa[[i]] <- waa.sub[[i]][samp.mx,]
  tmp.mat[[i]] <- mat.sub[[i]][samp.mx,]

  # Now I can add a climate effect here to change mx and nm.  Ideally it would be driven by stock or functional group, but for first pass
  # We'll just play about with this for everyone...
  if(c.effect < n.yrs.proj)
  {
    ce.mx <-  (climate.mx.effect/10)*(1:(n.yrs.proj-c.effect))
    tmp.mx[[i]][c.effect:n.yrs.proj,] <- tmp.mx[[i]][c.effect:n.yrs.proj,] - ce.mx*tmp.mx[[i]][c.effect:n.yrs.proj,]
    # Need to be careful here since we are on a proportional scale, lets go to instantaneous and convert back
    ce.nm <-  (climate.nm.effect/10)*(1:(n.yrs.proj-c.effect)) # So we'll say this is instantaneous rate.
    tmp.insta.nm <- -log(1-tmp.nm[[i]])
    tmp.nm[[i]][c.effect:n.yrs.proj,] <- 1-exp(-(tmp.insta.nm[c.effect:n.yrs.proj,] + ce.nm*tmp.nm[[i]][c.effect:n.yrs.proj,]))
    
  }
  
  
}  # end for(i in Stocks)

# OK, so I have everything I need now to run the simulations forward.

# So the key here is setting up a good fishing mortality scenario for each stock...
# But I have set it up so that if Abundance drops below 20% of K then F declines to 10% of the mean, basically
# allowing for incidental mortality from other fisheries.
fish.mort <- list(c(0.11,0.1),
                  c(0.35,0.1),
                  c(0.15,0.1),
                  c(0.3,0.1),
                  c(0.15,0.1),
                  c(0.175,0.1),
                  c(0.175,0.1),
                  c(0.06,0.1),
                  c(0.1,0.1),
                  c(0.175,0.1))
names(fish.mort) <- Stocks
res.ts <- NULL
res.r <- NULL
for(s in Stocks)
{
  # The last year of data for all the stocks is 2016
  years <- 2016:(2017+n.yrs.proj-2)
  prop.nat.mort <- tmp.nm[[s]] 
  age.mat <- tmp.mat[[s]]
  mx <- tmp.mx[[s]] 
  vpa.abund <- vpa.sub[[s]]
  weight.age <- tmp.waa[[s]]
  ages <- tmp.age[[s]]
  N.start <- vpa.abund[length(vpa.abund)]
  K <- max(vpa.abund)
  fm <- fish.mort[[s]] # Mean and variance, is log-normally distributed.
  
  tst <- for.sim(years,
                 mat.age = age.mat,
                 nm = prop.nat.mort,
                 w.age = weight.age,
                 ages = ages,
                 rems =  fm,
                 fecund = mx,
                 N.start = N.start,
                 pop.model = 'bounded_exp', K = K,
                 sim= "project",
                 n.sims = 1,
                 sd.mat = 0,
                 sd.nm = 0,
                 sd.wt = 0,
                 sd.fecund = 0)
  #ggplot(tst$Pop) + geom_line(aes(x=years,y=abund)) + geom_hline(yintercept = K)
  

res.ts[[s]] <- data.frame(tst$Pop[,-2],stock = s,sim= j)
res.r[[s]] <- data.frame(tst$r[,-3],stock=s,sim=j)

}

ts.unpack[[j]] <- do.call('rbind',res.ts)
r.unpack[[j]] <- do.call('rbind',res.r)

# Pop a note when done each simulation
timer <- Sys.time() - st.time
print(paste("Simulation ", j))
print(signif(timer,digits=2))

} # end n.sims

ts.final <- do.call("rbind",ts.unpack)
ts.final$fm <- ts.final$removals/ts.final$abund
r.final <- do.call("rbind",r.unpack)
windows(11,11)
ggplot(ts.final) + geom_line(aes(x=years,y=abund,group = sim,color=sim)) + facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) + scale_color_viridis_b()

quants <- ts.final %>%  dplyr::group_by(years,stock) %>% dplyr::summarise(L.50 = quantile(abund,probs=c(0.25),na.rm=T),
                                                                          med = median(abund,na.rm=T),
                                                                          U.50 = quantile(abund,probs=c(0.75),na.rm=T),
                                                                          fml.50 = quantile(fm,probs=c(0.25),na.rm=T),
                                                                          fm = median(fm,na.rm=T),
                                                                          fmu.50 = quantile(fm,probs=c(0.75),na.rm=T))
windows(11,11)
ggplot(quants) + geom_line(aes(x=years,y=med)) + facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) +
                  geom_ribbon(data=quants, aes(x=years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 

