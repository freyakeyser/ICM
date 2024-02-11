# OK, here's I'm going to try and develop a multi=species projection model for the North Sea, we have 10 stocks with good data, the retro fits aren't perfect
# But they mostly aren't a disaster either.  All have a reasonable length of time series as well.

#################  Section 1 Loading #################  Section 1 Loading #################  Section 1 Loading  ###############################################
library(tidyverse)
library(GGally)
library(cowplot)
library(ggthemes)
# Set the base plot theme
theme_set(theme_few(base_size = 14))

# Download the function to go from inla to sf
funs <- c("https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/Lotka_r.r",
          "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_sim.r",
          "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_project.r"
)
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}

#source("D:/Github/ICM/Scripts/functions/Lotka_r.r") # For testing purposes, delete when done
#load(file = "C:/Users/Owner/Documents/Github/ICM/Results/model_inputs.Rdata")
#load(file = "D:/Github/ICM/Results/model_inputs.Rdata")
loc <- 'D:/GitHub/ICM'
load(file = paste0(loc,"/Results/model_inputs.Rdata"))

load(file = paste0(loc,"/Results/all_cleaned_forward_tune_summaries_fec_nm.Rdata"))
Stocks <- Stocks[grep("WGNSSK",Stocks)] # The stocks we want
########################### End Section 1 Loading  ###################################### End Section 1 Loading ###############################################




########################## Section 2 Parameters ########################## Section 2 Parameters ########################## Section 2 Parameters

yrs.all <- 1990:2016 # These are the years we have data for all 10 stocks
n.yrs.proj <- 25
n.sims <- 5
sd.nm <- 0.1
sd.mx <- 0.1
# Now we can make a super simple or complex climate effect.  As proof of concept this will change either the fecundity
# or natural mortality terms by 1% or 0.5% a decade, based on how I code things below.  Lots of scope to
# do whatever the heck we wanted here.
c.effect <- round(0.25*n.yrs.proj) # The year in which the climate effect starts, 
climate.mx.effect <- 0
climate.nm.effect <- 0




########################## Section 3 Inter-Species Correlation ########################## Section 3 Inter-Species Correlation ########################## 


yrs.sub <- NULL
nm.sub <- NULL
mx.sub <- NULL
vpa.sub <- NULL
ages.sub <- NULL
waa.sub <- NULL
mat.sub <- NULL
mx.sum <- NULL
fm.sub <- NULL
for(i in Stocks)
{
  # Subset data
  yrs.t <- for.tune.all[[i]]$res$year
  nm.t <- 1-exp(-for.tune.all[[i]]$nm.opt)
  mx.t <- for.tune.all[[i]]$fecund.opt
  vpa.t <- for.tune.all[[i]]$res$est.abund
  ages.t <- ages.tmp[[i]]
  waa.t <- waa.tmp[[i]]
  mat.t <- am.tmp[[i]]
  rem.t <- rem.tmp[[i]]
  # Pick the right data to subset
  sel <- which(yrs.t %in% yrs.all)
  yrs.t <- yrs.t[sel]
  nm.t <- nm.t[sel,]
  mx.t <- mx.t[sel,]
  vpa.t <- vpa.t[sel]
  waa.t <- waa.t[sel,]
  mat.t <- mat.t[sel,]
  rem.t <- rem.t[sel,]
  fm.t <- data.frame(years = rem.t$Year,fm = rem.t$rem/vpa.t, stock =i)
  # Now pop back into an object
  yrs.sub[[i]] <- yrs.t
  nm.sub [[i]] <- nm.t
  mx.sub[[i]] <- mx.t
  vpa.sub[[i]] <- vpa.t
  ages.sub[[i]] <- ages.t
  waa.sub[[i]] <- waa.t
  mat.sub[[i]] <- mat.t
  fm.sub[[i]] <- fm.t
  #I might want to make this a weighted average at some point, but this is probably fine for the moment
  mx.sum[[i]] <- as.numeric(rowSums(mx.t,na.rm=T))
}  
  
# Now we can look for correlations between fecundity for the different stocks. Going to leave it at that for the moment because Nat.M is irregularly variable over time.
mx.cor.dat <- do.call('rbind',mx.sum)

mx.cor.dat <- t(mx.cor.dat)
# Make the names nice for the plot
colnames(mx.cor.dat) <- c("Melanogrammus aeglefinus","Trisopterus esmarkii",
                   "Pollachius virens", "Gadus morhua", "Scopthalmus maximus",
                   "Pleuronectes platessa 7d", "Pleuronectes platessa 4,20", "Solea solea 7d", "Solea solea 4",
                   "Merlangius merlangus")
mx.cor.dat <- data.frame(mx.cor.dat)
# Now with this we can set up a matrix to make all the stocks mx's depend on each other!
mx.cor.mat <- cor(mx.cor.dat)

p.mx.cor <- ggpairs(data = mx.cor.dat)
#save_plot("D:/Github/ICM/Figures/NS_sims/NS_fecundity_correlation_plot.png",p.mx.cor,base_width = 22,base_height = 22)

# So from correlation matrix Cod is positively correlated with virens and maximus (these 3 lump)
# platessa 7d is positively correlated with platessa 4,20, solea 4 and 7d ( these 4 lump)
# leaves a weak + correlation with aeglinfuns and esmarkii
# While merlengus isn't correlated with anyone.
# 

# Final thing here is getting fishing mortality estimates by stock, we can use these to project under a business as usual case
# Then tweak around these as we see fit.

fm.ts <- do.call("rbind",fm.sub) 
fm.ts$stock <- factor(fm.ts$stock,levels= Stocks)
# Now I wonder if fishing is correlated between some stocks...
fm.cor.ts <- fm.ts %>% tidyr::pivot_wider(names_from = stock,values_from = fm)
colnames(fm.cor.ts) <- c('years',"Melanogrammus aeglefinus","Trisopterus esmarkii",
                          "Pollachius virens", "Gadus morhua", "Scopthalmus maximus",
                          "Pleuronectes platessa 7d", "Pleuronectes platessa 4,20", "Solea solea 7d", "Solea solea 4",
                          "Merlangius merlangus")
fm.cor.mat <- cor(fm.cor.ts[,-1],use="complete.obs")

p.fm.cor <- ggpairs(data = fm.cor.ts[,-1])

save_plot(paste0(loc,"/ICM/Figures/NS_sims/NS_fm_correlation_plot.png"),p.fm.cor,base_width = 22,base_height = 22)



# OK, so I'm going to base past fishing mortality on what was observed in the past
ggplot(fm.ts) + geom_line(aes(x=years,y=fm,group=stock,color=stock))
# So get the mean/variance, on the log scale since we'll need to log-normal distro this
# because there are some 0's in the data we have to take mean the log of that, rather than vice versa.

fish.mort <- fm.ts %>% dplyr::group_by(stock) %>% dplyr::summarise(mn = mean(fm,na.rm=T),
                                                                   sd = sd(fm,na.rm=T))
fish.mort$stock <- Stocks
fish.mort <- fish.mort[,c(2,3,1)] # reorder for below



############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea



# These 4 stocks are the ones we use to get the correlation for each 'functional group' (which I define as the ones that are correlated with each other)
focal.stocks <- c("ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus", 
                  "ICES-WGNSSK_NS 4-6- 3a_Pollachius_virens",
                  "ICES-WGNSSK_NS 7d_Pleuronectes_platessa",
                  "ICES-WGNSSK_NS 4-7d_Merlangius_merlangus")

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



# So everything will need to get wrapped up in a simulation loop
for(j in 1:n.sims)
{
st.time <- Sys.time()

# Now get the fecundity and natural mortality matricies for the simulations
for(i in Stocks)
{
  # This picks the years of the mx and nm samples, note all stocks sample from the same year, tho it is separate for mx and nm
  samp.mx <- round(runif(n.yrs.proj,1,nrow(mx.sub[[i]])))
  samp.nm <- round(runif(n.yrs.proj,1,nrow(mx.sub[[i]])))
  
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
  
  # This is bit hacky, but these if's will get us the inter-stock fecundity correlations
  if(i == "ICES-WGNSSK_NS 4-3aN_Trisopterus_esmarkii")
  {
    row.sel <- which(row.names(mx.cor.mat) == "Trisopterus.esmarkii")
    col.sel <- which(row.names(mx.cor.mat) =="Melanogrammus.aeglefinus")
    corr <- mx.cor.mat[row.sel,col.sel]
    mx.d.tmp <- mx.dev[["ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus"]]-1
    # Toss in a little variability in the correlation, this weakens the correlation a little bit (approx 10%), but does a pretty solid job.
    mx.dev[[i]] <- rnorm(n.yrs.proj,mx.d.tmp,abs(corr))+1
  }# end Functional Group 1 if statement
  
  # More hacky
  if(i %in% c("ICES-WGNSSK_NS 4-7d,20_Gadus_morhua","ICES-WGNSSK_NS4 _Scopthalmus_maximus"))
  {
    if(i == "ICES-WGNSSK_NS 4-7d,20_Gadus_morhua") row.sel <- which(row.names(mx.cor.mat) == "Gadus.morhua")
    if(i == "ICES-WGNSSK_NS4 _Scopthalmus_maximus") row.sel <- which(row.names(mx.cor.mat) == "Scopthalmus.maximus")
    col.sel <- which(row.names(mx.cor.mat) =="Pollachius.virens")
    corr <- mx.cor.mat[row.sel,col.sel]
    mx.d.tmp <- mx.dev[["ICES-WGNSSK_NS 4-6- 3a_Pollachius_virens"]]-1
    # Toss in a little variability in the correlation, this weakens the correlation a little bit (approx 10%), but does a pretty solid job.
    mx.dev[[i]] <- rnorm(n.yrs.proj,mx.d.tmp,abs(corr))+1
  } # end Functional Group 2 if statement
  
  # Final hacky section of this hacky section
  if(i %in% c("ICES-WGNSSK_NS 4,20_Pleuronectes_platessa","ICES-WGNSSK_NS 7d._Solea_solea","ICES-WGNSSK_NS4 _Solea_solea"))
  {
    if(i == "ICES-WGNSSK_NS 4,20_Pleuronectes_platessa") row.sel <- which(row.names(mx.cor.mat) == "Pleuronectes.platessa.4.20")
    if(i == "ICES-WGNSSK_NS 7d._Solea_solea") row.sel <- which(row.names(mx.cor.mat) == "Solea.solea.7d")
    if(i == "ICES-WGNSSK_NS4 _Solea_solea") row.sel <- which(row.names(mx.cor.mat) == "Solea.solea.4")
    col.sel <- which(row.names(mx.cor.mat) =="Pleuronectes.platessa.7d")
    corr <- mx.cor.mat[row.sel,col.sel]
    mx.d.tmp <- mx.dev[["ICES-WGNSSK_NS 7d_Pleuronectes_platessa"]]-1
    # Toss in a little variability in the correlation, this weakens the correlation a little bit (approx 10%), but does a pretty solid job.
    mx.dev[[i]] <- rnorm(n.yrs.proj,mx.d.tmp,abs(corr))+1
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
  # if(c.effect < n.yrs.proj)
  # {
  #   ce.mx <-  (climate.mx.effect/10)*(1:(n.yrs.proj-c.effect))
  #   tmp.mx[[i]][c.effect:n.yrs.proj,] <- tmp.mx[[i]][c.effect:n.yrs.proj,] - ce.mx*tmp.mx[[i]][c.effect:n.yrs.proj,]
  #   # Need to be careful here since we are on a proportional scale, lets go to instantaneous and convert back
  #   ce.nm <-  (climate.nm.effect/10)*(1:(n.yrs.proj-c.effect)) # So we'll say this is instantaneous rate.
  #   tmp.insta.nm <- -log(1-tmp.nm[[i]])
  #   tmp.nm[[i]][c.effect:n.yrs.proj,] <- 1-exp(-(tmp.insta.nm[c.effect:n.yrs.proj,] + ce.nm*tmp.nm[[i]][c.effect:n.yrs.proj,]))
  # } # end if(c.effect < n.yrs.proj)
}  # end for(i in Stocks)

# OK, so I have everything I need now to run the simulations forward.

res.ts <- NULL
res.r <- NULL
for(s in Stocks)
{
  # The last year of data for all the stocks is 2016
  years <- 2016:(2017+n.yrs.proj-2)
  inst.nat.mort <- tmp.nm[[s]] # 
  age.mat <- tmp.mat[[s]]
  mx <- tmp.mx[[s]] 
  vpa.abund <- vpa.sub[[s]]
  weight.age <- tmp.waa[[s]]
  ages <- tmp.age[[s]]
  N.start <- vpa.abund[length(vpa.abund)]
  K <- max(vpa.abund,na.rm=T)
  fm <- unlist(fish.mort[fish.mort$stock == s,1:2]) # Mean and variance, is log-normally distributed.
  
  tst <- for.sim(years,
                 mat.age = age.mat,
                 nm = inst.nat.mort,
                 w.age = weight.age,
                 ages = ages,
                 rems =  list("R_based",0), #fm,
                 fecund = mx,
                 N.start = N.start,
                 pop.model = 'bounded_exp', K = rep(K,length(years)),
                 sim= "project",
                 n.sims = 1,
                 sd.mat = 0,
                 sd.nm = 0,
                 sd.wt = 0,
                 sd.fecund = 0
                 )
  #ggplot(tst$Pop) + geom_line(aes(x=years,y=abund)) + geom_hline(yintercept = K)
  for.proj(option = "logistic",pop.last = 10,K=20,r=0.6)

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

quants <- ts.final %>%  dplyr::group_by(years,stock) %>% dplyr::summarise(L.50 = quantile(abund,probs=c(0.25),na.rm=T),
                                                                          med = median(abund,na.rm=T),
                                                                          U.50 = quantile(abund,probs=c(0.75),na.rm=T),
                                                                          fml.50 = quantile(fm,probs=c(0.25),na.rm=T),
                                                                          fm = median(fm,na.rm=T),
                                                                          fmu.50 = quantile(fm,probs=c(0.75),na.rm=T))
# If happy save the 3 objects
# saveRDS(object = ts.final,file = paste0("D:/Github/ICM/Results/NS_climate_starts_at_year_",c.effect,
#                                           "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,
#                                          "_time_series_projections.Rds"))
# 
# saveRDS(object = quants,file = paste0("D:/Github/ICM/Results/NS_climate_starts_at_year_",c.effect,
#                                       "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,
#                                       "_time_series_quantiles.Rds"))
# 
# saveRDS(object = r.final,file = paste0("D:/Github/ICM/Results/NS_climate_starts_at_year_",c.effect,
#                                        "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,
#                                        "_r_projections.Rds"))


# Two simple plots.
p.sims <- ggplot(ts.final) + geom_line(aes(x=years,y=abund,group = sim),alpha=1) +
                             facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) + scale_x_continuous(breaks = seq(2010,2200,by=10)) 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_all_realizations_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims,base_height = 12,base_width = 20)

p.sims.quants <- ggplot(quants) + geom_line(aes(x=years,y=med)) + facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) +
                                  geom_ribbon(data=quants, aes(x=years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_quantiles_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims.quants,base_height = 12,base_width = 20)












### Section in progress, getting the K's and the bimoasses ##

load(file = "C:/Users/Owner/Documents/Github/ICM/Results/model_inputs.Rdata")


Stocks <- Stocks[grep("WGNSSK",Stocks)]
# Haddock in NS is busted before 1972, not sure why, but I'm chucking all that data here in a very sloppy way

# So lets look at total abundance and total biomass in the system by year...

years.ns <- NULL
vpa.ns <- NULL
bm.ns <- NULL
num.ns <- NULL
waa.ns <- NULL
pnm.ns <- NULL
rem.ns <- NULL
mx.ns <- NULL
am.ns <- NULL
for(i in  Stocks)
{
years.ns[[i]] <- years.tmp[[i]]
vpa.ns[[i]] <- vpa.tmp[[i]]
num.ns[[i]] <- ASR_long %>% dplyr::filter(Stock == i,type == "Num")
num.ns[[i]] <- num.ns[[i]] %>% dplyr::filter(age != "tot")
waa.ns[[i]] <- ASR_long %>% dplyr::filter(Stock == i,type == "WA")
bm.ns[[i]] <- data.frame(Year = num.ns[[i]]$Year,Stock = num.ns[[i]]$Stock,age = num.ns[[i]]$age,
                         bm = num.ns[[i]]$value*waa.ns[[i]]$value,
                         num = num.ns[[i]]$value)
pnm.ns[[i]] <- pnm.tmp[[i]]
rem.ns[[i]] <- rem.tmp[[i]]
mx.ns[[i]] <- mx.tmp[[i]]
am.ns[[i]] <- am.tmp[[i]]
}

ns.had <- "ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus"
years.ns[[ns.had]] <- years.ns[[ns.had]][-1:-7]
pnm.ns[[ns.had]] <- pnm.ns[[ns.had]][-1:-7,]
waa.ns[[ns.had]] <- waa.ns[[ns.had]][-1:-7,]
#ages.ns[[ns.had]] <- ages.ns[[ns.had]][-1:-7]
rem.ns[[ns.had]] <- rem.ns[[ns.had]][-1:-7,]
mx.ns[[ns.had]] <- mx.ns[[ns.had]][-1:-7,]
bm.ns[[ns.had]] <- bm.ns[[ns.had]] %>% filter(Year >= 1972)
vpa.ns[[ns.had]] <- vpa.ns[[ns.had]][-1:-7]
am.ns[[ns.had]] <- am.ns[[ns.had]][-1:-7,]

ns.sole <- "ICES-WGNSSK_NS4 _Solea_solea"
years.ns[[ns.sole]] <- years.ns[[ns.sole]][-1:-7]
pnm.ns[[ns.sole]] <- pnm.ns[[ns.sole]][-1:-7,]
waa.ns[[ns.sole]] <- waa.ns[[ns.sole]][-1:-7,]
#ages.ns[[ns.sole]] <- ages.ns[[ns.sole]][-1:-7]
rem.ns[[ns.sole]] <- rem.ns[[ns.sole]][-1:-7,]
mx.ns[[ns.sole]] <- mx.ns[[ns.sole]][-1:-7,]
bm.ns[[ns.sole]] <- bm.ns[[ns.sole]] %>% filter(Year >= 1971)
vpa.ns[[ns.sole]] <- vpa.ns[[ns.sole]][-1:-7]
am.ns[[ns.sole]] <- am.ns[[ns.sole]][-1:-7,]

bm.tst <- do.call("rbind",bm.ns)

bm.tot <- bm.tst %>% dplyr::group_by(Stock,Year) %>% dplyr::summarise(bm = sum(bm,na.rm=T),
                                                                      num = sum(num,na.rm=T))

eco.bm <- bm.tot %>% dplyr::group_by(Year) %>% dplyr::summarise(num = sum(num),bm = sum(bm))

# Now all stocks start in 1990 so need to consider that...
bm.final <- left_join(bm.tot,eco.bm,by=c("Year"))
names(bm.final) <- c("Stock","Year","bm.stock","num.stock","num.total","bm.total")
# Get proportions...
bm.final <- bm.final %>% dplyr::mutate(bm.prop = bm.stock/bm.total,
                                          num.prop = num.stock/num.total)
bm.final <- bm.final[bm.final$bm.stock > 0,]

what.year <- bm.final %>% dplyr::group_by(Stock) %>% dplyr::summarize(min = min(Year),
                                                                      max = max(Year))
first.year <- max(what.year$min)
last.year <- min(what.year$max)

bm.best <- bm.final %>% dplyr::filter(Year %in% first.year:last.year)
# So this is our target biomass.  Now we want to divy



ggplot(bm.best) + geom_line(aes(x=Year,y=bm.total),linewidth=2) + ylim(c(0,max(bm.final$bm.total))) + geom_hline(yintercept = mean(bm.best$bm.total))
ggplot(bm.best) + geom_line(aes(x=Year,y=num.total)) + ylim(c(0,max(bm.final$num.total))) 
# Autocorrelation in K.
pacf(bm.best$bm.total)

ggplot(bm.best) + geom_line(aes(x=Year,y=bm.prop,group = Stock,color=Stock),linewidth=2) #+ ylim(c(0,0.5))

ggplot(bm.best) + geom_bar(aes(x=Year,y=bm.prop,fill=Stock),position="stack",stat = 'identity') #+ ylim(c(0,0.5))

prop.acfs <- NULL
for(i in Stocks) prop.acfs[[i]] <- pacf(bm.best$bm.prop[bm.best$Stock == i])

catch <- ASR_long %>% dplyr::filter(Stock == "ICES-AFWG_NEA1-2_Gadus_morhua", type == "Catch")             
an.num <- num %>% dplyr::filter(age != 'tot') %>% dplyr::group_by(Year) %>% dplyr::summarise(nums = sum(value,na.rm=T))
num <- ASR_long %>% dplyr::filter(Stock == "ICES-AFWG_NEA1-2_Gadus_morhua", type == "Num")
