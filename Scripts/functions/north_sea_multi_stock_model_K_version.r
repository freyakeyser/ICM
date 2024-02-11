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
load(file = "d:/Github/ICM/Results/model_inputs.Rdata")
#load(file = "D:/Github/ICM/Results/model_inputs.Rdata")
loc <- 'D:/GitHub/ICM'

load(file = paste0(loc,"/Results/all_cleaned_forward_tune_summaries_fec_nm.Rdata"))

########################### End Section 1 Loading  ###################################### End Section 1 Loading ###############################################






########################## Section 2 Parameters ########################## Section 2 Parameters ########################## Section 2 Parameters

yrs.all <- 1990:2016 # These are the years we have data for all 10 stocks
n.yrs.proj <- 25
n.sims <- 5


# OK, so I'm going to base past fishing mortality on what was observed in the past
#ggplot(fm.ts) + geom_line(aes(x=years,y=fm,group=stock,color=stock))
# So get the mean/variance, on the log scale since we'll need to log-normal distro this
# because there are some 0's in the data we have to take mean the log of that, rather than vice versa.

# fish.mort$stock <- Stocks
# fish.mort <- fish.mort[,c(2,3,1)] # reorder for below

# SO now I want to get the carrying capacity and slice it up by stock

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
ages.ns <- NULL
for(i in  Stocks)
{
  years.ns[[i]] <- years.tmp[[i]]
  vpa.ns[[i]] <- vpa.tmp[[i]]
  ages.ns[[i]] <- ages.tmp[[i]]
  num.ns[[i]] <- ASR_long |> collapse::fsubset(Stock == i & type == "Num")
  num.ns[[i]] <- num.ns[[i]] |> collapse::fsubset(age != "tot")
  waa.ns[[i]] <- ASR_long|> collapse::fsubset(Stock == i & type == "WA")
  bm.ns[[i]] <- data.frame(Year = num.ns[[i]]$Year,Stock = num.ns[[i]]$Stock,age = num.ns[[i]]$age,
                           bm = num.ns[[i]]$value*waa.ns[[i]]$value,
                           num = num.ns[[i]]$value)
  pnm.ns[[i]] <- pnm.tmp[[i]]
  rem.ns[[i]] <- rem.tmp[[i]]
  rem.ns[[i]]$Stock <- i
  mx.ns[[i]]  <- for.tune.all[[i]]$fecund.opt
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

rem.tst <- do.call("rbind",rem.ns)
# This isn't perfect, it's the fishing removals in numbers not biomass, but that is 
# what we are modelling... should make biomass...
fm.dat <- left_join(bm.tot,rem.tst,by=c("Stock",'Year'))
# This seems ok, we want this as an exploitation rate
fm.dat$exploit <- fm.dat$rem/fm.dat$num
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

# # Autocorrelation in K.
K.cor <- pacf(log(bm.best$num.total))
# 

# OK so now make a fake K time series, note the stock is irrlevant here as we're looking at total bm...
N.target.sum <- bm.best %>% dplyr::filter(Stock == "ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus") %>% 
                           dplyr::summarise(sd = sd(log(num.total)),
                                            mn = mean(log(num.total)),
                                            med = median(log(num.total)))

K.devs <- NULL
cors <- NULL
K.sims <- NULL
for(i in 1:n.sims) 
{
   K.devs[[i]] <- as.data.frame(arima.sim(list(order = c(1,0,0), ar = K.cor$acf[1]), 
                                                                 n = n.yrs.proj,
                                                                 sd = N.target.sum$sd))
  #K.devs[[i]] <- data.frame(x = rep(0,n.yrs.proj)) # Lets see what happens when it's fixed at a mean level.
  K.sims[[i]] <- exp(N.target.sum$mn + K.devs[[i]]$x)
}

# For the moment, lets just cut the world up so each stock gets a fixed % of K bassed on their biomass
prop.stock <- bm.best %>% dplyr::group_by(Stock) %>% dplyr::summarise(prop = mean(num.prop))
# Now get this for each stock...
K.stock <- NULL
fm.stock <- NULL
tmp <- NULL
  for(i in 1:n.sims)
  {
    tmp <- NULL
    for(s in Stocks)
    {
      prop <- prop.stock %>% dplyr::filter(Stock == s)
      tmp[[s]] <- prop$prop *as.numeric(K.sims[[i]])
      if(i == 1) fm.stock[[s]] <- fm.dat %>% dplyr::filter(Stock ==s) %>% dplyr::summarise(mn = median(exploit,na.rm=T),
                                                                                           sd = sd(log(exploit[exploit > 0]),na.rm=T))
    }
    K.stock[[i]] <- tmp
    
  }

# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.total),linewidth=2) + ylim(c(0,max(bm.final$bm.total))) + geom_hline(yintercept = mean(bm.best$bm.total))
# ggplot(bm.best) + geom_line(aes(x=Year,y=num.total)) + ylim(c(0,max(bm.final$num.total))) 

# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.prop,group = Stock,color=Stock),linewidth=2) #+ ylim(c(0,0.5))
# 
# ggplot(bm.best) + geom_bar(aes(x=Year,y=bm.prop,fill=Stock),position="stack",stat = 'identity') #+ ylim(c(0,0.5))
# 
# prop.acfs <- NULL
# for(i in Stocks) prop.acfs[[i]] <- pacf(bm.best$bm.prop[bm.best$Stock == i])




############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea



# These 4 stocks are the ones we use to get the correlation for each 'functional group' (which I define as the ones that are correlated with each other)
tmp.mx <- NULL
tmp.nm <- NULL
tmp.mat <- NULL
tmp.age <- NULL
tmp.waa <- NULL
mx.dev <- NULL
nm.dev <- NULL
ts.unpack <- NULL
r.unpack <- NULL
Ks <- NULL

# So everything will need to get wrapped up in a simulation loop
for(j in 1:n.sims)
{
  st.time <- Sys.time()
  # This is getting the K for the stock for a particular simulation
  Kss <- K.stock[[j]]
  # Now get the fecundity and natural mortality matricies for the simulations

  res.ts <- NULL
  res.r <- NULL
  for(s in Stocks)
  {
    # The last year of data for all the stocks is 2016
    years <- 2016:(2017+n.yrs.proj-2)
    prop.nat.mort <- pnm.ns[[s]] 
    age.mat <- am.ns[[s]]
    mx <- mx.ns[[s]] 
    #vpa.abund <- vpa.ns[[s]]
    #weight.age <- waa.ns[[s]]
    ages <- ages.ns[[s]]
    Ks <- Kss[[s]]
    vpa.ns  <- bm.best$num.stock[bm.best$Stock == s]
    N.start <- vpa.ns[length(vpa.ns)]
    
    tst <- for.sim(years,
                   mat.age = age.mat,
                   nm = -(log(1-prop.nat.mort)),
                   w.age = NULL,
                   ages = ages,
                   rems =  list(fm.stock[[s]]$mn,fm.stock[[s]]$sd), #fm,
                   fecund = mx,
                   N.start = N.start,
                   pop.model = 'bounded_exp', 
                   sim= "project",
                   n.sims = 10,
                   sd.mat = 0,
                   sd.nm = 0,
                   sd.wt = 0,
                   sd.fecund = 0,
                   K = 5*Ks)
    
#ggplot(tst$Pop) + geom_line(aes(x=years,y=abund,color=sim,group=sim)) 


res.ts[[s]] <- data.frame(tst$Pop[,-2],stock = s,sim= j)
res.r[[s]] <- data.frame(tst$r[,-3],stock=s,sim=j)

  } # end stock look 
  
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
p.sims <- ggplot(ts.final) + geom_line(aes(x=years,y=abund,group = sim),alpha=0.8) +
  facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) + scale_x_continuous(breaks = seq(2010,2200,by=10)) 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_all_realizations_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims,base_height = 12,base_width = 20)

p.sims.quants <- ggplot(quants) + geom_line(aes(x=years,y=med)) + facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) +
  geom_ribbon(data=quants, aes(x=years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_quantiles_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims.quants,base_height = 12,base_width = 20)



### Section in progress, getting the K's and the bimoasses ##

catch <- ASR_long %>% dplyr::filter(Stock == "ICES-AFWG_NEA1-2_Gadus_morhua", type == "Catch")             
an.num <- num %>% dplyr::filter(age != 'tot') %>% dplyr::group_by(Year) %>% dplyr::summarise(nums = sum(value,na.rm=T))
num <- ASR_long %>% dplyr::filter(Stock == "ICES-AFWG_NEA1-2_Gadus_morhua", type == "Num")
