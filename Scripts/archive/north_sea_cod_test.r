# OK, so using the ICES assessments here's what we get for North Sea cod.
library(readxl)
library(tidyverse)
library(rio)
library(ggthemes)

# OK, lets give this a spin...
funs <- c("https://raw.githubusercontent.com/Dave-Keith/ICM/main/Scripts/functions/backwards_sim.r")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}

# Bring in the data
# Handy import function from the rio package... tricky part was getting the right filename here...
age.mat <- rio::import('https://github.com/Dave-Keith/ICM/blob/main/Data/Cod_dat.xlsx?raw=true',which = 'age_mat')
n.offspring <- rio::import('https://github.com/Dave-Keith/ICM/blob/main/Data/Cod_dat.xlsx?raw=true',which = "age_fecundity")
nat.mort <-  rio::import('https://github.com/Dave-Keith/ICM/blob/main/Data/Cod_dat.xlsx?raw=true',which = "age_nat_mort")
abund <- rio::import('https://github.com/Dave-Keith/ICM/blob/main/Data/Cod_dat.xlsx?raw=true',which = "abundance") # This is in 1000s
weight.age <- rio::import('https://github.com/Dave-Keith/ICM/blob/main/Data/Cod_dat.xlsx?raw=true',which = "age_weight") 
removals <- rio::import('https://github.com/Dave-Keith/ICM/blob/main/Data/Cod_dat.xlsx?raw=true',which = "age_removals")

# Or if you want to bring in the local data.
#age.mat <- read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "age_mat")
#n.offspring <- read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "age_fecundity")
##nat.mort <-  read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "age_nat_mort")
#abund <- read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "abundance") # This is in 1000s, I need to correct for that!
#weight.age <- read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "age_weight") 
#removals <- read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "age_removals")

# OK, so rejig to get age 0, use a super high M for that and then get recruits from that...
nat.mort[nrow(nat.mort)+1,] <- cbind(year = 2021,nat.mort[nrow(nat.mort),-1])

nm.tmp <- data.frame(year = nat.mort$year,`0` = 2)
nat.mort <- cbind(nm.tmp,nat.mort[,2:8])
# Using that M let's get a new abundance estimate...
# First turn M into a proportion
prop.nat.mort <- nat.mort
prop.nat.mort[,c(-1,-ncol(nat.mort))] <- 1-exp(-prop.nat.mort[,c(-1,-ncol(nat.mort))])

n.tmp <- data.frame(year = abund$Year,`0` = abund[,2]/(1-prop.nat.mort[,2]))
abund <- cbind(n.tmp,abund[,2:8])
# Need to add an age at maturity column...
am.tmp <- data.frame(year = age.mat$Year,`0` = 0)
age.mat <- cbind(am.tmp,age.mat[,2:8])
# And weight at age
wa.tmp <- data.frame(year = weight.age$Year,`0` = 0.01)
weight.age <- cbind(wa.tmp,weight.age[,2:8])


# What if we up the ante on the older age classes...
#prop.nat.mort[,3:7] <- 0.95

rem <- c(removals$total_catch,NA)
years <- age.mat$year
N.end <- rowSums(abund[nrow(abund),2:8])
vpa.abund <- rowSums(abund[,3:8]) # 
#dec.mod <- lm(log(vpa.abund)~years)
#summ.dm <- summary(dec.mod)
#dec.rate <- summ.dm$coefficients[2,1]

# The real mx matrix, recruits produced per individual in each age class... Not perfect as I need to offset recruits/ssb, but close enough for the moment..
recruits <- abund[,2]
# This is the number of spawners
ssn <- abund[,2:8]*age.mat[,c(-1,-ncol(age.mat))]
# To get biomass of spawners
ssb <- ssn*weight.age[,c(-1,-ncol(weight.age))]
# Total biomass by year in number/kg or 1000's per tonne
tot.ssb <- rowSums(ssb)
# The number of recruits per kilogram of SSB
r.p.ssb <- recruits/tot.ssb
# Number of recruits produced by everyone in each age class
recs.per.age <- cbind(r.p.ssb*ssb[,1],r.p.ssb*ssb[,2],r.p.ssb*ssb[,3],r.p.ssb*ssb[,4],r.p.ssb*ssb[,5],r.p.ssb*ssb[,6],r.p.ssb*ssb[,7])
# SO now we want the number of recruits produced by each mom.
mx <- recs.per.age/ssn/2 # Moms only! Dividing by around 8 really nails it for NS cod for whatever reason :-)
# Something is wrong with the decline rate method, but the exponential and logistic are working pretty... pretty... pretty good...


tst <- back.sim(years = years,mat.age = age.mat[,c(-1,-ncol(age.mat))],nm = prop.nat.mort[,c(-1,-ncol(prop.nat.mort))],
               w.age = weight.age[,c(-1,-ncol(weight.age))],ages = 0:6,rems = rem,
               fecund = mx,N.end =N.end,pop.model = 'exponential',
               n.sims = 5,sd.mat = 0.5,sd.nm = 0.5,sd.wt = 0.5,sd.fecund = 0.5)


#Combine the data 
did.it.work <- data.frame(abund = c(vpa.abund,tst$Pop$abund),years = c(years,tst$Pop$years),sim = c(rep('VPA',length(years)),tst$Pop$sim))
# Here get the Upper and lower 50% quantiles to make a functional boxplot
quants <- did.it.work %>% dplyr::group_by(years) %>% dplyr::summarise(L.50 = quantile(abund,probs=c(0.25)),
                                                                      med = median(abund),
                                                                      U.50 = quantile(abund,probs=c(0.75)))

# All the results in a plot
ggplot(did.it.work) + geom_line(aes(x=years,y=abund,color=sim)) +xlab("") + ylab("Abundance (1000s)") + 
                      geom_line(data=did.it.work %>% dplyr::filter(sim == "VPA"),aes(x=years,y=abund),color='black',size=2) +
                      scale_y_continuous(breaks = seq(0,3e6,by=5e5)) + scale_x_continuous(breaks = seq(1960,2025,by=5)) +  
                      theme_few() + theme(legend.position = 'none') + scale_color_viridis_d(end = 0.75)
# Same thing but functional boxplots
ggplot(quants) + geom_line(aes(x=years,y= med)) + geom_ribbon(aes(x=years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') +
                 geom_line(data=did.it.work %>% dplyr::filter(sim == "VPA"),aes(x=years,y=abund),color='black',size=2) +
                 xlab("") + ylab("Abundance (1000s)") + 
                 scale_y_continuous(breaks = seq(0,3e6,by=5e5)) + scale_x_continuous(breaks = seq(1960,2025,by=5)) +  
                 theme_few() + theme(legend.position = 'none')
# Plotting how the estimate of r changes over time (which you can do when we have time varying inputs (e.g. Wgt, M, or maturity at age)
ggplot(tst$r) + geom_line(aes(x=years,y=r,group=n.sims,color=n.sims)) + 
                xlab("") + ylab("Estimated growth rate (r)") + 
                scale_y_continuous(breaks = seq(0,1.5,by=0.1)) + scale_x_continuous(breaks = seq(1960,2025,by=5)) +  
                theme_few() + theme(legend.position = 'none') + scale_color_viridis_c(end = 0.75)
