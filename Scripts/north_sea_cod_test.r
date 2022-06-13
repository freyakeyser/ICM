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

# Tidy up the data for input...
prop.nat.mort <- nat.mort
prop.nat.mort[,c(-1,-ncol(nat.mort))] <- 1-exp(-prop.nat.mort[,c(-1,-ncol(nat.mort))])

rem <- c(removals$total_catch,NA)
years <- age.mat$Year
N.end <- rowSums(abund[59,2:7])
vpa.abund <- rowSums(abund[,2:7])
#dec.mod <- lm(log(vpa.abund)~years)
#summ.dm <- summary(dec.mod)
#dec.rate <- summ.dm$coefficients[2,1]

# The real mx matrix, recruits produced per individual in each age class... Not perfect as I need to offset recruits/ssb, but close enough for the moment..
recruits <- abund[,2]
ssn <- abund[,2:7]*age.mat[,c(-1,-ncol(age.mat))]
ssb <- ssn*weight.age[,c(-1,-ncol(weight.age))]
tot.ssb <- rowSums(ssb)
r.p.ssb <- recruits/tot.ssb
recs.per.age <- cbind(r.p.ssb*ssb[,1],r.p.ssb*ssb[,2],r.p.ssb*ssb[,3],r.p.ssb*ssb[,4],r.p.ssb*ssb[,5],r.p.ssb*ssb[,6])
mx <- recs.per.age/ssn/2 # Moms only! Dividing by around 8 really nails it for NS cod for whatever reason :-)
# Something is wrong with the decline rate method, but the exponential and logistic are working pretty... pretty... pretty good...


tst <- icm.sim(years = years,mat.age = age.mat[,c(-1,-ncol(age.mat))],nm = prop.nat.mort[,c(-1,-ncol(nat.mort))],
               w.age = weight.age[,c(-1,-ncol(weight.age))],ages = 1:6,rems = rem,
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
