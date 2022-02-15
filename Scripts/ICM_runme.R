##### Incidental Catch Model

### H.Bowlby April 2020
# Base version: incorporates code from POR.sim.R - NOTE that the u.calc and f.crit.calc functions in the original code were modified
## and MC simulation_matrix_POR_v1.R

### D. Keith Feb 2022
# DK revising to generalize the code.  Extracted functions so loading is automatic and each function can be revised independently as needed
# Goal here is to have this work with any stock and then we can start seeing how well this model works for 'everything'





setwd('c:/mydocs/ICCAT/POR2020')

## required libraries for Enric's code
library(MASS)
library(popbio)
library(tidyverse)
library(grid)
library(truncnorm)

library(ggplot2)


#################################################################

## Load functions 

# Bring in our in house functions. First combine them all in a vector
funs <- c("https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/plot_por_sim.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/por_sim.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/parms_calc.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/u_calc.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/F_crit.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/Lotka_r.r",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Maps/standard_SHK.r")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}



#################################################################
### note that the parms.calc function generates the information required for the main simulation model
## this function is simplified from standard.SHK - runs faster than calling the full standard.SHK function

### the standard.SHK function called by the plotting code is specifically used to do calculations
## relative to reference points. By separating the main simulations from the evaluation components
## we can run one and not the other if need be - speeds up overall evaluation.
## the functions could be better integrated with more time.






  

# Standard Leslie matrix models with MC simulations------------------------------------------------------------

# Runs matrix analysis with a post-breeding census
## n is number of MC simulations
## Age at first capture (AAFC) and Fishing mortality (F) can be specified
## repro.cycle refers to periodicity: 1 = annual, 1.5 = annual + biennial, 2 = biennial











########################################################

### start analysis

########################################################
## to actually run the analysis, load in all of the functions, including por.sim() below.
## then call: por.sim.result<-por.sim(n.sims=??, repro.cycle=1 or 2)
## then call: plot.por.sim()

### note that the number of simulations, the age at first selectivity to the fishery and 
# the removals scenarios that shoudl be evaluated need to be set in the por.sim() 
# and standard.SHK functions

### NOTE that I left the number of iterations hardcoded into standard.SHK so that the main 
## simulation model coudl be run with a different number of iterations than the life hsitory evaluation.






## read in removals data from Task 1 - Provided by Carlos JUNE 18 - incorporates modifications to address
## changes to catch series from the meeting - estimation of recent catches .
# see associated details - t1nc_POR-all_20200605_v3_FINALcatchSERIES.xlsx
calc.N<-read.csv('calc.N.NWupdated.csv')
## note - calculation of N in the data file still uses mean size from Canada 39.34

## note the assumption is a median length of 149 cm and median weight of 39.34kg to gets numbers from biomass
## from 147 and W_a and W_b above, get 39.35 rather than 
nafo<-read.csv('calc.N.NAFO.csv')
#### multiple options for removals series.

## option 1: using size frequency all years combined. 
### add size distribution from Rui - Provided June 3, 2020
size<-read.csv('NW.size.csv')
### use growth parameters above to transform this directly into biomass.
W_a <- 1.482E-5
W_b <- 2.9641

size$wt<-W_a*size$Size_cm_FL^W_b

junk<-hist(size$wt,breaks=43,plot=F)  ### every 5cm size class
##proportion of animals in each weight bin
split<-junk$counts/sum(junk$counts)

## biomass * split / weight in each bin
N.ratio<-c()
for (i in 1:length(calc.N$por.kg))
{
  ## weighted average in 5 cm length bins
  N.ratio[i]<-sum(calc.N$por.kg[i]*split/junk$breaks[2:length(junk$breaks)])
}


## option 2: using median from NW size frequency, all years combined
N.med<-calc.N$por.kg/median(size$wt)

## year-specific length-frequency or median size.

N.med.yr<-c()
N.ratio.yr<-c()
year<-calc.N$year

for (i in 1:length(year))
{
  if(i < 35)
  {
    junk<-hist(size$wt[size$Year==1994],breaks=43,plot=F)  ### every 5cm size class
    ##proportion of animals in each weight bin
    split<-junk$counts/sum(junk$counts)
       ## biomass * split / weight in each bin
       N.ratio.yr[i]<-sum(calc.N$por.kg[i]*split/junk$breaks[2:length(junk$breaks)])
       N.med.yr[i]<-calc.N$por.kg[i]/median(size$wt[size$Year==1994])
  }
  else
  {
    junk<-hist(size$wt[size$Year==year[i]],breaks=43,plot=F)  ### every 5cm size class
    ##proportion of animals in each weight bin
    split<-junk$counts/sum(junk$counts)
    ## biomass * split / weight in each bin
    N.ratio.yr[i]<-sum(calc.N$por.kg[i]*split/junk$breaks[2:length(junk$breaks)])
    N.med.yr[i]<-calc.N$por.kg[i]/median(size$wt[size$Year==year[i]])
  }
}

N.options<-data.frame(calc.N,N.ratio,N.ratio.yr,N.med,N.med.yr)

write.csv(N.options,'options.for.biomass.to.number.animals.NW.csv')

# average removals for the most recent 5 years (2014-2018) are 39 mt = 39000kg. 
# mean size of removals = 147cm for Canada; used 149 because listed in Figure 3 of Campana et al. 2012
# 149cm = average weight of 39.34 kg
## gives an average removals of 991 animals (39000/39.34) in 2014-2018: used below to represent recent time period in forward projections
# HOWEVER, average removals thorughout the time series are 1789mt; gives 45,475 animals (using Task 1 series)



  
  