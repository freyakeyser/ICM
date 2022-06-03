# So here I'm gonna take a stab at doing the backwards projections for North Sea Cod using what I can find from 
# the ICES assessments and Fish Base if necessary.
library(readxl)
# Bring in our in house functions. First combine them all in a vector
funs <- c("https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/plot_por_sim.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/main/Scripts/functions/por_sim.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/main/Scripts/functions/parms_calc.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/main/Scripts/functions/u_calc.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/main/Scripts/functions/F_crit.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/main/Scripts/functions/Lotka_r.r",
          "https://raw.githubusercontent.com/Dave-Keith/ICM/main/Scripts/functions/standard_SHK.r")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs)
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}

# Decided to slightly revise the Euler Lotka formulation, this is slightly easier to implement than the old code
# I checked and it gives the same answer :-)
eulerlotka <- function(r) (sum(lx * mx * exp(-r * x)) - 1)^2




source("D:/Github/ICM/Scripts/functions/Lotka_r.r")

# Bring in the data 
age.mat <- read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "age_mat")
n.offspring <- read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "age_fecundity")
nat.mort <-  read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "age_nat_mort")
abund <- read_xlsx("D:/Github/ICM/Data/Cod_dat.xlsx",sheet = "abundance") # This is in 1000s, I need to correct for that!
prop.nat.mort <- nat.mort
prop.nat.mort[,c(-1,-ncol(nat.mort))] <- 1-exp(-prop.nat.mort[,c(-1,-ncol(nat.mort))])

mean.age.mat <- colMeans(age.mat[,c(-1,-ncol(age.mat))])
mean.nm <- colMeans(prop.nat.mort[,c(-1,-ncol(nat.mort))])
mean.noff <- colMeans(n.offspring[,c(-1,-ncol(n.offspring))])

# So here we have the problem of just using fecundity for something highly fecund but high mortality
# Need to get this scaled to recruits...
# 
# tst <- lotka.r(age.mat = mean.age.mat,n.offspring = mean.noff,nat.mort = mean.nm,max.age=6,u=0,sel =1)
# tst$par

# So the easy way to scale n.off to recruits, is to use the numbers of recruits and number of spawners
# but that is using model output, which is probably fine as an 'option' for proof of concept
# but we need something far more generalizable for stocks without a model...

############################# r via method 1, using model output #############################

# So method 1 is getting r using model output....
# The first age class we have becomes our recruits
recruits <- as.data.frame(abund[,2])*1000
# Need to offset the spawner numbers by this many years
lag <- as.numeric(names(recruits))
# The number of female spawners in each as class is simply the age.mat matrix times the numbers matrix
# Not I'm assuming a 50:50 sex distribution here
n.spawners <- 1000*age.mat[,c(-1,-ncol(age.mat))] * abund[,c(-1,-ncol(age.mat))] /2

recruits <- as.vector(recruits)
# Now we can make the fecudity data relative to get the proportion of the total number of eggs that are produced by a particular age class
# and we will assume that every egg is created equal...
# First we get the total number of eggs produced in each age class in each year. Multiply number of offspring each individual in an age class produces by
# the number of individuals in that age class to get the egg production each year.
# Wright had the interesting idea of tossing in the survival of eggs to age 1, which you could put into the 
# Euler-Lotka equation, but I can't think what the advantage of doing that is since mx is 0, but let's test it out...
n.eggs <- n.offspring[,c(-1,-ncol(n.offspring))]*n.spawners/2
# Now we get the proportion of the total egg production attributed to each age class
prop.fec <- n.eggs/apply(n.eggs,1,sum)

# We then multiply this by the offset recruit numbers and we the number of recruits produced by each age class
rec.produced.by.age <- c(recruits[-(1:lag),],rep(NA,lag))*prop.fec
# Then we divide that by the number if spawners in each age class and we have a fecundity matrix.
# Did I just make the math way harder than it needed to be?
fec.mat <- rec.produced.by.age / n.spawners
mean.fec.mat <- colMeans(fec.mat,na.rm=T)

# So we can get 'age 0' survivorship here by taking ratio of recruits to total number of eggs in the previous year...
tot.eggs <- rowSums(n.eggs)

age.0.mort <- 1-c(recruits[-(1:lag),],rep(NA,lag))/tot.eggs
mn.age.0.nm <- mean(age.0.mort,na.rm=T)

# Let's just do this myself...

# here's the estiamte if we include year 0 and have a mortality from eggs to recruits, this is the highest estimate we have
lx <- c(mn.age.0.nm,mean.nm)
#lx <- rep(0.2,6)
si <- 1-lx
max.age <- 6

x<-0:(max.age)
lx<-1
for(i in 2:(length(si))) lx[i]<-lx[i-1]*si[i-1]


mx <- c(0,mean.noff)
#mx <- mean.noff
lxmx<-lx*mx
x <- 0:6
# Now we are cooking!
junk<-nlminb(start = 0.1, obj = eulerlotka)
junk$par


# Or we can ignore the eggs....

lx <- mean.nm
#lx <- rep(0.2,6)
si <- 1-lx
max.age <- 6

x<-1:(max.age)
lx<-1
for(i in 2:(length(si))) lx[i]<-lx[i-1]*si[i-1]


mx <- mean.fec.mat
#mx <- mean.noff
lxmx<-lx*mx
x <- 1:6
# Now we are cooking!
junk<-nlminb(start = 0.1, obj = eulerlotka)




# Or we can calculate this each year and see what happens.
n.years <- nrow(age.mat)
years <- age.mat$Year
res <- data.frame(years = NA, r = NA)
for(j in 1:(n.years-1))
{
  
  lx <- exp(-nat.mort[j,2:7])
  #lx <- rep(0.2,6)
  si <- 1-lx
  max.age <- 6
  
  x<-1:(max.age)
  lx<-1
  for(i in 2:(length(si))) lx[i]<-lx[i-1]*si[i-1]
  
  
  mx <- fec.mat[j,]
  #mx <- mean.noff
  #lxmx<-lx*mx
  x <- 1:6
  # Now we are cooking!
  junk<-nlminb(start = 0.1, obj = eulerlotka)
  res[j,] <- c(years[j],junk$par)
  
}

ggplot(res) + geom_line(aes(x=years,y=r)) + theme_bw()
mean(res$r)

# Method #2 using the SPR code to get our lotka r values
# Here we can directly use the FishMethods package to get what we need
