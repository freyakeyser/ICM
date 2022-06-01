# This is a function to estimate population history over time currently using logistic or exponential growth models
library(tidyverse)

# Load the backwards projection via github
fun <- c("https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/backwards_project.r")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
download.file(fun,destfile = basename(fun))
source(paste0(getwd(),"/",basename(fun)))
file.remove(paste0(getwd(),"/",basename(fun)))

# OK, so lets compare trajectories of the backwards logistic and backwards exponential
# Play around with pop.next r, and K, most interesting.  Note that if K <= our the value we start at (i.e. the most recent year) 
# we end up in some weird scenarios (if we are = K it's deterministic so we don't go anywhere (r effectively = 0), 
# and if we are above K, unless you really crank R up (to get into 'chaos' territory), the backwards calculation will 
# say that the population is > K for the whole time series, because the logistic with low r above K, just 
# declines to K smoothly
# If tossed these being deterministic and added some noise we could get some more entertaining results.
# Seems like pulling the 'min' value works just fine for the logistic as the max() solution is generally garbage when r is low
# but as r increases the 'max' solution becomes increasingly viable, just watch the logistic.unreal trend as you increase r.

years <- 1980:2020
n.years <- length(years)
pop.next <- 40000
K = 500000
r=0.1
eff = 0.1 # Proportional removals.
res <- data.frame(year = years, 
                  exponential = c(rep(NA,n.years-1),pop.next),
                  rem.exp = c(rep(NA,n.years)),
                  exp.change = c(rep(NA,n.years)),
                  logistic = c(rep(NA,n.years-1),pop.next),
                  logistic.change =  c(rep(NA,n.years)), 
                  rem.log = c(rep(NA,n.years)),
                  logistic.unreal = c(rep(NA,n.years-1),pop.next),
                  logistic.unreal.change = c(rep(NA,n.years)))

# So a really interesting problem for back calculation of the logistic model.  As you approach K, r becomes increasingly small
# so if you 'start' the model anywhere near enough to K then you are going to have very little growth.  Buuut if you have removals
# that are higher than the growth, because we are going 'backwards' that means the population was larger last year than it was this year
# so what happens in the logistic model is that the backwards trajectory suggests that the biomass was higher in the past
# and in short order the biomass estimate is bizzaro high.  Basically if you are within N = 0.8K and F > 0.5r the logistic is going to go off
# into weird space (and probably not work) rather quickly.

# Meanwhile, if you F is similar to R the models also behave very differently than you might expect.
# If F << r then you get nice clean exponential and (sometimes) logistic growth (depending on how close K is to N), which makes sense since F is low relative to R
# If F ??? r then you get patterns that do not at all resemble logistic or exponential growth (which is good!)
# If F >> r, the models start to look like exponential decline to current stock levels cause, well you are overfishing them. But the logistic model
# is especially sensistive here and can really go into bizarro population sizes quickly.  Worth noting that even with a very high K
# the logistic model and the exponential models can end up in very different places

for(i in n.years:2)
{
  if(i == n.years) removal.log = rlnorm(1,log(eff*pop.next),sd=0.1)
  if(i == n.years) removal.exp = rlnorm(1,log(eff*pop.next),sd=0.1)
  if(i < n.years) removal.log = rlnorm(1,log(eff*pop.log.next),sd=0.1)
  if(i < n.years) removal.exp = rlnorm(1,log(eff*pop.exp.next),sd=0.1)
  
  # Start with the exponential model
  if(i != n.years) pop.next = pop.exp.next
  exp.res <- back.proj(option = "exponential",pop.next = pop.next,K=K,r=r,removals = removal.exp)
  pop.exp.next <- exp.res$Pop.current
  res$exponential[i-1] <- pop.exp.next
  res$rem.exp[i] <- removal.exp
  res$exp.change[i] <- exp.res$Pop.ops-res$exponential[i]
  res$exp.F <- res$rem.exp[i]/ res$exponential
  # Run backwards through the logistic model
  if(i != n.years) pop.next = pop.log.next
  log.res <- back.proj(option = "logistic",pop.next = pop.next,K=K,r=r,removals = removal.log)
  pop.log.next <- min(log.res$Pop.current)
  res$logistic[i-1] <- pop.log.next
  res$logistic.change[i] <- min(log.res$Pop.ops) - res$logistic[i]
  res$rem.log[i] <- removal.log
  res$logistic.unreal[i-1] <- max(log.res$Pop.current)
  res$logistic.unreal.change[i] <- max(log.res$Pop.ops) - res$logistic.unreal[i]
}

res.long <- pivot_longer(res,cols =c('logistic','exponential','logistic.unreal'),names_to = "Model",values_to = "Abundance")
# Interesting how an F of even 5% seems to obliterate the shape of the logistic curve and makes the logistic model just look like an exponential model
# I haven't wrapped my head around what that is...
ggplot(res.long %>% dplyr::filter(Model != 'logistic.unreal')) + geom_line(aes(x=year,y=Abundance,color=Model),size=2) + scale_color_manual(values = c("blue","orange"))

# The unreal solution, but check out what happens as you increase r, the minimum solution increasingly gets silly instead of the max, and in reality the time series
# could be made up of a hodge-podge of mins and maxs, but I don't think we'll need to worry about it at r values we'll encounter...
# but still stewing on it! Curious that at high r values the population stabilizes at a N that is much higher than K, what that N level is, is a function of r.. 
# so N = [(1/r)+1]K when we pull the max values only.... there's a proof I should know in there somewhere isn't there :-D
ggplot(res.long) + geom_line(aes(x=year,y=Abundance,color=Model),size=2) + 
  scale_color_manual(values = c("blue","orange",'red'))

#r = 2 was 46500 or 1.5K
#r = 1 was 62000 or 2K
#r = 0.5 was 93000 or 3K
#r = 0.25 was 155000 or 5K
# N = (1/r +1)K





