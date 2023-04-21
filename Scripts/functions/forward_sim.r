# here we project forward, should be fairly straightforward tweak of the backward simulation function

# All the data you need to make it dance...
# years:          The years you are running the backwards calculation for
# n.sims:         The number of simulations you are running
# mat.age         Age at maturity if one value it is the age at 50% maturity, if a vector it is the age at maturities for each age(or age/year), if a matrix we want this to be unique for each simulation
#                 we could also let this vary by year, but we'll need to go 3-D array or something for that# ages:           What are the ages you are using this is used to calculate max age, which may not be ideal with plus group stocks.

# nm:             Natural mortality, this wants the instantaneous natural mortality, not the proportional. There are several options here....
#                 You can enter a vector of instantaneous natural moralities, just put a single number of all ages, or put in a character
#                 string that will tell the function to calculate the natural mortality based on life history data. If a matrix the rows are different simulations
#                 we could also let this vary by year, but we'll need to go 3-D array or something for that

# sel             The selectivity of the stock.  Currently this is set up as a single number and it is tweaked in the code, we'll need to make this more complex.

# rems            Projected removals.  Currently this is set up so that you provide the mean and sd of fishing mortality for the forward projections.  Will need to generalize this.

# N               Population size in a given year

# u:              Exploitation rate (annual not instantaneous), currently set up to be 1 value.

# pop.model:      What method you going to use to get the population growth modeled for the backwards model.  You can use exponential model, logistic model, 
#or a dec.rate (decline rate) model.

# What we need for the Lotka.r function...
# fecund:      How are you estimating fecundity. If a vector it is the number of recruits produced by the average individual in the age classes. 
#              If a matrix the rows are different simulations. Also have options to set fecund = 'eggs' which we can use if we can get an estimate of egg mortality
#              fecund = 'SPR' uses the spawner per recruit metric, to get this to work we need to know how many recruits there are in a given year
#              We could also let this vary by year, but we'll need to go 3-D array or something for that
# mat.ogive.K 
# L.inf 
# K 
# t0 
# a.len.wgt 
# b.len.wgt 
# a.fec.len 

# w.age           The weight of individuals for each age, if a matrix the rows are different simulations
#                 we could also let this vary by year, but we'll need to go 3-D array or something for that
# K
# dec.rate
# sim
# proj.sim

for.sim<-function(years,n.sims=1,mat.age = NULL,nm=NULL,w.age = NULL,ages =NULL,fecund = NULL,N.start = NULL,
                  sel,rems = c(0.1,0.2),N,u,pop.model = "exponential", sim = 'retro', proj.sim='dist',
                  dec.rate = NULL,
                  L.inf = NULL,K = NULL,t0 = NULL, a.len.wgt = NULL, b.len.wgt = NULL, a.fec.len = NULL, b.fec.len = NULL,
                  sd.mat = 0,sd.nm = 0,sd.wt = 0,sd.fecund = 0)
{
  # Download the function to go from inla to sf
  funs <- c("https://raw.githubusercontent.com/freyakeyser/ICM/main/Scripts/functions/Lotka_r.r",
            "https://raw.githubusercontent.com/freyakeyser/ICM/main/Scripts/functions/backwards_project.r"
  )
  # Now run through a quick loop to load each one, just be sure that your working directory is read/write!
  for(fun in funs) 
  {
    download.file(fun,destfile = basename(fun))
    source(paste0(getwd(),"/",basename(fun)))
    file.remove(paste0(getwd(),"/",basename(fun)))
  }
  
  source("D:/Github/ICM/Scripts/functions/Lotka_r.r") # For testing purposes, delete when done
  
  st.time <- Sys.time()
  require(optimx)  || stop("Please load the 'optimx' package which you'll need for the optimations to run")
  
  #Initialize a bunch of objects
  n.years<-length(years)
  Pop<-data.frame(abund = rep(NA,n.years*n.sims),sim = rep(1:n.sims,n.years),years = sort(rep(years,n.sims)))    #matrix(data=NA, nrow=<<see below>>, ncol=<<see below>>,byrow=F, dimnames=NULL) 
  Pop.vec<-rep(NA,n.years)
  rem<-rep(NA,n.years)
  r.vec<-NULL
  B1.vec<-rep(0,n.sims) ### backwards projections 
  #rems<-matrix(rep(0,n.years*n.sims),n.sims,n.years)

  #Calculate r for your example
  for(i in 1:n.sims)
  {
    # For the first run we use the mean estimate, for runs after that we add in all the uncertainty specified
    # For now I've only tested this using the stock assessment data, needs cleaned up for other options.
      junk<-lotka.r(yrs = years,age.mat = mat.age,nat.mort = nm,ages=ages,wt.at.age=w.age,fecund=fecund,sim =sim,proj.sim = proj.sim,
                    L.inf = L.inf,K = K,t0 = t0, a.len.wgt = a.len.wgt, b.len.wgt = b.len.wgt, a.fec.len = a.fec.len, b.fec.len = b.fec.len,
                    sd.mat = sd.mat,sd.nm = sd.nm,sd.wt = sd.wt,sd.fecund = sd.fecund)   
    #browser()
    tmp <-junk$res[,2] 
    if(length(tmp) == 1) tmp <- rep(tmp,n.years)
    r.vec[[i]] <- data.frame(r = tmp,years = years[],n.sims=i)

  }
  print(Sys.time() - st.time)
  st.time2 <- Sys.time()
  #unwrap your r vector
  r.vec <- do.call('rbind',r.vec)
  #temp.r.vec<-r.vec[r.vec>0 & r.vec<r.cutoff]
  #r.vec<-temp.r.vec[1:n.sims]
  ############################################################
  
  #then do for each population simulation
  
  for(i in 1:n.sims)
  {
    
    # Get the final year estimate of your population
    Pop.vec[1] <- pop.last <- N.start  # Assuming this is our known starting point, could add uncertainty to this as well if we want to.
    r.tmp <- r.vec[r.vec$n.sims == i,]
    #Now run your model backwards with the r from the Lotka function and
    # if you use the logistic growth model the K estimated for the population.
    for(y in 2:n.years)
    {
      #browser()
      # DK Note: So removals 
      if(sim == 'retro') removals.next <- rems[y]
      if(sim == 'project') 
      {
        fm <- summary(rlnorm(1,log(rems[1]),rems[2]))
        removals.next <- fm*pop.last
      }
      #browser()
      r.up <- r.tmp$r[y] # 
      # The exponential model
      if(pop.model == 'exponential') 
      {
        exp.res <- for.proj(option = "exponential",pop.last = pop.last ,r=r.up,removals = removals.next,direction,fishery.timing = 'beginning')
        pop.last <- exp.res$Pop.current
      }
      # If you are running the logistic model
      if(pop.model == 'logistic')
      {
        log.res <- for.proj(option = "logistic",pop.last = pop.last,K=K,r=r.up,removals = removals.next)
        pop.last <- min(log.res$Pop.current)
      }
      if(pop.model == 'dec.rate')
      {
        # This needs way more thought than it's been given by DK!
        #removals[i,]<-(-1+exp(r.vec[i])+dec.rate[i]) ## redone with Jamie.
        pop.last <- pop.last + pop.last*(r.up+dec.rate[i])
      }
      #browser() 
      Pop.vec[y] <- pop.last

    } #Loop through all the years.
    
    Pop$abund[Pop$sim == i]<-Pop.vec
    
    ### calc removals (increase due to r plus amount population declined by):   
    # removals[i,]<-Pop.vec-(Pop.vec/exp(r.vec[i]))+(Pop.vec*(decl.rate[i])) WRONG
    
    ## removals in final year of population decline to calculate forward exploitation rate   
    #rem<- removals[n.years]
    
    #print(Pop.vec)
    #B1.vec[i]<-lm(log(Pop.vec)~years)$coef[2]  ## logic check
    
    
    #calculate exploitation
    #junk<-c()
    #for(y in 1:n.years) junk[y]<-u.calc(nat.mort.vec[i],juv.mult.vec[i],max.age.vec[i],age.mat.vec[i],sel,rem*0.7,Pop.vec[y])
    #mean.u.vec[i]<-mean(junk) 
  } #end backwards projection
  print(Sys.time() - st.time2)
    return(list(Pop=Pop,r = r.vec))
}
  
  