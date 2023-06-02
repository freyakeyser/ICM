# This is a function to estimate population history over time currently using logistic or exponential growth models
# But it starts at the end of the time series and projects backwards to the start of the time series. Requires 
# and estimate of r, K(logistic option only), intrinsic rate of growth (r) and the removals for the year

# Arguments
#option:         So far we can do exponential or logistic growth, other options are on the table if someone has a good one
#pop.next:       The population size in year t+1, which is used to estimate the size of the population in year t
#r:              The intrinsic rate of growth for the population
#removals:       The removals between t+1 and t
#K:              The carrying capacity, only necessary for logistic model
#fishery.timing: When do you want the removals by the fishery to occur, at the 'beginning' of the year, or at the end of the year.  
#                   At the 'beginning' excludes them from the population dynamics, at the 'end' includes them in the population dynamics (because we are going backwards).

back.proj <- function(option = "exponential",pop.next,r,removals.next,K,fishery.timing = 'beginning')
{

  
  ### So rolling with the easy peasy continuous time exponential
  ### simple exponential population growth, in this formulation with the removals I think I'm saying the removals happen right at the outset
  if(option == 'exponential') 
  {
    if(fishery.timing == 'beginning')
    {
    # Going to 1+r model since we are discrete time...
    #Pop.current <- (pop.next)/(exp(r))
    Pop.current <- (pop.next)/(1+r)
    Pop.current<- Pop.current+removals.next 
    }
    if(fishery.timing =='end') Pop.current <- (pop.next+removals.next)/((1+r)) #  (pop.next+removals.next)/(exp(r))
  }
  # Logistic growth model, 
  # The maths be.... N(t+1) = N(t) + rN(t)(1-N(t)/K).... now we are solving for N(t) not N(t+1), so we rearrange terms
  # to get it in the quadratic formula....  (r/K)N^2 -(1+r)N  + N(t+1) = 0... we have r, K, and N(t+1) so we can solve for N(t)
  # But #1, this is the discrete time version, and that's not useful
  # But thanks to my longtime biological math hero Sarah Otto and her book on Mathematical modelling
  # This will solve the backwards logistic and is easy peasy to use for forwards logistic
  # Hurray for math!
  if(option == "logistic")
  {
    # Get our parameter for the logistic model
    if(fishery.timing == 'beginning') par <- pop.next
    if(fishery.timing == 'end') par <- pop.next <- pop.next + removals.next
    # Ha, we can use this to solve for N in the previous year 
    #I'm hung up on which of these to use, they give almost the same answer, but not quite, but going with the 1+r version since
    # we are in discrete time here.
    #logistic.n.last <- function(N.last) sum((exp(r)* ((N.last)/(1-(N.last/K))) / (1 + ( ((N.last)/(1-(N.last/K)))*exp(r))/K) - pop.next))^2
    logistic.n.last <- function(N.last) sum(((1+r)* ((N.last)/(1-(N.last/K))) / (1 + ( ((N.last)/(1-(N.last/K)))*(1+r))/K) - pop.next))^2
    
    # Solve the above function
    Pop.current <- optimx(par, logistic.n.last, method = "BFGS")$p1
    # If the fishery removes individuals at the start of the time series add them in here.
    if(fishery.timing == 'beginning') Pop.current <- Pop.current+removals.next
  }
  return(list(Pop.current = Pop.current))
}
