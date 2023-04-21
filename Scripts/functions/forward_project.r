# The analogue of the backword projection function.... but to forward project

# Arguments
#option:         So far we can do exponential or logistic growth, other options are on the table if someone has a good one
#pop.last:       The population size in year t-1, which is used to estimate the size of the population in year t
#r:              The intrinsic rate of growth for the population
#removals:       The removals between t+1 and t
#K:              The carrying capacity, only necessary for logistic model
#fishery.timing: When do you want the removals by the fishery to occur, at the 'beginning' of the year, or at the end of the year.  
#                   At the 'beginning' excludes them from the population dynamics, at the 'end' includes them in the population dynamics (because we are going backwards).
 
for.proj <- function(option = "exponential",pop.last,r,removals.next,K,fishery.timing = 'beginning')
{

  
  ### So rolling with the easy peasy continuous time exponential
  ### simple exponential population growth, in this formulation with the removals I think I'm saying the removals happen right at the outset
  if(option == 'exponential') 
  {
    if(fishery.timing == 'end')
    {
    Pop.current <- (pop.last)*((1+r))
    Pop.current<- Pop.current-removals.next 
    }
    if(fishery.timing =='beginning') Pop.current <- (pop.last-removals.next)*((1+r))
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
    if(fishery.timing == 'end') par <- pop.last
    if(fishery.timing == 'beginning') par <- pop.last <- pop.last - removals.next
    # Get the current population size using logistic model
    Pop.current <-  pop.last + (r* pop.last)*(1-(pop.last/K)) 
    # If the fishery removes individuals at the start of the time series add them in here.
    if(fishery.timing == 'end') Pop.current <- Pop.current - removals.next
  }
  return(list(Pop.current = Pop.current))
}
