# The analogue of the backword projection function.... but to forward project

# Arguments
#option:         So far we can do exponential or logistic growth, other options are on the table if someone has a good one
#pop.last:       The population size in year t-1, which is used to estimate the size of the population in year t
#r:              The intrinsic rate of growth for the population
#K:              The carrying capacity, only necessary for logistic model
                 
 
for.proj <- function(option = "exponential",pop.last,r,K)
{
  ### So rolling with the easy peasy continuous time exponential
  ### simple exponential population growth, in this formulation with the removals I think I'm saying the removals happen right at the outset
  if(option == 'exponential') Pop.current <- (pop.last)*((exp(r)))

  # Logistic growth model, 
  # The maths be.... N(t+1) = N(t) + rN(t)(1-N(t)/K).... now we are solving for N(t) not N(t+1), so we rearrange terms
  # to get it in the quadratic formula....  (r/K)N^2 -(1+r)N  + N(t+1) = 0... we have r, K, and N(t+1) so we can solve for N(t)
  # But #1, this is the discrete time version, and that's not useful
  # But thanks to my longtime biological math hero Sarah Otto and her book on Mathematical modelling
  # This will solve the backwards logistic and is easy peasy to use for forwards logistic
  # Hurray for math!
  if(option == "logistic")   Pop.current <-  (exp(r)* pop.last)*(1-(pop.last/K)) 
  
  return(list(Pop.current = Pop.current))
}
