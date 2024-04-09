# Very simple model that projects the population forward using exponential or logistic models

# Arguments
#option:         So far we can do exponential or logistic growth, other options are on the table if someone has a good one
#pop.last:       The population size in year t-1, which is used to estimate the size of the population in year t
#r:              The rate of growth for the population
#K:              The carrying capacity, only necessary for logistic model
                 
 
for.proj <- function(option = "exponential",pop.last,r,K)
{
  ### So rolling with the easy peasy continuous time exponential
  if(option == 'exponential') Pop.current <- (pop.last)*((exp(r)))
  # Logistic growth model, 
  if(option == "logistic")   Pop.current <-  (exp(r)* pop.last)*(1-(pop.last/K)) 
  return(list(Pop.current = Pop.current))
}
