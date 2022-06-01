# This is a function to estimate population history over time currently using logistic or exponential growth models
# But it starts at the end of the time series and projects backwards to the start of the time series. Requires 
# and estimate of r, K(logistic option only), intrinsic rate of growth (r) and the removals for the year

# Arguments
#option:        So far we can do exponential or logistic growth, other options are on the table if someone has a good one
#pop.next:      The population size in year t+1, which is used to estimate the size of the population in year t
#r:             The intrinsic rate of growth for the population
#removals:      The removals between t+1 and t
#K:             The carrying capacity, only necessary for logistic model
back.proj <- function(option = "exponential",pop.next,r,removals.next,K)
{

  
  ### So first we have the discrete time exponential model, n(t+1) = R*n(t) because we only have data at discrete time steps to update our model.
  ### Where r = R-1, so rearranging r+1 is what we want for this model, I think?
  ### simple exponential population growth, in this formulation with the removals I think I'm saying the removals happen right at the outset
  if(option == 'exponential') 
  {
    pop.ops <- (pop.next)/(1+r)
    Pop.current<- pop.ops+removals.next 
  }
  # Logistic growth model, 
  # The maths be.... N(t+1) = N(t) + rN(t)(1-N(t)/K).... now we are solving for N(t) not N(t+1), so we rearrange terms
  # to get it in the quadratic formula....  (r/K)N^2 -(1+r)N  + N(t+1) = 0... we have r, K, and N(t+1) so we can solve for N(t)
  # But of course that means there are two possible solutions using the logistic model, the minimum solution from the quadratic
  # is the biologically reasonable solution, the maximum becomes more reasonable as you increase r, but is generally biologically irrelvant
  # at least in the parameter space we are interested in.
  # As with the exponential model, I think doing the removals like this means that they are happening right at the start of the year and don't 
  # get to contributed to growth/reproduction of the population
  if(option == "logistic")
  {
    if(!exists('quadratic.solver')) # Just so you only download this once

    {
      # Get the quadratic solver function
      fun <- c("https://raw.githubusercontent.com/Dave-Keith/ICM/master/Scripts/functions/quadratic_solver.r")
      # Now run through a quick loop to load each one, just be sure that your working directory is read/write!
      download.file(fun,destfile = basename(fun))
      source(paste0(getwd(),"/",basename(fun)))
      file.remove(paste0(getwd(),"/",basename(fun)))
    }
    
    pop.ops <- as.numeric(quadratic_solver(a = (r/K),b =-(1+r), c = pop.next))
    Pop.current<-c(pop.ops)+removals.next
  }
  return(list(Pop.current = Pop.current,Pop.ops = pop.ops))
}
