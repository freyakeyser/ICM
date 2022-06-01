# Create quadratic formula function modified from some code I found on here https://dk81.github.io/dkmathstats_site/index.html
# For backwards calculating the population size from a logistic growth model a = r/K, b = -(1+r), and c = the population size you are back calculating from.
#
# a:  For our logistic example a = r/K
# b:  For our logistic example b = -(1+r)
# c:  The population size 'next year'


quadratic.solver <- function(a, b, c) 
{
  discriminant <- (b^2) - (4*a*c)
  if(discriminant < 0) return(paste0("This quadratic equation has no real numbered roots."))
  if(discriminant > 0) x.int <- c(format(round((-b - sqrt(discriminant)) / (2*a)), nsmall = 5),format(round((-b + sqrt(discriminant)) / (2*a)), nsmall = 5))
  if(discriminant == 0) format(round(x.int <- (-b) / (2*a)), nsmall = 5)
  return(x.int)
}
