# Two functions to start the festivities...
# The eulerlotka function
eulerlotka <- function(r) (sum(lx * mx * exp(-r * ages)) - 1)^2
# And it's friend to get the eigen values from a Leslie Matrix

eigenfun<- function(A)
{
  # Get the eigen values
  lambdas <- eigen(A)
  lmax <- which(Re(lambdas$values) == max(Re(lambdas$values)))
  lmax <- lmax[1]
  # Grab the largest eigenvalue
  lambda1 <- Re(lambdas$values[lmax])
  return(list(lambda1=lambda1,lambdas = lambdas)) # this output the variables of interest for future fun
}

ages <- 1:10 # ages
mx <- c(0,0,rep(0.25,8)) # fecundity
nat.mort <- rep(0.89,10) # natural mortality

si <- 1-nat.mort # Survivorship
# Turn this into a Leslie matrix...
Leslie <- matrix(data = 0,nrow=10,ncol=10)
Leslie[1,] <- mx
for(i in 1:9) Leslie[i+1,i] <- si[i]
#Leslie[10,10] <- si[10] 
# If you set it to 0 then everyone dies at age 10, which is fine too, doesn't make much difference here.
Leslie[10,10] <- 0
# Get the lx matrix for Euler-Lotka
lx <- 1
for(s in 2:(length(si))) lx[s]<-lx[s-1]*si[s-1]

# Now solve the Euler-lotka add 1 to make it lambda
junk<-optimize(lower=-15,upper=15,f = eulerlotka)
res <- c(junk$minimum)
lambda.euler <- exp(res) # exp you fool!!
lambda.euler

# Now solve your Leslie matrix
res.leslie <- eigenfun(Leslie)
lambda.leslie <- res.leslie$lambda1
lambda.leslie

ages <- 0:8
nat.mort <- c(2.4,2.4,2.4,1.4,0.85,0.75,0.51,0.43,0.39)
mx <- c(0,0,0,9,14,23,22,23,37)
si <- exp(-nat.mort) # Survivorship
lx <- 1
for(s in 2:(length(si))) lx[s]<-lx[s-1]*si[s-1]

# Now solve the Euler-lotka add 1 to make it lambda
junk<-optimize(lower=-15,upper=15,f = eulerlotka)
res <- c(junk$minimum)
res