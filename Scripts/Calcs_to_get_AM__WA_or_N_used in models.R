
# Pacific Cod in Bering Sea from Stark 2007.
logis.mod <- function(alpha,beta,x){1/(1+exp((alpha + beta*x)))}

x <- 0:10
alpha = 4.7143
beta = -0.9654

res <- logis.mod(alpha,beta,x)
res

#GOA Gadus macrocephalus, in Assessment document and from Stark 2007 as well.

x2 <- 0:20
alpha2 = 4.3
beta2 = -0.963 # There is a typo in the assessment document, should be 0.963 not 1.963...

res2 <- logis.mod(alpha2,beta2,x2)
res2

# AI Gadus Macrocephalus paramters, in assessment document

x3 <- 0:20
alpha3 = 4.883
beta3 = -0.965

res3 <- logis.mod(alpha3,beta3,x3)
res3
