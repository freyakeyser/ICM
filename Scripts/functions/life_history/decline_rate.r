# A function to caclucate the rate of decline, currently it just is sampling from estimates from Baum and Curtis, which are of course irrelvant for anything other than white shark... so that needs to change

#option:  Which method do you want to use, "Curtis" or "Baum", both pointless for anything but white shark, but want to see if I can getter to work.
#n.sims:  Number of simulations you are running
dec.fun <- function(option = "Baum",n.sims = n.sims)
{
  if(option == "Baum")
  {
    years<<-1986:2000 #15 years
    decl.junk<-seq(1.058,1.137,by=0.0001)  #  200 values - note that this is the old parameterization of decline
    decl.rate<<-sample(decl.junk,n.sims*4,replace=T) 
    ## decline rate 59% over 15 years to 89% over 15 years -Baum et al. 2003)
    #to calculate annual rate from total decline:
    #1-.59 = 0.41 to 1-.89 = 0.11  total decline (beginning value = 1; ending value = 0.41) ending/beginning
    #0.41^(1/14)-1 = -0.05770781 per year ### NOTE - 14 time steps over 15 years
    #0.11^(1/14)-1 = -0.1368369 per year  
  } # End Baum Option
  
  if(option == "Curtis")
  {
    years<<-1970:1989 ## 20 years
    decl.junk<-seq(0.05098,0.06659,by=0.0001)  #  approx 150 values
    decl.rate<<-sample(decl.junk,n.sims*4,replace=T) 
    ## 63% to 73% decline reported in Curtis et al. (2014)
    #to calculate annual rate from total decline:
    #1-.63 = 0.37 to 1-.73 = 0.27  total decline
    #0.37^(1/19)-1 =  -0.05098347 per year ## 19 timesteps over 20 years.
    #0.27^(1/19)-1 =  -0.06659144 per year 
  }# End Curtis option
}  
