setwd("D:/Github/ICM/Data/")

if(!require(stringr)){install.packages("stringr"); library(stringr)}
#devtools::install_github("james-thorson/FishLife")
install.packages("remotes")
if(!require(Rcpp)){install.packages("Rcpp"); library(Rcpp)}
remotes::install_github("James-Thorson/FishLife",force=TRUE)
if(!require(FishLife)){install.packages("FishLife"); library(FishLife)}
#vignette("tutorial","FishLife")

library(tidyverse)
library(FishLife)
# species list
species<-read.table("specieslist.txt",header=T)
head(species)
speciesclear<-data.frame(str_split(species$genus.species,"_", simplify = TRUE))[,1:2];colnames(speciesclear)<-c("genus","species")
speciesclear$genus<-as.character(speciesclear$genus)
speciesclear$species<-as.character(speciesclear$species)
dim(speciesclear)

listspecies<-seq(1,nrow(speciesclear),1)
output<-matrix(0,nrow(speciesclear),20)
for (i in listspecies) 
{ 
  Predict <-try(Plot_taxa(Search_species(Genus=as.character(speciesclear$genus[i]),
                                         Species=as.character(speciesclear$species[i]))$match_taxonomy),silent = TRUE)
                         #plot(seq(1,3,1),seq(1,3,1),main=i)
                         #graphics.off()
                         output[i,]<-try(Predict[[1]]$Mean_pred,silent = TRUE)
                         if(output[i,1]==c("Error in Predict[[1]]$Mean_pred : \n  $ operator is invalid for atomic vectors\n")){output[i,]<-"NA"}
                         else {output[i,]<-output[i,]}
                         print(paste(round((100*i)/length(listspecies),digits=2),"%"))
                          }
                         output<-data.frame(output)
                         colnames(output)<-names(Predict[[1]]$Mean_pred)
                         rownames(output)<-species$genus.species
                         
                         colnames(output)<- c("Asymptotic_length","Brody_growth_coefficient","Asymptotic_mass_cm","Maximum_age_year","Age_at_maturity_year",
                                              "Mortality_rate",
                                              "Length_at_maturity_cm",
                                              "Average_Temperature",
                                              "Conditional_recruitment_variance",
                                              "Recruitment_autocorrelation",
                                              "Maximum_annual_spawners_per_spawner", 
                                              "SD_of_recruitment",
                                              "Steepness_h",
                                              "logit_Steepness_h",
                                              "Ratio_of_F_msy_and_M",
                                              "Fishing_mortality_rate_at_MSY",
                                              "intrinsic_growth_rate_pop_growth",
                                              "intrinsic_growth_rate",
                                              "ln_Generation_time",
                                              "Geration_time")
                         
                         output[,1:7]<-exp(output[,1:7])
                         head(output)
                         write.table(output,file="output.txt")    
                         
        