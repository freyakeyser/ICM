# OK, so using the ICES assessments here's what we get for North Sea cod.
library(readxl)
library(tidyverse)
library(rio)
library(ggthemes)

# OK, lets give this a spin...
funs <- c("https://raw.githubusercontent.com/Dave-Keith/ICM/main/Scripts/functions/backwards_sim.r")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}

# Choose 5 ICES stocks that we have the necessary data for
ASR <- read_xlsx("./Data/ASR_2018.xlsx", sheet = "ICES")
datatypes <- unique(gsub(x = names(ASR), pattern = "[^a-zA-Z]", replacement=""))
# we want:
# Year, Num, WA, Catch, AM, NM, StockID, Management, Area, Order, Family, Genus, Species
ASRdat <- ASR[,c(grep(x=names(ASR), "Num"), 
                 grep(x=names(ASR), "WA"), 
                 grep(x=names(ASR), "Catch"),
                 grep(x=names(ASR), "AM"),
                 grep(x=names(ASR), "NM"))]
ASRdat <- apply(X = ASRdat, 2, as.numeric)
ASRsp <- ASR[, which(!1:length(names(ASR)) %in% grep(x=names(ASR), ".", fixed=T))]
ASR_trim <- cbind(ASRsp, ASRdat)

# need a unique ID for stock
table(ASR_trim$Management, ASR_trim$Species)
ASR_trim$Stock <- paste0(ASR_trim$Management, "_", ASR_trim$Area, "_", ASR_trim$Genus, "_", ASR_trim$Species)

ASR_long <- ASR_trim %>%
  pivot_longer(!c("Management", "Area", "Order", "Family", "Genus", "Species", "Stock", "Year")) %>% 
  separate(col=name, into=c("type", "age"), sep = "\\.")

Stocks <- ASR_long %>%
  filter(!is.na(value)) %>%
  group_by(Stock, type) %>%
  summarize(count=length(unique(value))) %>%
  group_by(Stock) %>%
  summarize(types=length(unique(type))) %>%
  filter(types==5) %>%
  select(Stock)

ASR_stocks <- ASR_long %>%
  filter(Stock %in% Stocks$Stock) %>%
  filter(!is.na(value)) %>%
  arrange(Stock, Year, type, as.numeric(age))

Stocks <- ASR_stocks %>%
  group_by(Stock, Species, type) %>%
  summarize(ages=length(unique(age)),
            years=length(unique(Year))) %>%
  arrange(-years, -ages) %>%
  filter(!Species=="morhua") %>%
  distinct(Stock) %>%
  pull(Stock)

print(Stocks)

# Bring in the data
# Handy import function from the rio package... tricky part was getting the right filename here...
for(i in Stocks){
  print(i)
  ASR_sub <- ASR_long %>%
    select(Year, Stock, type, age, value) %>%
    filter(Stock==i)
  
  age.mat <- ASR_sub %>% filter(type=="AM") %>% rename(AM=value) %>% select(-type)
  nat.mort <-  ASR_sub %>% filter(type=="NM") %>% rename(NM=value) %>% select(-type)
  abund <- ASR_sub %>% filter(type=="Num") %>% rename(Num=value) %>% select(-type)
  weight.age <- ASR_sub %>% filter(type=="WA") %>% rename(WA=value) %>% select(-type)
  removals <- ASR_sub %>% filter(type=="Catch") %>% rename(Catch=value) %>% select(-type)
  
  data <- age.mat %>%
    full_join(nat.mort) %>%
    full_join(abund) %>%
    full_join(weight.age) %>%
    full_join(removals)
  
  data$available <- apply(is.na(data[, c("AM", "Num", "WA")]), 1, function(x) all(!x==T))
  data <- data[data$available==T,]
  
  # Tidy up the data for input...
  data$prop.nat.mort <- 1-exp(-data$NM)
  #prop.nat.mort[,-which(names(prop.nat.mort) %in% c("Year", "Stock", "type"))] <- 1-exp(-prop.nat.mort[,-which(names(prop.nat.mort) %in% c("Year", "Stock", "type"))])
  
  rem <- data %>% group_by(Year) %>% summarize(rem=sum(Catch,na.rm=T)) %>% pull(rem)
  #rowSums(removals[,-which(names(removals) %in% c("Year", "Stock", "type"))], na.rm=T)
  years <- data %>% pull(Year) %>% unique() %>% sort()
  N.end <- sum(data[data$Year==max(years),]$Num)
  #N.end <- rowSums(abund[nrow(abund),-which(names(abund) %in% c("Year", "Stock", "type"))], na.rm=T)
  #vpa.abund <- rowSums(abund[,-which(names(abund) %in% c("Year", "Stock", "type"))], na.rm=T)
  vpa.abund <- data %>% group_by(Year) %>% summarize(vpa=sum(Num)) %>% pull(vpa)
  
  # The real mx matrix, recruits produced per individual in each age class... Not perfect as I need to offset recruits/ssb, but close enough for the moment..
  #minage
  minage <- as.numeric(min(data$age))
  #maxage
  maxage <- as.numeric(max(data$age))
  #recruits
  annual <- data.frame(Year=data$Year[data$age==minage], recruits=data$Num[data$age==minage])
  #ssn
  data$ssn <- data$Num * data$AM
  #ssb
  data$ssb <- data$ssn * data$WA
  #tot.ssb
  annual <- data %>% group_by(Year) %>% summarize(tot.ssb = sum(ssb)) %>% left_join(annual)
  #r.p.ssb
  annual$r.p.ssb <- annual$recruits/annual$tot.ssb
  # recs.per.age
  data <- left_join(data, annual)
  data$recs.per.age <- data$ssb*data$r.p.ssb
  # mx
  data$mx <- data$recs.per.age/data$ssn/2 # Moms only! Dividing by around 8 really nails it for NS cod for whatever reason :-)
  # Something is wrong with the decline rate method, but the exponential and logistic are working pretty... pretty... pretty good...
  
  age.mat <- data %>% select("Year", "age", "AM") %>% pivot_wider(names_from=age, values_from = AM) %>% select(-Year)
  prop.nat.mort <- data %>% select("Year", "age", "prop.nat.mort") %>% pivot_wider(names_from=age, values_from = prop.nat.mort) %>% select(-Year)
  weight.age <- data %>% select("Year", "age", "WA") %>% pivot_wider(names_from=age, values_from = WA) %>% select(-Year)
  mx <- data %>% select("Year", "age", "mx") %>% pivot_wider(names_from=age, values_from = mx) %>% select(-Year) %>% as.data.frame()
  
  tst <- icm.sim(years = years,
                 mat.age = age.mat,
                 nm = prop.nat.mort,
                 w.age = weight.age,
                 ages = minage:maxage,
                 rems = rem,
                 fecund = mx,
                 N.end = N.end,
                 pop.model = 'exponential',
                 n.sims = 5,
                 sd.mat = 0.5,
                 sd.nm = 0.5,
                 sd.wt = 0.5,
                 sd.fecund = 0.5)
  
  #Combine the data 
  did.it.work <- data.frame(abund = c(vpa.abund,tst$Pop$abund),years = c(years,tst$Pop$years),sim = c(rep('VPA',length(years)),tst$Pop$sim))
  # Here get the Upper and lower 50% quantiles to make a functional boxplot
  quants <- did.it.work %>% dplyr::group_by(years) %>% dplyr::summarise(L.50 = quantile(abund,probs=c(0.25)),
                                                                        med = median(abund),
                                                                        U.50 = quantile(abund,probs=c(0.75)))
  
  # All the results in a plot
  print(ggplot() + 
          geom_line(data=did.it.work, aes(x=years,y=abund,color=sim)) +
          xlab("") + 
          ylab("Abundance (1000s)") + 
          geom_line(data=did.it.work %>% dplyr::filter(sim == "VPA"),aes(x=years,y=abund),color='black',size=2) +
          #scale_y_continuous(breaks = seq(0,3e6,by=5e5)) + scale_x_continuous(breaks = seq(1960,2025,by=5)) +  
          theme_few() + theme(legend.position = 'none') + scale_color_viridis_d(end = 0.75) +
          ggtitle(i))
  # Same thing but functional boxplots
  print(ggplot() + 
          geom_line(data=quants, aes(x=years,y= med)) + 
          geom_ribbon(data=quants, aes(x=years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') +
          geom_line(data=did.it.work %>% dplyr::filter(sim == "VPA"),aes(x=years,y=abund),color='black',size=2) +
          xlab("") + ylab("Abundance (1000s)") + 
          #scale_y_continuous(breaks = seq(0,3e6,by=5e5)) + scale_x_continuous(breaks = seq(1960,2025,by=5)) +  
          theme_few() + theme(legend.position = 'none') +
          ggtitle(i))
  # Plotting how the estimate of r changes over time (which you can do when we have time varying inputs (e.g. Wgt, M, or maturity at age)
  print(ggplot(tst$r) + 
          geom_line(aes(x=years,y=r,group=n.sims,color=n.sims)) + 
          xlab("") + 
          ylab("Estimated growth rate (r)") + 
          #scale_y_continuous(breaks = seq(0,1.5,by=0.1)) + 
          #scale_x_continuous(breaks = seq(1960,2025,by=5)) +  
          theme_few() + theme(legend.position = 'none') + 
          scale_color_viridis_c(end = 0.75) +
          ggtitle(i))
}
