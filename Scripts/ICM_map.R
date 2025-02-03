# Map for ICM

require(sf)
require(tidyverse)
require(ggplot2)
require(rnaturalearth)
require(ggplot2)

stocks <- read.csv("C:/users/keyserf/Documents/Github/ICM-Dave/Data/metadata.csv")
world <- ne_countries(scale = "large", returnclass = "sf")
ggplot() + geom_sf(data=world) + coord_sf(expand=F)

unique(stocks$Area.long)

bering <- c(60.60575, -169.41152)
alaska <- c(58.58383, -146.74756)
GB <- c(41.58437, -67.25976)
NL <- c(45.87445, -58.24261)
north <- c(56.33245, 3.46180)
celtic <- c(49.11854, -8.31070)
baltic <- c(58.49575, 20.07792)
faroe <- c(61.84574, -7.63399)
ice <- c(63.35810, -11.95852)
irish <- c(53.83017, -5.28181)
biscay <- c(45.49964, -3.48611)
maine <- c(43.30045, -68.46818)
us <- c(37.88852, -70.45552)
arctic <- c(72.13541, -2.38604)
source("C:/Users/keyserf/Documents/Github/Assessment_fns/Maps/github_spatial_import.R")
nafo <- github_spatial_import("NAFO", "Divisions.zip")
ggplot() + geom_sf(data=nafo, aes(fill=Division))

nafo_sub <- nafo[nafo$Division %in% c("4T", "5Y", "5Z"),]
nafo_sub <- nafo_sub %>% group_by(Division) %>%
  summarize() %>%
  st_centroid()
nafo_sub$lon <- st_coordinates(nafo_sub)[,1]
nafo_sub$lat <- st_coordinates(nafo_sub)[,2]

coords <- data.frame(Area.long = unique(stocks$Area.long),
           lat = c(bering[1],
                   bering[1],
                   bering[1],
                   alaska[1],
                   GB[1],
                   NL[1],
                   NL[1],
                   arctic[1],
                   north[1],
                   celtic[1],
                   baltic[1],
                   ice[1],
                   ice[1],
                   irish[1],
                   celtic[1],
                   biscay[1],
                   GB[1],
                   GB[1],
                   GB[1],
                   GB[1]),
           lon = c(bering[1],
                   bering[1],
                   bering[1],
                   alaska[1],
           GB[2],
           NL[2],
           NL[2],
           arctic[2],
           north[2],
           celtic[2],
           baltic[2],
           ice[2],
           ice[2],
           irish[2],
           celtic[2],
           biscay[2],
           GB[2],
           GB[2],
           GB[2],
           GB[2]),
           ID = c(1, 1, 1, 2,
                  3,
                   4,
                   4,
                  5,
                   6,
                   7,
                  8,
                   9,
                  9,
                   10,
                   7,
                   11,
                   3,
                   3,
                   3,
                   3))

stocks <- left_join(stocks, coords)
stocks <- st_as_sf(stocks, coords=c("lon", "lat"), crs=4326, remove=F)

ggplot()+ geom_sf(data=nafo, aes(fill=SubArea)) +
  geom_sf_text(data=stocks, aes(label=ID))



sf_use_s2(FALSE)
worldcrop <- st_crop(world, xmin=-180, ymin=30, xmax=30, ymax=75)

locs <- stocks %>% group_by(ID,lat, lon) %>%
  summarize(nstocks=n())

# locs %>% group_by(ID) %>%
#   summarize(sum(nstocks))

cols <- c('#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4')
names(cols) <- sort(unique(locs$Order))
require(magick)
barfig<- NULL
for (i in 1:11){
  barfig[[i]] <- magick::image_graph(width = 50, height = 50, res = 72, bg="transparent")
  print(ggplot() + geom_bar(data=locs[locs$ID==i,], aes(x=1, y=nstocks),
                      position="stack", stat="identity") +
    #scale_fill_manual(values=cols, guide=F) +
    ylim(0,15)+ # take this lim from the locs summary above
    xlab("") + ylab("") +
    theme_bw() +theme(axis.text=element_blank(), axis.ticks = element_blank(), panel.background=element_blank(), panel.grid=element_blank()))
  dev.off()
  barfig[[i]] <- image_transparent(barfig[[i]], 'white')
}

basemap <- image_graph(width = 850, height = 550, res = 96)
ggplot() + geom_sf(data=worldcrop) +
  #geom_sf_text(data=locs, aes(label=ID)) +
  coord_sf(expand=F) + theme_minimal() + theme(panel.grid = element_blank()) + xlab("")+ ylab("")
dev.off()

# all <- ggplot() + geom_bar(data=locs, aes(x=1, y=nstocks, fill=Order),
#                                  position="stack", stat="identity") +
#                scale_fill_manual(values=cols,name=NULL) +
#                xlab("") + ylab("") +
#                theme_bw() +theme(axis.text=element_blank(), axis.ticks = element_blank(), panel.background=element_blank())
# legend <- cowplot::get_legend(all)
# magickleg <-  magick::image_graph(width = 140, height = 150, res = 96, bg="transparent")
# plot(legend)
# dev.off()
# magickleg <-image_transparent(magickleg, "white")

final <- image_composite(basemap, barfig[[1]], offset = "+65+190")
final <- image_composite(final, barfig[[2]], offset = "+160+205")
final <- image_composite(final, barfig[[3]], offset = "+455+305")
final <- image_composite(final, barfig[[4]], offset = "+485+275")
final <- image_composite(final, barfig[[5]], offset = "+690+120")
final <- image_composite(final, barfig[[6]], offset = "+710+210")
final <- image_composite(final, barfig[[7]], offset = "+660+260")
final <- image_composite(final, barfig[[8]], offset = "+770+200")
final <- image_composite(final, barfig[[9]], offset = "+655+175")
final <- image_composite(final, barfig[[10]], offset = "+680+230")
final <- image_composite(final, barfig[[11]], offset = "+685+280")
# final <- image_composite(final, magickleg, offset = "+60+250")
png("C:/Users/keyserf/Documents/github/ICM/Figures/map.png", width = 850, height=550, res=96)
print(final)
dev.off()
