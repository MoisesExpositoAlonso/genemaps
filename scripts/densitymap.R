################################################################################
# Script to generate a denstiy map from geographic points, example niche model, 
# and example high quality raster plot
# author: moisesexpositoalonso@gmail.com
################################################################################

library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(raster)
library(devtools)
library(RColorBrewer)
devtools::install_github('MoisesExpositoAlonso/moiR')
library(moiR)
load_all('.')

################################################################################
# Several techniques to modify raster resolutions, calculate densities, and 
# apply masks

#### load bioclim ####
world<-raster::getData('worldclim', var='bio', res=10)[[1]]
world<-cropenvironment(world)
plot(world)

#### get the arabidopsis data ####
ara<-read.table('data/Arabidopsis_thaliana_world_accessions_list.tsv',header = T,fill=T)
ara$longitude <- fn(ara$longitude)
ara$latitude <- fn(ara$latitude)


#### estimate density ####
coords = cbind(ara$longitude,ara$latitude)
coords<-na.omit(coords)
sp = SpatialPoints(coords)

worlddiluted<-aggregate(x = world,fact=10)
den <- rasterize(sp, worlddiluted, fun='count')
plot(den)
denrescale<-resample(den, world, method='bilinear') 
plot(denrescale)

# Generate a mask below certain density
denrescale[denrescale <= 10]<-NA # generate a maks density
worldmask<-world
worldmask[is.na(denrescale)] <-  NA # mask in original environment
plot(worldmask)

################################################################################
# Example niche model
library(randomForest)

world<-raster::getData('worldclim', var='bio', res=10)
world<-cropenvironment(world)

ara<-read.table('data/Arabidopsis_thaliana_world_accessions_list.tsv',header = T,fill=T)
ara$longitude <- fn(ara$longitude)
ara$latitude <- fn(ara$latitude)
ara <- dplyr::filter(ara, latitude > 32, latitude <65, longitude> -10, longitude<53)
dim(ara)

 # add the climate
extractedenv<- lapply(names(world), FUN=function(x){
  mytmp2<-sapply( names(world[[x]] ) , FUN= function(y) {
  mytmp<-data.frame( extract(world[[x]][[y]],ara[,c("longitude","latitude")] ) )
  names(mytmp)<-paste0(x)
  return(mytmp)} )
  mytmp2<- do.call(cbind,mytmp2)
  return(mytmp2)
  } ) %>% do.call(cbind,.)
head(extracted)
colnames(extracted)<-paste0("bio",1:19)

# dummy prediction (environmental niche model of who is RegMap accession vs 1001)
# This prediction is obviously not so "biological", it just would tell us what 
# environments Regmap is biased compared to 1001g. You want to subsitute that 
# response vairable for something like microbial load
master<-data.frame(regmap=(((ara$RegMap) %>% as.character() ) ==TRUE)*1, 
                   # to make TRUE FALSE -> 1/0
                   extracted)

mod<-randomForest(data=master,  regmap ~ .)

pred<-raster::predict(world,mod)
plot(pred)


################################################################################
# Ways to plot any raster a bit fancier than the default raster "plot"
example=pred
# Fancy plot removing geographic edges and nicer palette 
moiR::envirplot(example) 

# Fancy plot with ggplot
require(ggthemes)

test_spdf <- as(example, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")


myfancierplot<-ggplot() +
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8)+
  scale_fill_gradient2("Annual precipitation change",  low='black',mid='lightgrey',high='blue')+
  # geom_point()+ # you can add dots on top of it
  coord_equal() +
  theme_map()+  theme(legend.position = "bottom")
myfancierplot

