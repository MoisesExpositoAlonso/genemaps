################################################################################
## Compare present and past ifnerences of degree of adaptaiton
################################################################################

devtools::load_all('.')
p<-argumentsgws(clim='MP8550')

####************************************************************************####
#### Check for consistency first ####
presadapt<- readRDS("maps/degreeadapt_gwaclim--spager.rda")
presadaptsix<- readRDS("maps/degreeadapt_gwaclim--sixenv.rda")

cor.test(presadapt$degreeadapt$meanabssel,presadaptsix$degreeadapt$meanabssel)
####************************************************************************####
 #### Change present to future ####

# futadapt<-readRDS(file = list.files(path='maps' ,pattern=paste0('degreeadapt_', predi,"_",clim,"_"),full.names=T  ))
# presadapt<-readRDS(file = list.files(path='maps' ,pattern=paste0('degreeadapt_', predi,"_",'pres',"_") , full.names=T)[1] ) # I know that the first is the one with more SNPs
presadapt<- readRDS(paste0("maps/degreeadapt_gwaclim--",myrfname,".rda"))
futadapt<- readRDS(paste0("maps/degreeadapt_gwaclim-MP8550-",myrfname,".rda"))

all(dim(presadapt) == dim(futadapt))

# degreeadapt<-presadapt$

mydiff<- futadapt$degreeadapt[,-c(1:3)]- presadapt$degreeadapt[,-c(1:3)]
mydiff<-cbind(futadapt$degreeadapt[,c(1:3)], mydiff)


# library(moiR)
mymap<-ggplot_world_map(xlim=c( -15, 100), ylim=c(25,65),projection = "perspective")
# mymap


pdf(file=paste0("figs/maps/changeadapt-",myrfname,".pdf"), height =  4,width = 12)

for(i in colnames(mydiff)[-c(1:3)]){
  message(i)
  p<-mymap+
            geom_point(aes(y=fn(mydiff[,'latitude']), x=fn(mydiff[,'longitude']),
                           color=mydiff[,i]),size=0.6)+
             # scale_color_gradientn(colours=make.pallete.contrast(mypalettes('moispectral'),contrast=1,raw = T))+
             scale_color_gradientn(colours=brewer.pal(5,"RdYlGn"))+
            ggtitle(i)
  ggplot2:::plot.ggplot(p)
}

dev.off()

####************************************************************************####
#### direct quantification ####


####  get the genome matrix
sg<-read_n_subset_genome(selectedsnps = coincidingSNPs)[euroacc,]
sg <- (sg- 0.5) *2

loadlocchange<- futadapt$raw  - presadapt$raw

loadlocchange[loadlocchange>=0.05]<-1
loadlocchange[loadlocchange<= -0.05]<- -1
loadlocchange[loadlocchange<0.05 & loadlocchange> -0.05 ]<-0

loadlocchange %>% fn %>% summary
mul<-loadlocchange * sg

res<-apply(mul,1, sum)

ca<-data.frame(latitude=acceuro$latitude,
                        longitude=acceuro$longitude,
                        kgroup=acceuro$kgroup,
                       changepoly=res
)


####  plot
present<-stack('dataint/newclim2.grd')
fut<-stack('dataint/MP8550newclim2.gri')

habitatdifference<-fut$bio12 - present$bio12
broadextremes<-broadeuroextents()
habitatdifference<-cropenvironment(habitatdifference, xlim=broadextremes$xlim, ylim=broadextremes$ylim +c(5,0))
# habitatdifference<-fut$gs_sum - present$gs_sum

example<-habitatdifference
test_spdf <- as(example, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

test_df<-test_df %>%
                  mutate(value= logistic(value,k=0.1) )


require(ggthemes)
png(paste0("figs/maps/mapsel-futurechange",'sumchangeadaptation' ,".png"),
          height=8,width=12,units="in", res=800)

myplot3<-ggplot() +
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8)+
  geom_point(aes(y=fn(ca[,'latitude']), x=fn(ca[,'longitude'])),color="black",size=1.5)+
  geom_point(aes(y=fn(ca[,'latitude']), x=fn(ca[,'longitude'])),color="white",size=1.1)+
  geom_point(aes(y=fn(ca[,'latitude']), x=fn(ca[,'longitude']),color=ca[,"changepoly"] ),size=1)+
  # scale_color_gradientn("Change in absolute selection" , colours=make.pallete.contrast(rev(brewer.pal(name='RdGy',5)) ,raw=T))+
  scale_color_gradient2("# local alleles changing selection" , low = 'red',high = 'green', mid='white',midpoint = 0)+
  # scale_fill_gradientn("Annual precipitation change",  colours=c('black', rev(brewer.pal(name='Greys',4)),brewer.pal(name='Greys',4)[1]))  +
  scale_fill_gradientn("Annual precipitation change",  colours=c('black',rev(brewer.pal(name='Greys',4)) ) )  +

  # scale_fill_gradient2("Annual precipitation change",  colours=c(rev(brewer.pal(name='Greys',4)[-1]),brewer.pal(name='PuBu',3)[-c(1:2)] ) ) +
  # scale_fill_gradient2("Annual precipitation change",  low='black',mid='lightgrey',high='lightgrey', midpoint=0)+
  # scale_fill_gradient2("Annual precipitation change",  low='black',mid='lightgrey',high='lightgrey', midpoint=0.5)+
  coord_equal() +
  theme_map()+  theme(legend.position = "bottom")
myplot3
d()



#####**********************************************************************#####
#### Wave of selection change ####

write.pdf(paste0("figs/change-wave",'sumchangeadaptation' ,".pdf"),heightmm = 90, widthmm = 180)

ggplot(ca) +
  geom_point(aes(x=latitude,y=changepoly), color='black', size=2)+
  geom_point(aes(x=latitude,y=changepoly), color='white',size=1.6)+
  geom_point(aes(x=latitude,y=changepoly,color=changepoly),size=1.4)+
                scale_color_gradient2("Change in absolute selection" , low = 'red',high = 'green', mid='white',midpoint = 0) +
  stat_smooth(aes(x=latitude,y=changepoly) ,se=F,col=transparent('black'))



dev.off()


####************************************************************************####
#### PLot against climate ####
accpres<-readRDS('dataint/acc515envextended.rda')
accfut<-readRDS('dataint/MP8550acc515envextended.rda')

diffloc<-accfut -accpres
diffloc<-diffloc[euroacc,]

summary(diffloc[,'bio1']/10)

require(randomForest)

toforest<-cbind(mydiff[,"meanselection"], diffloc)
toforest<-na.omit(toforest)
mod<-randomForest(toforest[,1], x= toforest[,-1],na.rm=T)
myimportance<-importance(mod) %>% data.frame()
myimportance$vars <- row.names(myimportance)

arrange(myimportance, IncNodePurity)

cormod<-cor(toforest)
cormod[1,] %>% sort


write.pdf(paste0("figs/change-",'polygenicscore' ,".pdf"))

(ggdotscolor(x=fn(mydiff[,'weightedsum']),
           # x=fn(diffloc[,'bio12']),
           y=fn(presadapt$degreeadapt[,'latitude']),
           varcol= fn(presadapt$degreeadapt[,"weightedsum"]),
           mycolors=rev(mypalettes('moispectral') ),
           ylab="Change in annual precipitation",
           xlab="Change in polygenic score") #+
  # geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0))
  ) %>%
  addggregression(docorrelation=T) %>%
  nolegendgg()
d()


plotbaseraster()
points( y=fn(presadapt$degreeadapt[,'latitude']),
         x=fn(presadapt$degreeadapt[,'longitude']))
d()

#####**********************************************************************#####
#### Change in precipiattion and adaptability ####


present<-stack('dataint/newclim2.grd')
fut<-stack('dataint/MP8550newclim2.gri')

habitatdifference<-fut$bio12 - present$bio12
broadextremes<-broadeuroextents()
habitatdifference<-cropenvironment(habitatdifference, xlim=broadextremes$xlim, ylim=broadextremes$ylim +c(5,0))
# habitatdifference<-fut$gs_sum - present$gs_sum

example<-habitatdifference
test_spdf <- as(example, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

require(ggthemes)
png(paste0("figs/maps/mapsel-futurechange",'polygenicscore' ,".png"),
          height=8,width=12,units="in", res=800)

myplot3<-ggplot() +
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8)+
  geom_point(aes(y=fn(mydiff[,'latitude']), x=fn(mydiff[,'longitude'])),color="black",size=1.2)+
  geom_point(aes(y=fn(mydiff[,'latitude']), x=fn(mydiff[,'longitude'])),color="white",size=1.1)+
  geom_point(aes(y=fn(mydiff[,'latitude']), x=fn(mydiff[,'longitude']),color=mydiff[,"weightedsum"]),size=1)+
  scale_color_gradientn("Polygenic score change" , colours=make.pallete.contrast(brewer.pal(name='RdGy',5) ,raw=T))+
  # scale_fill_gradientn("Annual precipitation change",  colours=c(rev(brewer.pal(name='Greys',4)[-1]),brewer.pal(name='PuBu',3)[-c(1:2)] ) ) +
  # scale_fill_gradient2("Annual precipitation change",  colours=c(rev(brewer.pal(name='Greys',4)[-1]),brewer.pal(name='PuBu',3)[-c(1:2)] ) ) +
  scale_fill_gradient2("Annual precipitation change",  low='black',mid='lightgrey',high='blue')+
  coord_equal() +
  theme_map()+  theme(legend.position = "bottom")
myplot3
d()

