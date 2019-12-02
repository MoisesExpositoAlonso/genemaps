#######################################################################################################################################
####  MOTIVATION ####
# The aim is to model environmentally the presence and absence of alleles, and use the fitted models to predict the future change of
# alleles.
# There are several implementations and options, but the most important is that several assumptions of migration can be made:
# no migration - it fits lat long in the model
# intermediate - controled by the values in 3 genetic PCA axis that have been previously modeled environmentally
# universal migration - no controls, thus alleles can be predicted anywhere as long as the environmental change is large enough.

#######################################################################################################################################

##PARSE COMMAND LINE
library(argparser)

cat("\n")
p <- arg_parser("Genome Environment Model (GEM) and prediction of alleles based on random forests")
p <- add_argument(p, "gen", help="The subset of SNPs should be analysed")
p <- add_argument(p, "nameanalysis", help="The name of the analysis")
p <- add_argument(p, "--runtype", help="Random Forest should be of type classification or regression", default="classification")
p <- add_argument(p, "--makeparallel", help="Do you want to parallel:", default=TRUE)
p <- add_argument(p, "--PCAcontrol", help="Limit predictions of alleles using PCA", default=FALSE)
p <- add_argument(p, "--geocontrol", help="Limit predictions of alleles using latitude and longitude", default=FALSE)
p <- add_argument(p, "--limitlayers", help="Predict over limited number of climate layers?", default=TRUE)
argv<-parse_args(p)

# argv <-list(gen="gene",nameanalysis="run2_gene_pca",PCAcontrol=T,limitlayers=T,geocontrol=F,runtype="classification")

cat(paste("\n ----------------------------------------------------- \n"))
cat("\n Interpreting arguments:\n")
genetics=argv$gen
  cat(paste("\n  --> The SNPs subset modeled will be of type: ",genetics,"\n"))
nameanalysis=argv$nameanalysis
  cat(paste("\n  --> The run name is: ",nameanalysis,"\n"))
makeparallel=argv$makeparallel
  cat(paste("\n    Computation in parallel?",makeparallel,"\n"))
PCAcontrol=argv$PCAcontrol
  cat(paste("\n    Control predictions by PCA structure?",PCAcontrol,"\n"))
geocontrol=argv$geocontrol
  cat(paste("\n    Control predictions by geographic limits?",geocontrol,"\n"))
runtype=argv$runtype
  cat(paste("\n    The modeling will be of type: ",runtype,"\n"))
limitlayers=argv$limitlayers
cat(paste("\n    The prediction will be limited in only important layers: ",limitlayers,"\n"))
cat(paste("\n ----------------------------------------------------- \n"))


#######################################################################################################################################
source('loadfunctions.R')

#######################################################################################################################################
# Get accmaster table and accmaster_add with alleles
cat("\n read the genome matrix ... ")

accmaster_add=read_alleles(name=genetics)

snpnames<-names(accmaster_add)[grep("chr",names(accmaster_add))]

# snpnames=snpnames[1:5] ## need to change this, only for profiling

print(paste("working with a total of:",length(snpnames),"SNPs"))


#######################################################################################################################################
### BIOCLIM DATA!

cat("\n load climate database ... ")

EuropeClim_all=load_rbioclim()

layernames=names(EuropeClim_all)

if(limitlayers==T){ layernames=c("bio_2-5m_bil","cc26bi70","cc85bi70","cclgmbi_2-5m") }

#######################################################################################################################################
if(PCAcontrol==T){
cat("\n load current PCA predicted rasters ... ")

load("results/PC1.RObject")
load("results/PC2.RObject")
load("results/PC3.RObject")

}
#######################################################################################################################################
if(geocontrol ==T){
  message("including geographic layers in environmental rasters")
  georaster=make_georaster(EuropeClim_all[[1]][[1]] ) #the europe raster is just to have same dimensions
}

#######################################################################################################################################
### MODEL PRESENT ##
cat("\n start modeling current climate ... \n")
layername=layernames[1] # the first is current climate

EuropeClim=EuropeClim_all[[layername]]
if(PCAcontrol==T){  EuropeClim<-morelayers_2_bioclim(bioclim = EuropeClim,newlayerslist = stack(list(PC1=PC1,PC2=PC2,PC3=PC3)) ) } # IMPORTANT THAT PCA EXTENTS NEED TO BE EQUAL AS THE ENVIRONMENT EXTENTS
if(geocontrol==T){  EuropeClim<-morelayers_2_bioclim(bioclim = EuropeClim,newlayerslist = georaster ) } # IMPORTANT THAT PCA EXTENTS NEED TO BE EQUAL AS THE ENVIRONMENT EXTENTS


modstack<-model_alleles(EuropeClim = EuropeClim,snpnames = snpnames,layername=layername,endingname = paste0("models_",nameanalysis) ,cluster = makeparallel,geo=geocontrol,PCA = PCAcontrol,runtype = runtype, accmaster_add=accmaster_add)

# dum<-model_alleles(EuropeClim = EuropeClim,snpnames = snpnames[1],layername=layername,endingname = paste0("models_","dum") ,cluster = F,geo=geocontrol,PCA = PCAcontrol,runtype = runtype, accmaster_add=accmaster_add)
# dum$chr1_3472775$custom.err.rate

#######################################################################################################################################
#### PREDICT IN THE REST OF THE CLIMATES
cat("\n predict alleles in all climates ... ")

for (layername in layernames){
print(paste("working in layer:", layername))

# EuropeClim<-cropenvironment(mybioclim =bioclim[[layername]],xlim =  c(-10.5,+ 53),ylim=c(32,65),replace = T, addPCA = PCAcontrol)
# print(EuropeClim)
EuropeClim<-EuropeClim_all[[layername]]
if(PCAcontrol==T){  EuropeClim<-morelayers_2_bioclim(bioclim = EuropeClim,newlayerslist = stack(list(PC1=PC1,PC2=PC2,PC3=PC3)) ) }
if(geocontrol==T){  EuropeClim<-morelayers_2_bioclim(bioclim = EuropeClim,newlayerslist = georaster ) } # IMPORTANT THAT PCA EXTENTS NEED TO BE EQUAL AS THE ENVIRONMENT EXTENTS

## predictions

allelstack<-predict_alleles(EuropeClim,modstack,snpnames,layername=layername,endingname=paste0("allelstack_",nameanalysis) ,cluster= makeparallel)

}
