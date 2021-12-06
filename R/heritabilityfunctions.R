
#### START FUNCTIONS ####
# RsqGLM <- function(obs = NULL, pred = NULL, model = NULL) {
  # version 1.2 (3 Jan 2015)
#
#   model.provided <- ifelse(is.null(model), FALSE, TRUE)
#
#   if (model.provided) {
#     if (!("glm" %in% class(model))) stop ("'model' must be of class 'glm'.")
#     if (!is.null(pred)) message("Argument 'pred' ignored in favour of 'model'.")
#     if (!is.null(obs)) message("Argument 'obs' ignored in favour of 'model'.")
#     obs <- model$y
#     pred <- model$fitted.values
#
#   } else { # if model not provided
#     if (is.null(obs) | is.null(pred)) stop ("You must provide either 'obs' and 'pred', or a 'model' object of class 'glm'")
#     if (length(obs) != length(pred)) stop ("'obs' and 'pred' must be of the same length (and in the same order).")
#     if (!(obs %in% c(0, 1)) | pred < 0 | pred > 1) stop ("Sorry, 'obs' and 'pred' options currently only implemented for binomial GLMs (binary response variable with values 0 or 1) with logit link.")
#     logit <- log(pred / (1 - pred))
#     model <- glm(obs ~ logit, family = "binomial")
#   }
#
#   null.mod <- glm(obs ~ 1, family = family(model))
#   loglike.M <- as.numeric(logLik(model))
#   loglike.0 <- as.numeric(logLik(null.mod))
#   N <- length(obs)
#
#   # based on Nagelkerke 1991:
#   CoxSnell <- 1 - exp(-(2 / N) * (loglike.M - loglike.0))
#   Nagelkerke <- CoxSnell / (1 - exp((2 * N ^ (-1)) * loglike.0))
#
#   # based on Allison 2014:
#   McFadden <- 1 - (loglike.M / loglike.0)
#   Tjur <- mean(pred[obs == 1]) - mean(pred[obs == 0])
#   sqPearson <- cor(obs, pred) ^ 2
#
#   return(list(CoxSnell = CoxSnell, Nagelkerke = Nagelkerke, McFadden = McFadden, Tjur = Tjur, sqPearson = sqPearson))
# }



extractheritability<-function(h2fit){
  mean(h2fit)
  plot(h2fit)
  autocorr(h2fit)
  formateo<-(c(mean(h2fit),HPDinterval(h2fit)[1],HPDinterval(h2fit)[2]))
  formateo<-round(formateo,digits=4)
  toexport<-(paste(formateo[1]," (",formateo[2]," , ",formateo[3],")",sep=""))
  return(toexport)
}

h2modelperexp<-function(dataexp,chainlength=10000){
  chainlength=chainlength
  prior<-  list(R=list(V=diag(2)*(0.002/1.002),nu=1.002),
                G=list(G1=list(V=diag(2)*(0.002/1.002),nu=1.002)) )
  model<-MCMCglmm(data=dataexp,verbose=F,
                           fixed = cbind(scale(Flowering),scale(Germination)) ~ trait  - 1,
                           random= ~ us(trait):Genotype ,
                           family=c("gaussian","gaussian"),
                           rcov=~us(trait):units,
                          prior=prior,
                           ginverse=list(Genotype=Ai),
                           pl=T,
                           pr=T,nitt=chainlength,burnin= (chainlength*0.1) )

  herit1<-model$VCV[,"Flowering:Flowering.Genotype"]/
    (model$VCV[,"Flowering:Flowering.Genotype"]+model$VCV[,"Flowering:Flowering.units"] )
  herit2<-model$VCV[,"Germination:Germination.Genotype"]/
    (model$VCV[,"Germination:Germination.Genotype"]+model$VCV[,"Germination:Germination.units"] )

  h2_1<-extractheritability(herit1)
  h2_2<-extractheritability(herit2)
  corr.gen<-model$VCV[,"Flowering:Germination.Genotype"]/
    sqrt(model$VCV[,"Flowering:Flowering.Genotype"]*model$VCV[,"Germination:Germination.Genotype"])
  corr_12<-extractheritability( corr.gen)
  listtoreturn=list(h2_1,corr_12,h2_2)
  return(listtoreturn)
}

h2modelperexpNUMERIC<-function(dataexp,chainlength=10000){
  chainlength=chainlength
  prior<-  list(R=list(V=diag(2)*(0.002/1.002),nu=1.002),
                G=list(G1=list(V=diag(2)*(0.002/1.002),nu=1.002)) )
  model<-MCMCglmm(data=dataexp,verbose=F,
                  fixed = cbind(scale(Flowering),scale(Germination)) ~ trait  - 1,
                  random= ~ us(trait):Genotype ,
                  family=c("gaussian","gaussian"),
                  rcov=~us(trait):units,
                  prior=prior,
                  ginverse=list(Genotype=Ai),
                  pl=T,
                  pr=T,nitt=chainlength,burnin= (chainlength*0.1) )

  herit1<-model$VCV[,"Flowering:Flowering.Genotype"]/
    (model$VCV[,"Flowering:Flowering.Genotype"]+model$VCV[,"Flowering:Flowering.units"] )
  herit2<-model$VCV[,"Germination:Germination.Genotype"]/
    (model$VCV[,"Germination:Germination.Genotype"]+model$VCV[,"Germination:Germination.units"] )

  h2_1<-posterior.mode(herit1)
  h2_2<-posterior.mode(herit2)
  corr.gen<-model$VCV[,"Flowering:Germination.Genotype"]/
    sqrt(model$VCV[,"Flowering:Flowering.Genotype"]*model$VCV[,"Germination:Germination.Genotype"])
  corr_12<-posterior.mode( corr.gen)
  listtoreturn=list(h2_1,corr_12,h2_2)
  return(listtoreturn)
}




multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

summary_env_pot<-function(datefile=dates_flowering, climafile=grasnclimate,startexp=startexp){

# med_t: median day temp
# mean_max_t: mean maximum temperature
# mean_min_t: mean minimum temperature
# max_max_t: maximum maximum temperature
#min_min_t: minimum minimum temperature
#cumul_t: accumulated degrees temperature
#range_t: mean diurnal range

# mean_p: mean day precipitation
# cumul_p: accumulated precipitation
#var_p: variance in precipitation

# last_15_m_t: 15 days before flowering, mean, sd and cumulatie
# last_15_sd_t: 15 days before flowering, mean, sd and cumulatie
# last_15_cumul_t: 15 days before flowering, mean, sd and cumulatie
# last_15_m_p: 15 days precipitation mean
# last_15_sd_p: 15 days precipitation sd
# last_15_cumul_p: 15 days precipitation cumulative

# after_15_m_t: 15 days before flowering, mean, sd and cumulatie
# after_15_sd_t: 15 days before flowering, mean, sd and cumulatie
# after_15_cumul_t: 15 days before flowering, mean, sd and cumulatie
# after_15_m_p: 15 days precipitation mean
# after_15_sd_p: 15 days precipitation sd
# after_15_cumul_p: 15 days precipitation cumulative

as.Date(startexp[1,1])+1

datefile$med_t =NA
datefile$mean_max_t =NA
datefile$mean_min_t =NA
datefile$max_max_t =NA
datefile$min_min_t =NA
datefile$cumul_t =NA
datefile$range_t =NA

datefile$mean_p=NA
datefile$cumul_p=NA
datefile$var_p = NA

datefile$last_15_m_t =NA
datefile$last_15_sd_t=NA
datefile$last_15_cumul_t =NA
datefile$last_15_m_p=NA
datefile$last_15_sd_p=NA
datefile$last_15_cumul_p=NA

datefile$after_15_m_t=NA
datefile$after_15_sd_t=NA
datefile$after_15_cumul_t=NA
datefile$after_15_m_p=NA
datefile$after_15_sd_p=NA
datefile$after_15_cumul_p=NA


climafile$DATE <- as.Date( climafile$DATE, '%m/%d/%Y')
datefile$DATE<-as.Date( datefile$DATE, '%m/%d/%Y')

explocations<-c("GRA","GRA","GRA","GRA","GRA","GRA","SN","SN","SN")


for(r in 1:dim(datefile)[1] ){
row=datefile[r,]

enddate<- as.Date( row$DATE, '%m/%d/%Y')
startdate<- as.character(startexp[row$EXP,1])
expsite<-explocations[row$EXP]

subclim<-subset(climafile, climafile$SITE==expsite & climafile$DATE >startdate & climafile$DATE< enddate)
# head(subclim)
# tail(subclim)

datefile$med_t [r]  = mean( ( subclim$TMIN + subclim$TMAX) /2 ,na.rm=T)
datefile$mean_max_t[r] = mean(subclim$TMAX)
datefile$mean_min_t [r]=mean(subclim$TMIN,na.rm=T)
datefile$max_max_t [r]= max(subclim$TMAX,na.rm=T)
datefile$min_min_t [r]=min(subclim$TMIN,na.rm=T)
datefile$cumul_t [r]= sum(( subclim$TMIN + subclim$TMAX) /2 ,na.rm=T)
datefile$range_t [r]= mean( ( subclim$TMAX - subclim$TMIN) ,na.rm=T)

datefile$mean_p [r] =mean(subclim$RAIN,na.rm=T)
datefile$cumul_p [r]=sum(subclim$RAIN,na.rm=T)
datefile$var_p [r] = sd(subclim$RAIN,na.rm=T)

subclim<-subset(climafile, climafile$SITE==expsite & climafile$DATE >enddate-15 & climafile$DATE< enddate)

datefile$last_15_m_t [r] = mean((subclim$TMAX- subclim$TMIN) /2 )
datefile$last_15_sd_t [r]=sd((subclim$TMAX- subclim$TMIN) /2 )
datefile$last_15_cumul_t [r] =sum((subclim$TMAX- subclim$TMIN) /2 )
datefile$last_15_m_p [r]= mean(subclim$RAIN)
datefile$last_15_sd_p [r]=sd(subclim$RAIN)
datefile$last_15_cumul_p [r]=sum(subclim$RAIN)

subclim<-subset(climafile, climafile$SITE==expsite & climafile$DATE >enddate & climafile$DATE< enddate +15)

datefile$after_15_m_t [r] = mean((subclim$TMAX- subclim$TMIN) /2 )
datefile$after_15_sd_t [r]=sd((subclim$TMAX- subclim$TMIN) /2 )
datefile$after_15_cumul_t [r] =sum((subclim$TMAX- subclim$TMIN) /2 )
datefile$after_15_m_p [r]= mean(subclim$RAIN)
datefile$after_15_sd_p [r]=sd(subclim$RAIN)
datefile$after_15_cumul_p [r]=sum(subclim$RAIN)

}
return(datefile)

}

myglmmMCMC<-function(data=data_final,response,fixpredictors=c("FRI * FLC"),randompredictors) {
  library(MCMCglmm)
#   fix<-paste(fixpredictors,collapse=" * ")
#   rand<-paste(randompredictors,collapse=" * ")
  resp=paste("scale(",response,")")
  formfix=formula(paste (resp , " ~ ", fixpredictors ) )

  formran=formula(paste( " ~ ", randompredictors, " + Genotype") )

lmfull<-MCMCglmm(data=data_final,
      formfix , random= formran  ,
#     pl=T,   pr=T,
      nitt=10000, burnin=1000,
      ginverse=list(Genotype=Ai))

  summary(lmfull)

}

# myglmmMCMC.f<-function(data,formfix,formran) {
#   library(MCMCglmm)
# lmfull<-MCMCglmm(data=data,
#       formfix , random= formran  ,
#       pl=T,   pr=T,
#       nitt=50000, burnin=5000,
#       ginverse=list(Genotype=Ai))
#
#   tr<-summaryglmmMCMC(lmfull)
#   print(tr)
#   print("Deviance Information Criterion:")
#   print(lmfull$DIC)
#   return(list(summary=tr,DIC=lmfull$DIC))
#
# }

designglmmMCMC<-function(data,moddesign) {
  library(MCMCglmm)
#   fix<-paste(fixpredictors,collapse=" * ")
#   rand<-paste(randompredictors,collapse=" * ")
  # resp=paste("scale(",response,")")
  # formfix=formula(paste (resp , " ~ ", fixpredictors ) )
  #
  # formran=formula(paste( " ~ ", randompredictors, " + Genotype") )

nummodels<-dim(moddesign)[1]

for ( m in 1:nummodels){
mod=moddesign[m,]
formfix<-formula( as.character(moddesign[m,]$V1))
formran<-formula( as.character(moddesign[m,]$V2) )

print(paste(formfix,formran))

lmfull<-myglmmMCMC.f(data,formfix ,formran )

} #end loop
} #end function

extractheritability<-function(h2fit){
  # mean(h2fit)
  # plot(h2fit)
  # autocorr(h2fit)
  # formateo<-(c(mean(h2fit),HPDinterval(h2fit)[1],HPDinterval(h2fit)[2]))
  # plot(proph2[,1],type="l")
  # plot(density(proph2[,1]))
  plot(h2fit,type="l")

  formateo<-c(mean(h2fit),quantile(h2fit,p=0.025),quantile(h2fit,p=0.975))
  formateo<-round(formateo,digits=4)
  toexport<-(paste(formateo[1]," (",formateo[2]," , ",formateo[3],")",sep=""))
  return(toexport)
 }

exh2<-function(modelrandom=model$VCV){

   allrandom<-apply(modelrandom,1,FUN=sum)
   proph2<-apply(modelrandom,2,FUN=function(x){x / allrandom})
   apply (proph2,2, FUN=extractheritability)

 }

summaryglmmMCMC<-function(model){
  fixres<-summary(model)$solutions
      term<-rep("fixed",dim(fixres)[1])
      fixres<-cbind(fixres,term )
  randres<-summary(model)$Gcovariances
  errorres<-summary(model)$Rcovariances

      propVar<-randres[,"post.mean"]/sum(randres[,"post.mean"] + errorres[,"post.mean"] )
      propVarerror<-errorres[,"post.mean"]/sum(randres[,"post.mean"] + errorres[,"post.mean"] )

term<-rep("random",length(propVar))
randres<-cbind(randres,propVar, term)
term<-"residual"
errorres<-cbind(errorres,propVarerror,term)

toreturn<-rbind(fixres,randres,errorres)
return(toreturn)

}


myglmmMCMC.f<-function(data,formfix,formran,uncorrelated=F) {
  library(MCMCglmm)

if(uncorrelated==F){
lmfull<-MCMCglmm(data=data,
      formfix , random= formran  ,
      pl=T,   pr=T,
      nitt=50000, burnin=5000,
      ginverse=list(Genotype=Ai))
}
if(uncorrelated==T){
lmfull<-MCMCglmm(data=data,
      formfix , random= formran  ,
      pl=T,   pr=T,
      nitt=50000, burnin=5000)
}

  tr<-summaryglmmMCMC(lmfull)
  print(tr)
  print("Deviance Information Criterion:")
  print(lmfull$DIC)
  return(list(summary=tr,DIC=lmfull$DIC))

}

myglmmMCMC.fmulti<-function(data,formfix,formran,uncorrelated=F) {
  library(MCMCglmm)

if(uncorrelated==F){
lmfull<-MCMCglmm(data=data,
      formfix , random= formran  ,
      pl=T,   pr=T,
      nitt=50000, burnin=5000,
      ginverse=list(Genotype=Ai),
      family=c("gaussian","gaussian"))
}
if(uncorrelated==T){
lmfull<-MCMCglmm(data=data,
      formfix , random= formran  ,
      pl=T,   pr=T,
      nitt=50000, burnin=5000)
}

  tr<-summaryglmmMCMC(lmfull)
  print(tr)
  print("Deviance Information Criterion:")
  print(lmfull$DIC)
  return(list(summary=tr,DIC=lmfull$DIC))

}
