

################################################################################

extinction.gridplot<-function(sim.d,name=deparse(substitute(sim.d))){

  gridplot<-ggplot(sim.d)+
    # geom_jitter(aes(x=d,y=X,color=evotype,size=No,alpha=No,shape=evotype)) +
    geom_jitter(aes(x=d,y=X,color=evotype,size=No,alpha=evotype)) +
    labs(y= 'Ecotype extinction',x= 'Individual mortality')+
    ggtitle(name)+
    ylim(c(0.5,1))+
    scale_size('N0',range = c(0.1, 3),
               labels =c(100,1000,10000,100000,1000000),breaks = c(100,1000,10000,100000,1000000))+
    scale_color_manual("",values=extinctioncolors()) +
    # scale_shape_manual("",values=c(extinction=19,rescue=19,growth=1))
    # scale_alpha('N0',range = c(0.8,0.1))+
    scale_alpha_manual("",values=c(extinction=0.8,rescue=0.3,growth=0.1))

  return(gridplot)

}

# The bidirectional plot
# p<-ggplot(rx_rescue,aes(x=Go,y=No , z= mind,colour=mind ))+ geom_jitter(size=2) +
#   xlab('Genotypes')+ylab('Starting populations size')
# p<- p+  scale_colour_gradientn('minimum drift',colours=moiR::mypalettes('jetb') )
# save_plot(file='figs/popX_grid.pdf',plot = p,base_height = 5,base_width =8)

# # The unidirectional plot
# Go<-ggdotscolor(rx_rescue$mind , x=rx_rescue$No,xlab = 'Starting population size',ylab='Minimum drift for evolutionary rescue',
#                  varcol = rx_rescue$mind,mycolors = moiR::mypalettes('jetb') ) %>% addggregression(se=FALSE) %>% nolegendgg()
# No<-ggdotscolor(rx_rescue$mind , x=rx_rescue$Go,xlab = 'Genotypes',ylab='Minimum drift for evolutionary rescue',
#                  varcol = rx_rescue$mind,mycolors = moiR::mypalettes('jetb') )%>% addggregression(se=FALSE)%>% nolegendgg()
# panel=plot_grid(Go,No,ncol=2)
# panel
# save_plot(file='figs/popX_regression.pdf',plot = panel ,base_height = 4.5,base_width =9 )


extinctioncolors<-function(){
c(extinction='#b2182b',growth='#2166ac',rescue='black')
}

findTmin<-function(cleansim){
  apply(cleansim,1,function(row){
    head(which(row == min(row)),n=1)
  })
}

findNmin<-function(cleansim){
  apply(cleansim,1, function(row){
    head(row[which(row == min(row))],n=1)
  })
}

make.datadrift<-function(cleansim){
  tmp=data.frame(d=fn(rownames(cleansim$N)), No=cleansim$N[,1],Go=cleansim$G[,1] , Nmin= findNmin(cleansim$N) ,tmin=findTmin(cleansim$N),simID=1:nrow(cleansim$N))

   tmp$evotype<-'extinction'
  tmp$evotype[tmp$tmin ==1]<-"growth"
  tmp$evotype[tmp$tmin !=1 & tmp$Nmin !=0]<-"rescue"

  tmp$X<-1-(tmp$Go / max(tmp$Go,na.rm = TRUE))

  tmp<-filter(tmp, No!=0 | Go!=0)

  return(tmp)
}

extinction.plink.freq<-function(cluster, name='extinctionfreq', plink='../gemma/515g'){


  write.table(quote = FALSE,row.names = FALSE,col.names = FALSE,sep='\t',
            cluster,file = 'tmpfreq.txt')

  cmd=paste('plink --bfile',plink,' --freq --within tmpfreq.txt --out ', name )
  message('running plink command: ',cmd)
  # system(cmd)
}

extinction.genplink.gwa<-function(cluster,
                               out='extinctionfreq',
                               path='../gemma',
                               famfile='515g.fam',
                               bimfile='515g.bim',
                               bedfile='515g.bed',
                               cleardir=F
                               ){

  # load Fam
  data(fam)

  # merge
  myplink<-merge(fam, cluster[,c('id','fate01')] , by.x='V1',by.y='id' ,all.y = T)

  # output name
  fampath=file.path(path, out)
  famfile=file.path(fampath,'/515g.fam')

  # create directory
  catchexit=system(paste('mkdir',fampath))

  if(catchexit !=0){
  message("problem creating directory!")
    if(cleardir==TRUE){
      message('removing directory')
      clear_gemma_dir(out)
      system(paste('mkdir',fampath))
    }else{
      stop('provide permission to remove directory with cleardir=TRUE')
    }
  }
  # hard link the genome matrix files for gemma
  system(paste('ln', file.path(path, bedfile),fampath))
  system(paste('ln', file.path(path, bimfile),fampath))


  # out
  write.table(myplink,file=famfile,col.names=F,row.names=F,quote=F,sep=" ")
  print(head(myplink))

}

ids.extinction.rescue<-function(sel, sim){

  res=unlist(lapply(sel$simID.x, function(i) (rownames(sim[[i]]))))
  ext=unlist(lapply(sel$simID.y, function(i) (rownames(sim[[i]]))))

  pop=c(
    rep('rescue',length(res)),
    rep('extinction',length(ext))
        )

  tmp<-data.frame(id=c(res,ext), pop=pop)

  tmp$fate01 <- ifelse(tmp$pop=='rescue',1,0)

  return(tmp)
}

select.extinction.rescue<-function(sim.d, sim){

# select those evo types wanted
selX<-filter(sim.d,evotype=='extinction') %>%
      mutate(X=round(X,digits = 2),
             No=ceiling(No/100)*100) %>%
      mutate(conditions=paste(No,X,d,sep='_'))
selER<-filter(sim.d,evotype=='rescue')%>%
      mutate(X=round(X,digits = 2),
             No=ceiling(No/100)*100) %>%
      mutate(conditions=paste(No,X,d,sep='_'))
# Find
selected<-merge(selER, selX, by='conditions' )

return(selected)
}

################################################################################

clean.popXsim<-function(sim){

Nt=lapply(sim,function(inds){
  N= apply(inds,2,function(i) sum(i,na.rm=TRUE) )
  colnames(N) = names(inds)
  return(N)
}) %>% do.call(rbind,.)

Gt=lapply(sim,function(inds){
  N= apply(inds,2,function(i) length(i[i!=0 & !is.na(i)] ) )
  colnames(N) = names(inds)
  return(N)
}) %>% do.call(rbind,.)



return(list(N=Nt,G=Gt))
}

plot.popXsim<-function(cleansim, name='',smooth=FALSE,ylim=cleansim[1,1]){
  t = ncol(cleansim)
  reps=nrow(cleansim)

  plot(ylab= "Total population (N)", xlab='Generations', y=1,x=1,cex=0,ylim= c(0,2*ylim), xlim=c(1,t),main=name)

  if(smooth==FALSE){
    for(i in 1:reps){
      if(cleansim[i,t] ==0){
        colo='darkred'
      }else if(cleansim[i,1] <cleansim[i,2] ){
        colo='navy'
      }else{
        colo='black'
      }

      lines(y=cleansim[i,] ,x= 1:t, col=transparent(colo))
    }
  }else{
    require(splines)
    for(i in 1:reps){
      # stmp<-smooth.spline(y=cleansim[i,],x= 1:t, spar=0.1)
      # stmp<-predict( lm(cleansim[i,] ~ bs(1:t, df = 2)))
      # stmp<-predict( loess(cleansim[i,] ~ c(1:t) ) )
      # stmp<-predict( loess(cleansim[i,] ~ c(1:t) + cleansim[1,1] ) )
      # stmp<-predict( lm(cleansim[i,] ~ bs(c(1:t)) + cleansim[1,1] ) )

      # lines(y=stmp,x=1:t, col=transparent('black'))
    }
  }
}


################################################################################

.extinctionprob<-function(ids=ecotypes,type='-latitude',dat=dry){

  Xprob= if(substring(type,1,1)=='-'){
              1-normalizeC( fn(dat[,substring(type,2,nchar(type))]) )
       }else{normalizeC(fn(dat[,type])) }
  Xprob=Xprob+0.000001

  return(Xprob)
}

.tellseedpop<-function(ids,N){

}

.tellseeds<-function(ids,code='mlp'){
  d=table( dplyr::filter(fieldfilter(field.c,code),id %in%ids)$id ) *30
  return(d)
}

.tellfitness<-function(ids,code='mlp'){
  d=dry[,paste0('Fitness_',code)]
  return(d)
}

################################################################################

# wrap.extinction.populationsimC<-function(ecotypes,t=20,dseq=c(0,seq(0.0001,0.0005,by=0.0001),seq(0.001,0.1,by=0.001)), code='mlp', name=NULL,
#                               ...){
# library(splines)
#
# data(dry)
# data(field.c)
#   sim=lapply(dseq, function(d) {
#       inds<-populationsimC(types= ecotypes,
#                            individuals= table(fieldfilter(field.c, code, id %in% ecotypes)$id) *30,
#                            fitness= dry[,paste0("Fitness_",code)],
#                            t=t,d = d)
#       N=apply(inds,2,function(i) sum(i,na.rm = TRUE))
#       })
#
#     sim= do.call(rbind,args = sim)
#
#     plot(ylab= "Total population (N)", xlab='Generations', y=1,x=1,cex=0,ylim= c(0,2*sim[1,1]), xlim=c(1,t),main=name)
#
#     for(i in 1:nrow(sim)){
#       lines(y=sim[i,] ,x= 1:t, col=transparent('black'))
#
#     }
#
#     returnlist=list(sim=sim)
#
#         return(returnlist)
# }

################################################################################
# l<-length

# extinct<-function(pop, frac=NULL, weights=NULL){
#
#   if(is.null(frac)) frac=0.5
#
#   if(is.null(weights)){
#     weights= rep( 1/length(pop), length(pop))
#   }
#
#   survivors<-sample(x=pop,
#          size = ceiling(frac*length(pop)),
#          prob =weights,
#          replace = FALSE)
#
#   return(survivors)
#
# }
#
#
# extinct.replicates<-function(pop,frac,rep=1,...){
#
# lapply(1:rep, function(r){
#
#   survivors<-
#     lapply(frac, FUN = function(f) {
#    survivors=extinct(pop,f,...)
#    return(survivors)
#    })
#
#   return(survivors)
# })
#
#
# }


#
# cleansim=function(sim,name){
# require(dplyr)
#
# lapply(sim,FUN =function(pop){
#
#   subdry<-subset(dry, id %in% pop)
#
#   returnlist<-list()
#   returnlist$numgeno<-517-length(pop)
#   returnlist$diversity=mean(fn(subsetmat(gendist,pop)), na.rm=TRUE)
#   returnlist$fitness= max(subdry$Fitness_mli,na.rm=TRUE)
#   returnlist$survival= max(subdry$Survival_fruit_mli,na.rm=TRUE)
#
#   return(unlist(returnlist))
#
# }) %>% do.call(.,what = rbind) %>% data.frame %>% mutate(name=name)
#
# }
#
# run_extinction_design<-function(designsim,dat){
# require(dplyr)
#
#   apply(designsim, 1, function(r){
#     print(r)
#     f=fn(r[2])
#     type=fc(r[1])
#     reps=fn(r[3])
#
#     if(type=='uniform'){
#       s=cleansim(extinctCrep(pop = dat$id,frac = f, rep = reps),
#                  name = 'uniform')
#     }else if(type=='bio18dif'){
#       s=cleansim(extinctCrep(pop = dat$id,frac = f, weights = 0.01 +normalizeC(dat$ccbio18 - dat$bio18),useweights = 1, rep = reps),
#                  name = 'bio18dif')
#     }else{
#
#       wei= if(substring(type,1,1)=='-'){
#                   1-normalizeC( fn(dat[,substring(type,2,nchar(type))]) )
#            }else{normalizeC(fn(dat[,type])) }
#       wei=wei+0.01
#       s=cleansim(extinctCrep(pop = dat$id,frac = f, weights =wei, useweights = 1,rep = reps),
#            name = type)
#     }
#
#     return(s)
#
#   }) %>% do.call(.,what=rbind)
#
# }


