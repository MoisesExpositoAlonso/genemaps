## NEED TO INCLUDE THIS SOMEWHERE
cnapplot<-function(x, y, xpval, ypval,
                   apthreshold=cnthreshold*2,
                   cnthreshold=0.01,xlab='',ylab='' ,
  # ap='#8856a7' , cn='#2166ac' ,rest='grey',
  ap='black' , cn='darkgrey' ,rest='white',
  alpha=1, totsnps=1353386){

    type<- x
    type[] <- 1
    type[(x>0 & y > 0) | (x < 0 & y < 0)] <- 1 # The same values


    cn1<-(xpval<cnthreshold & ypval>cnthreshold)
    cn2<-(ypval<cnthreshold & xpval>cnthreshold)
    type[ cn1| cn2 ] <- 2 # conditional neutrality

    type[
      (xpval<apthreshold & ypval<apthreshold) &
        (
          (x <0 & y>0) | (x>0 & y<0)
        )
      ] <- 3 #antagonistic pleiotropy

    tabtmp<-matrix(c(
                  0 , fn(table(cn2)['TRUE']),
                  fn(table(cn1)['TRUE']), fn(table(type==3)['TRUE'] )
                  ), ncol=2,nrow=2)
    tabtmp[is.na(tabtmp)] <-0
    tab<-tabtmp
    tab[1,1] = totsnps - sum(tabtmp)
    tab<-as.table(tab)

    tab

    message("Total SNPs cnap:")
    print(table(type>1))

    fisher<- fisher.test(tab)
    fisher
    ods<-fisher$estimate %>% format(digits=3)
    pval<-fisher$p.value %>% format(digits=3,scientific=T)
    res<-paste('odds ratio = ',ods, ', p=', pval)
    message(res)

    resultsplot<-ggplot(data.frame(tab)) +
      geom_text(aes(x = Var1, y=Var2,label=Freq), size=5) +
      labs(x='Significant for environment 1', y='Significant for environment 2')+
          theme(
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank()) +
      ggtitle(res)



    tp<-data.frame(x=x,y=y,z=type)
    tmp<-ggplot() +
      xlab(xlab)+ ylab(ylab)+
      # geom_point(data=filter(tp,z <2), aes(y=y,x=x) ,
      #   color=cn,
      #   shape=1)+
      geom_point(data=filter(tp,z==2), aes(y=y,x=x) ,
        color=transparent(cn,alpha),
        shape=19)+
      geom_point(data=filter(tp,z==3), aes(y=y,x=x) ,
        color=transparent(ap,alpha),
        shape=19)+

      geom_hline(yintercept = 0, lty='dashed', col='darkgrey')+geom_vline(xintercept = 0,lty='dashed', col='darkgrey')

return( list(plot=plot_grid(tmp,resultsplot ),
             fisher=c(OR=ods,p=pval),
              ap=which(type==3),
              cn=which(type==2),
              threshold=cnthreshold
               )
     )

}

## NEED TO INCLUDE THIS SOMEWHERE
dirindir<-function(x, y, xpval, ypval, xlab='',ylab='' ,
  # ap='#8856a7' , cn='#2166ac' ,rest='grey',
  ap='black' , cn='darkgrey' ,rest='white',
  alpha=1, totsnps=1353386){

    type<- x
    type[] <- 1
    # type[(x>0 & y > 0) | (x < 0 & y < 0)] <- 1 # The same values


    cn1<-(xpval<0.01 & ypval==0)
    cn2<-(ypval>0 & xpval>0.01)
    type[ cn1| cn2 ] <- 2 # conditional neutrality

    type[
      (xpval<0.01 & ypval>0) ] <- 3 #antagonistic pleiotropy

    tabtmp<-matrix(c(
                  0 , fn(table(cn2)['TRUE']),
                  fn(table(cn1)['TRUE']), fn(table(type==3)['TRUE'] )
                  ), ncol=2,nrow=2)
    tabtmp[is.na(tabtmp)] <-0
    tab<-tabtmp
    tab[1,1] = totsnps - sum(tabtmp)
    tab<-as.table(tab)

    tab

    fisher<- fisher.test(tab)
    fisher
    ods<-fisher$estimate %>% format(digits=3)
    pval<-fisher$p.value %>% format(digits=3,scientific=T)
    res<-paste('odds ratio = ',ods, ', p=', pval)
    message(res)

    resultsplot<-ggplot(data.frame(tab)) +
      geom_text(aes(x = Var1, y=Var2,label=Freq), size=5) +
      labs(x='Significant differential (marginal)', y='Non-zero inclusion probability (conditional)')+
          theme(
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank()) +
      ggtitle(res)



    tp<-data.frame(x=x,y=y,z=type)
    tmp<-ggplot() +
      xlab(xlab)+ ylab(ylab)+
      # geom_point(data=filter(tp,z <2), aes(y=y,x=x) ,
      #   color=cn,
      #   shape=1)+
      geom_point(data=filter(tp,z==2), aes(y=y,x=x) ,
        color=transparent(cn,alpha),
        shape=19)+
      geom_point(data=filter(tp,z==3), aes(y=y,x=x) ,
        color=transparent(ap,alpha),
        shape=19)+

      geom_hline(yintercept = 0, lty='dashed', col='darkgrey')+geom_vline(xintercept = 0,lty='dashed', col='darkgrey')

return(list(
            plot=tmp,
             fisher=c(OR=ods,p=pval)
            ))

}

#

cnaptest<-function(n=100){
  x=rnorm(n,0,1)
  y=rnorm(n,0,2)
  xpval=rnorm(n,0,1)
  xpval[xpval<0] <-1e-6
  ypval=rnorm(n,0,1)
  ypval[ypval<0] <-1e-6

  print(cnapplot(x,y,xpval,ypval))
}
