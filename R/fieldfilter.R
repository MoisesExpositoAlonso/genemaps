
fieldfilter<-function(d,code='mli',...){

     if(code=='mhi') return(  dplyr::filter(d, site=='madrid', water=='h', indpop=='i',...) )
else if(code=='mhp') return(  dplyr::filter(d, site=='madrid', water=='h', indpop=='p',...) )
else if(code=='mli') return(  dplyr::filter(d, site=='madrid', water=='l', indpop=='i',...) )
else if(code=='mlp') return(  dplyr::filter(d, site=='madrid', water=='l', indpop=='p',...) )

else if(code=='thi') return(  dplyr::filter(d, site=='tuebingen', water=='h', indpop=='i',...) )
else if(code=='thp') return(  dplyr::filter(d, site=='tuebingen', water=='h', indpop=='p',...) )
else if(code=='tli') return(  dplyr::filter(d, site=='tuebingen', water=='l', indpop=='i',...) )
else if(code=='tlp') return(  dplyr::filter(d, site=='tuebingen', water=='l', indpop=='p',...) )
else stop('code not recognized!')


}
