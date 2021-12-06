prop_data_na<-function(x, cr= 2, countna=TRUE){

    tmp=apply(x,cr, function(x){
      table(is.na(x))[as.character(countna)] / length(x)

    } )
    as.matrix(tmp)
}
