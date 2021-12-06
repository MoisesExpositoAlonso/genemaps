
empiricalp_<-function(value,distribution){
  pemp<-as.numeric( table(distribution<value)['TRUE'] / length(distribution) )
  if(is.na(pemp)){pemp=0}
  return(pemp)
}

empiricalp<-function(values,distribution){
  sapply(values,function(value)empiricalp_(value,distribution) )
}

####************************************************************************####

freqchange<-function(cod, field=field.c, map=genomes$map,fam=genomes$fam, pdist=NULL){
	
	message("Analysing experiment ",cod)
	f<-fieldfilter(field.c,cod) %>% 
	              filter( id %in% fam$sample.ID)
	haplotypes<-fn(sapply(f$id, function(i) which( fam$sample.ID %in% i)))
	fit<-f$Fitness

	fres<-realfreqchange(Go@address,
	           fit,
	           haplotypes-1,
	           1:ncol(Go)-1)

	tmp<-dplyr::select(map, chromosome, physical.pos,allele1,allele2) %>%
			rename(chr=chromosome, pos=physical.pos) %>%
			mutate(SNP=paste0(chr, "_",pos))
	tmp$env<-cod
	tmp$fq0<-fn(fres$fq0)
	tmp$fq1<-fn(fres$fq1)
	tmp$dQ<-fn(fres$dQ)
	tmp$p_binom<-fn(fres$p_binom)

	tmp$FDR<- tmp$p_binom < fdr_level(tmp$p_binom)
  tmp$Bonf<- tmp$p_binom < 1/nrow(tmp)
 
  if(!is.null(pdist)){
    tmp$p_empir<- empiricalp(tmp$p_binom,pdist)
  }
	return(tmp)
}

####************************************************************************####
runfreqempirical<-function(cod, field=field.c, map=genomes$map,fam=genomes$fam, 
                           iter=100, sampleSNPs=10){
	message("runing permutations for empirical p-values ")
	f<-fieldfilter(field.c,cod) %>% 
	              filter( id %in% fam$sample.ID)
	haplotypes<-fn(sapply(f$id, function(i) which( fam$sample.ID %in% i)))
	fit<-f$Fitness
  
	pvals<-fn(freqempirical(Go@address,fit,haplotypes-1,sample(1:ncol(Go)-1,sampleSNPs,replace=F), iter=iter))

	return(fn(pvals))
}