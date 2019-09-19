
singleLRs_vec_product <- c()
LR_vector <- c()
while (k < 14){
  
  single_LR <- LR(Repliste=c(sim.mix$mix.all[[1]]), #k         
                                                              #K
                  Tp=c(as.numeric(strsplit(sim.mix$mix.prof[1,1], "/")[[1]])), 
                                                                        #K
                  Td=0,Vp=0,Vd=c(as.numeric(strsplit(sim.mix$mix.prof[1,1], "/")[[1]])),
                  
                  xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                                                              #k
                  prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[1,1]])
  
  LR_vector <- c(LR_vector,single_LR$LR)
  k = k + 1
}
  singleLRs_vec_product[j] <<- log10(prod(LR_vector))
  
  j = j + 1
  
  return(log10_totalLRs_vec)
  
  
  k=1
  list_Alleles = c()
  while (k < 14){  
    ## Calculated 13 times for 13 loci ##
    #c = num_contribs - 1  
   
    for(i in 1:num_contribs){
      ### QUESTION: why do we need to use strsplit ? Rename these alleles so they make sense 
      list_Alleles = c(list_Alleles, as.numeric(strsplit(sim.mix$mix.prof[i,k], "/")[[1]]))
      print(paste("loci",k,list_Alleles))
      k = k + 1
    }
  }
  