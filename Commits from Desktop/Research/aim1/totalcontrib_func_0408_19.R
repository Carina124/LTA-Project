total_contrib_function <- function(file_name,pop_size,num_contribs,iterations,non_or_true_contrib) 
{
  set.seed(657)
  afri.gene <- read.csv(file_name)    
  afripop <- tabfreq(tab=afri.gene, pop.names=as.factor('Afri'))
  afrigeno <- simugeno(tab=afripop, n=pop_size, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                              "TPOX", "D3S1358", "D7S820", "D16S539",
                                                              "D18S51", "D21S11", "D2S1338", "D19S433",
                                                              "vWA"))
  
  
  noncontrib_Afri <- simugeno(tab=afripop, n=1, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                              "TPOX", "D3S1358", "D7S820", "D16S539",
                                                              "D18S51", "D21S11", "D2S1338", "D19S433",
                                                              "vWA"))
  #print(afrigeno@tab.geno)
  noncon.sus <- simumix(noncontrib_Afri, ncontri = 1)
  
  
  
  fga.2 <- afripop$tab$Afri$FGA           ## Inititialize variables for vector ##
  csf.2 <- afripop$tab$Afri$CSF1PO
  th01.2 <- afripop$tab$Afri$TH01
  tp0x.2 <- afripop$tab$Afri$TPOX
  vwa.2 <- afripop$tab$Afri$vWA
  d3s.2 <- afripop$tab$Afri$D3S1358
  d7s.2 <- afripop$tab$Afri$D7S820
  d8s.2 <- afripop$tab$Afri$D8S1179
  d16.2 <- afripop$tab$Afri$D16S539
  d18.2 <- afripop$tab$Afri$D18S51
  d21.2 <- afripop$tab$Afri$D21S11
  d2s.2 <- afripop$tab$Afri$D2S1338
  d19.2 <- afripop$tab$Afri$D19S433
  
  freq_vector <- matrix(list(), nrow=13, ncol=1)  ## Creating Vector ##
  freq_vector[[1,1]] <- c(d8s.2)
  freq_vector[[2,1]] <- c(th01.2)
  freq_vector[[3,1]] <- c(fga.2)
  freq_vector[[4,1]] <- c(csf.2)
  freq_vector[[5,1]] <- c(tp0x.2)
  freq_vector[[6,1]] <- c(d3s.2)
  freq_vector[[7,1]] <- c(d7s.2)
  freq_vector[[8,1]] <- c(d16.2)
  freq_vector[[9,1]] <- c(d18.2)
  freq_vector[[10,1]] <- c(d21.2)
  freq_vector[[11,1]] <- c(d2s.2)
  freq_vector[[12,1]] <- c(d19.2)
  freq_vector[[13,1]] <- c(vwa.2)
  
  LR_vector_true <- c()
  total_LR <- c()
  plot_vector <- c()
  
  if (non_or_true_contrib == 0){
    
    j = 1         
    while (j < iterations){
      mix.func <- simumix(afrigeno,ncontri = num_contribs)
      
      LR_vector <- c()
      total_LR_vector <- c()
      
      if (num_contribs == 1){
        
        k = 1
        while (k < 14){
          single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                          Tp=c(as.numeric(strsplit(mix.func$mix.prof[1,k], "/")[[1]])),  
                          Td=0,Vp=0,Vd=c(as.numeric(strsplit(mix.func$mix.prof[1,k], "/")[[1]])),
                          xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                          prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
          #print(single_LR$LR)
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
        
      } else {  
        k = 1                             ## LR Calculation. ##
        while (k < 14){  
          ## Calculated 13 times for 13 loci ##
          c = num_contribs - 1  
          listOfAlleles = c()
          for(i in 1:num_contribs){
            listOfAlleles = c(listOfAlleles, as.numeric(strsplit(mix.func$mix.prof[i,k], "/")[[1]]))
            listofAllels.Defense = c()
            for (i in 1:c){
              listofAllels.Defense = c(listofAllels.Defense, as.numeric(strsplit(mix.func$mix.prof[i,k], "/")[[1]]))   
              
            }
          }
          single_LR <- LR(Repliste=c(mix.func@mix.all[[k]]),
                          Tp=c(listOfAlleles),
                          Td=c(listofAllels.Defense), ## Vd = Non contrib under Hd - Suspect ##
                          Vp=0,Vd=c(as.numeric(tail(listOfAlleles,2), "/")[[1]]), ## Vp = Non contrib under Hp - 0 ##
                          xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                          prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
          while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0) {
            single_LR <- LR(Repliste=c(mix.func@mix.all[[k]]),
                            Tp=c(listOfAlleles),
                            Td=c(listofAllels.Defense), ## Vd = Non contrib under Hd - Suspect ##
                            Vp=0,Vd=c(as.numeric(tail(listOfAlleles,2), "/")[[1]]), ## Vp = Non contrib under Hp - 0 ##
                            xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                            prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
            
           # print(afrigeno@tab.geno)
            #LR_vector <- c(LR_vector,single_LR$LR)
            #if(single_LR$LR>0)
            #if(is.nan(single_LR$LR)==FALSE){
            #  print(paste("Lr for iteration ", j, "allele", k, "is not nan it is ", single_LR$LR))
            #break
            #LR_vector <- c(LR_vector,single_LR$LR)
            
          } 
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      }
      LR_vector_true[j] <- prod(LR_vector)
      
      j = j + 1
      total_LR_vector <- c(log10(LR_vector_true))
      total_LR_plot_vectornav6 <<- total_LR_vector 
      g = 1
      while (g < iterations){
        afri.gene <- read.csv(file_name)    
        afripop <- tabfreq(tab=afri.gene, pop.names=as.factor('Afri'))
        afrigeno <- simugeno(tab=afripop, n=pop_size, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                                    "TPOX", "D3S1358", "D7S820", "D16S539",
                                                                    "D18S51", "D21S11", "D2S1338", "D19S433",
                                                                    "vWA"))
        g = g + 1
      }
    }
  } else if (non_or_true_contrib == 1) {
    j = 1
    #set.seed(657)
    while (j < iterations){ #loop over simulations
      mix.func <- simumix(afrigeno,ncontri = num_contribs)
      
      LR_vector <- c()  
      total_LR_vector <- c()
      
      if (num_contribs == 1){
        k = 1
        while (k < 14){
          single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                          Tp=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),  
                          Td=0,Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                          xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                          prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
          #print(listOfAlleles)
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      } else {
        k = 1
        while (k < 14){  
          c = num_contribs - 1
          listOfAlleles = c()
          for(i in 1:num_contribs){
            listOfAlleles = c(listOfAlleles, as.numeric(strsplit(mix.func$mix.prof[i,k], "/")[[1]]))
            listofAllels.Defense = c()
            for (i in 1:c){
              listofAllels.Defense = c(listofAllels.Defense, as.numeric(strsplit(mix.func$mix.prof[i,k], "/")[[1]]))
              
              
            }
          }  
          repeat{
            single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                            Tp=c(listOfAlleles,  
                                 as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])), 
                            Td=c(listofAllels.Defense),   
                            Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                            xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                            prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
            
            #LR_vector <- c(LR_vector,single_LR$LR)
            if(is.nan(single_LR$LR)==FALSE){
              break
            }
          }  
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      }     
      total_LR[j] <- (LR_vector[1]*LR_vector[2]*LR_vector[3]*
                        LR_vector[4]*LR_vector[5]*LR_vector[6]*
                        LR_vector[7]*LR_vector[8]*LR_vector[9]*
                        LR_vector[10]*LR_vector[11]*LR_vector[12]*
                        LR_vector[13])
      
      total_LR_vector <- c(log10(total_LR))
      
      total_LR_vector[total_LR_vector < 0.0] <- -1
      
      total_LR_plot_vector_nonnav8 <<- total_LR_vector
      j = j + 1
    }
  }
  return(total_LR_vector)
}
LRs = total_contrib_function('Trinidadian.csv', 50, 4, 11, 0)
LRs


g = 1
while (g < 10){
  afri.gene <- read.csv('Trinidadian.csv')    
  afripop <- tabfreq(tab=afri.gene, pop.names=as.factor('Afri'))
  afrigeno <- simugeno(tab=afripop, n=50, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                        "TPOX", "D3S1358", "D7S820", "D16S539",
                                                        "D18S51", "D21S11", "D2S1338", "D19S433",
                                                        "vWA"))
  g = g + 1
  #print(afrigeno@tab.geno)
}

