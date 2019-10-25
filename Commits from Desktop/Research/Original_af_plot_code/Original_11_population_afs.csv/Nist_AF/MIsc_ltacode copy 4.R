#histogram code from MAtt
h3 = total_slr_vector.4
h4 = total_slr_vector.5
hist(h3, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,0.9),
     main="False Positive (1 Person Mixture)", xlab="Likelihood Ratio (Log10)")
hist(h4, breaks=seq(from=-1,to=25,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
###############################################################################################
########noncntrib########
afri.gene <- read.csv('afri.gene.csv')    ## Open file holding allele data ##
afripop <- tabfreq(tab=afri.gene, pop.names=as.factor('Afri'))
afrigeno <- simugeno(tab=afripop, n=100, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                       "TP0X", "D3S1358", "DS7820", "D16S539",
                                                       "D18S51", "D21S11", "D2S1338", "D19S433",
                                                       "VWA"))
noncontrib_Afri <- simugeno(tab=afripop, n=1, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                            "TP0X", "D3S1358", "DS7820", "D16S539",
                                                            "D18S51", "D21S11", "D2S1338", "D19S433",
                                                            "VWA"))
noncon.sus <- simumix(noncontrib_Afri, ncontri = 1)
fga.2 <- afripop$tab$Afri$FGA           ## Inititialize variables for vector ##
csf.2 <- afripop$tab$Afri$CSF1PO
th01.2 <- afripop$tab$Afri$TH01
tp0x.2 <- afripop$tab$Afri$TP0X
vwa.2 <- afripop$tab$Afri$VWA
d3s.2 <- afripop$tab$Afri$D3S1358
d7s.2 <- afripop$tab$Afri$DS7820
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


j = 1
set.seed(657)
while (j < 11){
  mix.2 <- simumix(afrigeno,ncontri = 2)
  
  slr_vector <- c()  
  total_slr_vector <- c()
  
  k = 1
  while (k < 14){
    slr <- LR(Repliste=c(mix.2$mix.all[[k]]),                         
              Tp=c(as.numeric(strsplit(mix.2$mix.prof[1,k], "/")[[1]]),
                   as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])), 
              Td=c(as.numeric(strsplit(mix.2$mix.prof[1,k], "/")[[1]])),   
              Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
              xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
              prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
    
    
    slr_vector <- c(slr_vector,slr$LR)
    k = k + 1
  }
  print(j)
  
  slr_total.3[j] <- (abs(slr_vector[1])*abs(slr_vector[2])*abs(slr_vector[3])*
                       abs(slr_vector[4])*abs(slr_vector[5])*abs(slr_vector[6])*
                       abs(slr_vector[7])*abs(slr_vector[8])*abs(slr_vector[9])*
                       abs(slr_vector[10])*abs(slr_vector[11])*abs(slr_vector[12])*
                       abs(slr_vector[13]))
  
  total_slr_vector.3 <- c(log10(slr_total.3))
  
  total_slr_vector.3[total_slr_vector.3 < 0.0] <- -1
  
  
  print(total_slr_vector.3)
  j = j + 1
}



