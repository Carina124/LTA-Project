####################################################################AFRICAN BASE###################################################

library(readr)
black_nist <- read_csv("black_nist .csv", 
                        skip = 1)
black.matrix <- as.matrix(black_nist)
black.gene= black.matrix[ ,-31:-91]
black.gene
#View(black_nist_)

## Open file holding allele data ##
afripop <- tabfreq(tab=black.gene, pop.names=as.factor('Afri'))
afrigeno <- simugeno(tab=afripop, n=100, which.loc = c("CSF1PO","TH01" ,"FGA", "D8S1179",
                                                       "TPOX", "D3S1358", "D7S820", "D16S539",
                                                       "D18S51", "D21S11", "D2S1338", "D19S433",
                                                       "vWA"))
noncontrib_Afri <- simugeno(tab=afripop, n=1, which.loc = c("CSF1PO","TH01" ,"FGA", "D8S1179",
                                                            "TPOX", "D3S1358", "D7S820", "D16S539",
                                                            "D18S51", "D21S11", "D2S1338", "D19S433",
                                                            "vWA"))
noncon.sus <- simumix(noncontrib_Afri, ncontri = 1)

fga.2 <- afripop$tab$Afri$FGA           ## Inititialize variables for vector ##
csf.2 <- afripop$tab$Afri$CSF1PO
th01.2 <- afripop$tab$Afri$TH01
tpox.2 <- afripop$tab$Afri$TPOX
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
freq_vector[[5,1]] <- c(tpox.2)
freq_vector[[6,1]] <- c(d3s.2)
freq_vector[[7,1]] <- c(d7s.2)
freq_vector[[8,1]] <- c(d16.2)
freq_vector[[9,1]] <- c(d18.2)
freq_vector[[10,1]] <- c(d21.2)
freq_vector[[11,1]] <- c(d2s.2)
freq_vector[[12,1]] <- c(d19.2)
freq_vector[[13,1]] <- c(vwa.2)

########################################Asian AF #################################################################
library(readr)
asian_nist <- read_csv("asian_nist.csv", 
                       skip = 1)
asian.matrix <- as.matrix(asian_nist)
asian.gene= asian.matrix[ ,-31:-91]


   ## Open file holding allele data ##
asianpop <- tabfreq(tab=asian.gene, pop.names=as.factor('asian'))
asiangeno <- simugeno(tab=asianpop, n=100, which.loc = c("CSF1PO","TH01" ,"FGA", "D8S1179",
                                                       "TPOX", "D3S1358", "D7S820", "D16S539",
                                                       "D18S51", "D21S11", "D2S1338", "D19S433",
                                                       "vWA"))
#noncontrib_asian <- simugeno(tab=asianpop, n=1, which.loc = c("CSF1PO","TH01" ,"FGA", "D8S1179",
                                                             "TPOX", "D3S1358", "D7S820", "D16S539",
                                                             "D18S51", "D21S11", "D2S1338", "D19S433",
                                                             "vWA"))
#noncon.sus <- simumix(noncontrib_asian, ncontri = 1)

########################################################################################
fga.2 <- asianpop$tab$asian$FGA           ## Inititialize variables for vector ##
csf.2 <- asianpop$tab$asian$CSF1PO
th01.2 <- asianpop$tab$asian$TH01
tpox.2 <- asianpop$tab$asian$TPOX
vwa.2 <- asianpop$tab$asian$vWA
d3s.2 <- asianpop$tab$asian$D3S1358
d7s.2 <- asianpop$tab$asian$D7S820
d8s.2 <- asianpop$tab$asian$D8S1179
d16.2 <- asianpop$tab$asian$D16S539
d18.2 <- asianpop$tab$asian$D18S51
d21.2 <- asianpop$tab$asian$D21S11
d2s.2 <- asianpop$tab$asian$D2S1338
d19.2 <- asianpop$tab$asian$D19S433
freq2_vector <- matrix(list(), nrow=13, ncol=1)  ## Creating Vector ##
freq2_vector[[1,1]] <- c(d8s.2)
freq2_vector[[2,1]] <- c(th01.2)
freq2_vector[[3,1]] <- c(fga.2)
freq2_vector[[4,1]] <- c(csf.2)
freq2_vector[[5,1]] <- c(tpox.2)
freq2_vector[[6,1]] <- c(d3s.2)
freq2_vector[[7,1]] <- c(d7s.2)
freq2_vector[[8,1]] <- c(d16.2)
freq2_vector[[9,1]] <- c(d18.2)
freq2_vector[[10,1]] <- c(d21.2)
freq2_vector[[11,1]] <- c(d2s.2)
freq2_vector[[12,1]] <- c(d19.2)
freq2_vector[[13,1]] <- c(vwa.2)





j = 1
set.seed(657)
while (j < 101){
  mix.2 <- simumix(afrigeno,ncontri = 2)
  
  slr_vector <- c()  
  total_slr_vector <- c()
  total_slr_vector.3 <- c()
  slr_total.3 <- c()
  z = 1
  while (z < 13){
    slr <- LR(Repliste=c(mix.2$mix.all[[z]]),                         
              Tp=c(as.numeric(strsplit(mix.2$mix.prof[1,z], "/")[[1]]),
                   as.numeric(strsplit(noncon.sus$mix.prof[1,z], "/")[[1]])), 
              Td=c(as.numeric(strsplit(mix.2$mix.prof[1,z], "/")[[1]])),   
              Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,z], "/")[[1]])),
              xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
              prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[z,1]])
    
    
    slr_vector <- c(slr_vector,slr$LR)
    z = z + 1
  }
  print(j)
 #[j] here 
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
total_slr_vector.3
length(slr_total.3)
#################################black/black histogram#################################################
h1 =total_slr_vector.3
h1

hist(h1, breaks=seq(from=-1,to=5,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),xlim=c(-1,1.0),
     main="A. American AF, A. American NonContri", xlab="Likelihood Ratio (Log10)", ylab = "frequency")
##########################################Asian/Asian Hist#############################################
h1 =total_slr_vector.3
h1

hist(h1, breaks=seq(from=-1,to=5,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),xlim=c(-1,1.0),
     main="Asian AF, A. Asian NonContri", xlab="Likelihood Ratio (Log10)", ylab = "frequency")
######################################Black/Asian######################################################
h1 =total_slr_vector.3
h1

hist(h1, breaks=seq(from=-1,to=5,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),xlim=c(-1,1.0),
     main="A. American AF, A. American NonContri", xlab="Likelihood Ratio (Log10)", ylab = "frequency")
#####################################Asian/Black#######################################################
h1 =total_slr_vector.3
h1

hist(h1, breaks=seq(from=-1,to=5,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),xlim=c(-1,1.0),
     main="A. American AF, A. American NonContri", xlab="Likelihood Ratio (Log10)", ylab = "frequency")