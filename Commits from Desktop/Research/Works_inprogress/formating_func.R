formatingfunc <- function(x){
  library(readr)
  population <- read_csv(x,
                         skip = 1)
  pop_mat <- as.matrix(population)
  head(pop_mat)
  matlesscol = pop_mat[ ,-27:-53]
  finalfreqtbl <- c()
  for(i in 1:nrow(matlesscol)){
    if(matlesscol[i] == "Allele"){
      rowstop = i
      print(rowstop)
      finalfreqtbl <-matlesscol[-rowstop:-nrow(matlesscol),]
      ##return(finalfreqtbl)
      
    }
    return(finalfreqtbl)
  }}

formatingfunc("black_nist .csv")
finalfreqtbl
######################################test############################
library(readr)
population <- read_csv("Navajo-Table 1.csv", na = "0", 
                     trim_ws = FALSE, skip = 1)
View(population)
pop_mat <- as.matrix(population)
head(pop_mat)
matlesscol <- pop_mat[ , -27:-53]
matlesscol
for(i in 1:nrow(matlesscol)){
  if(matlesscol[i] == "Allele"){
    rowstop = i
    print(rowstop)
    #finalfreqtbl <-matlesscol[-rowstop:-nrow(matlesscol),]
    ##return(finalfreqtbl)
    
  }
  print(finalfreqtbl)
}
##
#######################################################################
for(i in 1:nrow(matlesscol)){
  print(1)
}

  
  if(matlesscol[i] == "Allele")
    rowstop = i
    print(rowstop)
######################################################################
library(readr)
population <- read_csv("Filipino-Table 1.csv", skip =1)  

matlesscol <-population[ ,-27:-53]
typeof(matlesscol)

if(matlesscol[i] == "Allele")
  rowstop = i
print(rowstop)


pop_mat <- as.matrix(population)
head(pop_mat)
matlesscol <- pop_mat[ , -27:-53]

for(i in 1:nrow(matlesscol)){
  print(i)
}
######################### 1/19/2019
afri.gene <- read.csv(file_name, skip =1)    
afripop <- tabfreq(tab=afri.gene, pop.names=as.factor('Afri'))
afrigeno <- simugeno(tab=afripop, n=pop_size, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                            "TPOX", "D3S1358", "D7S820", "D16S539",
                                                            "D18S51", "D21S11", "D2S1338", "D19S433",
                                                            "vWA"))
