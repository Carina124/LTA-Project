#################################### Working Matrix (Friday w/ Maria) 2/22/19 ################################################################
library(readr)

#filter any pop that has less than 13 loci #keep track of the data sets we dropped #whats the minimum amount of people that 
#these allele freqs are here

#vector that holds the 11 populations in it 
freq_vec <- c("Afri_paper.csv","Apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
              "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")



#loci_vec <- c("FGA","CSF1PO","TH01","D8S1179","D3S1358","D7S820","D16S539",
              #"D21S11","D2S1338","D19S433",
              #"vWA","D18S51","TPOX")
loci_vec <- c("CSF1PO","D3S1358","D5S818","D7S820","D8S1179","D13S317","D16S539","D18S51","D21S11",
            "FGA","TH01","TPOX","VWA") #"D1S165","D2S441","D2S1338","D10S1248","D12S391","D19S433","D22S1045")
#Vectors for world data
file_names




#intiates the vector that holds all  heterozygosity of each loci  
exp_hvec <- c()
exp_hmat.loci <- matrix(data = NA, nrow = length(file_names), ncol = length(loci_vec))

#assigns the col name to the matrix from the freq vector (changes)
colnames(exp_hmat.loci) <- loci_vec

#asigns the row name to the matrix (changes)
rownames(exp_hmat.loci) <- file_names
#
pops.out.of.bounds <- c()
#loop that iterates through the 13 loci in Loci_vec
for (locus_i in 1:length(loci_vec)){
  current_locus <- loci_vec[locus_i]
  
  
  #loop iterates through each CSV file   
  for( pop_i in 1:length(file_names)){
   #fixed the col spec outputs
    afs_matrix <- as.matrix(read.csv(file_names[pop_i]))
   
    
    #loop iterates through colums of file and checks if it matches current loci 
    for (i in 2:dim(afs_matrix)[2]){
      name =  colnames(afs_matrix)[i]
      
      if(name == current_locus){
        cur_afs <- afs_matrix[,i]
        
        if (sum(as.numeric(cur_afs)) >= 1.01 || sum(as.numeric(cur_afs))<=.99){
          print(paste("Sum is out of bounds",loci_vec[locus_i],file_names[pop_i],sum(as.numeric(cur_afs))))
          pops.out.of.bounds <- c(pops.out.of.bounds,file_names[pop_i])
        }
       
        afs_sqrd <- as.numeric(cur_afs)^2
        #loci_sum <- sum(afs_sqrd)
        exp_hetero <- (1 -(sum(afs_sqrd)))
        #where you are putting it/what you are putting in it 
        exp_hmat.loci[pop_i,locus_i] <- exp_hetero
        #print(paste("heterozygosity for locus  ",current_locus,"for pop ",j, "is ",exp_hetero))
      }
    }
  }
}

#these populations dont add up to one so they will be removed 
pops.tobe.removed <-unique(pops.out.of.bounds)
for (i in length(file_names):1){
  if (file_names[i] %in% pops.tobe.removed){
    file_names = file_names[-i]
    print(paste(i,"is a match"))
  }
  
}



write.csv(exp_hmat.loci,"Exp_hetero_wrld.csv")


#Calculate average expected heterozygosity for each population 

mean_vec <- c()

for( i in 1:length(exp_hmat.loci)){
mean_exp <- mean(exp_hmat.loci[i,])
mean_vec <-c(mean_vec,mean_exp)
}


# turns the vector into a matrix 
average_exphetero_mat <- as.matrix(mean_vec, ncol=2)

colnames(average_exphetero_mat)[1] <- "Average Expected Heterozygosity" 



#asigns the row name to the matrix (changes)
rownames(average_exphetero_mat) <- file_names



###### Removing populations without the STR core ############
colSums(is.na(average_exphetero_mat))
sum(is.na(average_exphetero_mat))
average.eh.matrix <- na.omit(average_exphetero_mat)
write.csv(average.eh.matrix,"Average_Expected_Hetero_9.19.csv")

#THESE ARE the populations that we will use for this experiment 
pops.for.sims <- rownames(average.eh.matrix)
write.csv(pops.for.sims,"Names of 273 populations for sims")

save(list = character(),file = "Names of 273 populations for sims.RData")
save(pops.for.sims)





