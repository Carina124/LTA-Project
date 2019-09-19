#################################### Working Matrix (Friday w/ Maria) 2/22/19 ################################################################
library(readr)

#vector that holds the 11 populations in it 
freq_vec <- c("Afri_paper.csv","Apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
              "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")



#loci_vec <- c("FGA","CSF1PO","TH01","D8S1179","D3S1358","D7S820","D16S539",
              #"D21S11","D2S1338","D19S433",
              #"vWA","D18S51","TPOX")
loci_vec <- c("CSF1PO","D3S1358","D5S818","D7S820","D8S1179","D13S317","D16S539","D18S51","D21S11",
            "FGA","TH01","TPOX","VWA","D1S165","D2S441","D2S1338","D10S1248","D12S391","D19S433","D22S1045")
#Vectors for world data
file_names

for (i in file_names) {
  print(i)
}
#intiates the vector that holds all exp heterozygosity values 
exp_hvec <- c()
exp_hmat <- matrix(data = NA, nrow = length(file_names), ncol = length(loci_vec))
#

#loop that iterates through the 13 loci in Loci_vec
for (locus_i in 1:length(loci_vec)){
  current_locus <- loci_vec[locus_i]
  #print(current_locus)
  
  
  #loop iterates through each CSV file   
  for( pop_i in 1:length(file_names)){
    #print(file_names[pop_i])
    #see if you can stop printing col specifications
    afs_matrix <- as.matrix(read.csv(file_names[pop_i]))
    #replace 
    #file_mat <- as.matrix(file_names)
    
    #loop iterates through colums of file and checks if it matches current loci 
    for (i in 3:dim(afs_matrix)[2]){
      #print(paste("in 3rd loop"))
      name =  colnames(afs_matrix)[i]
      
      if(name == current_locus){
        cur_afs <- afs_matrix[,i]
        if (sum(cur_afs) >= 1.01 || sum(cur_afs)<=.99){
          print(paste("Sum is too high",locus_vec[locus_i],file_names[pop_i]))
        }
        #isthere = !is.na(loci_vector)
        #isthere_loci <- loci_vector[isthere]
        afs_sqrd <- as.numeric(cur_afs)^2
        #loci_sum <- sum(afs_sqrd)
        exp_hetero <- (1 -(sum(afs_sqrd)))
        #where you are putting it/what you are putting in it 
        exp_hmat[pop_i,locus_i] <- exp_hetero
        #print(paste("heterozygosity for locus  ",current_locus,"for pop ",j, "is ",exp_hetero))
      }
    }
  }
}



# empty matrix that we will use to 
hetero_mat <- matrix(exp_hvec,nrow = file_names, ncol = length(loci_vec))

#assgns the col name to the matrix from the freq vector (changes)
colnames(hetero_mat) <- loci_vec


#asigns the row name to the matrix (changes)
rownames(hetero_mat) <- pop_vec2

write.csv(hetero_mat,"ExpectedHeteroMay.csv")


#Calculate average expected heterozygosity for each population 
mean_vec <- c()

for( i in 1:length(hetero_mat)){
mean_exp <- mean(hetero_mat[i,])
print(mean_exp)
mean_vec <-c(mean_vec,mean_exp)
}


# turns the vector into a matrix 
exphetero_mat <- as.matrix(mean_vec)

colnames(exphetero_mat) <- "Average Expected Heterozygosity"


#asigns the row name to the matrix (changes)
rownames(exphetero_mat) <- pop_vec2

write.csv(exphetero_mat,"Average_Expected_HeteroMay.csv")

