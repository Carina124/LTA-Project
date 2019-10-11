#########This for loop plots Allele vs. Frequency from multiple CSV files

library(readr)

#contains all the files we want to plot 
freq_vec <- c("afri_NIST.csv","Apache.csv","Bahamian.csv","Caucasian.csv","Chamorro.csv","Filipino.csv",
              "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trini.csv")

#starts the graphics device driver for PDF graphics 
pdf("11_pop_plot2", width=20, height=20)

#This is how many plots per row,column 
par(mfrow=c(5,5))


#nested for look will turn all the csv's into a matrix
##will remove NA's and plot the alleles vs. Frequency & write them to a PDF 
for( j in 1:length(freq_vec)){
  file_name <- read_csv(freq_vec[j])
  #file_name <- read_csv("afri.gene.csv")
  file_mat <- as.matrix(file_name)
  #dev.off()
  
  #start of the nested for loop
  for (i in 2:dim(file_mat)[2]){
    freq <- file_mat[,i]
    isthere = !is.na(freq)
    #print(isthere)
    #print(freq[isthere])
    alleles <- file_mat[,1]
    #isthereallele <-!is.na(alleles)
    #print(alleles[isthere])
    #length(isthereallele)
    name =  colnames(file_mat)[i]
    #print(i)
    plot(x= alleles[isthere], y = freq[isthere], type = "h",lwd = 4,col="red", xlab = "Alleles",
         ylab = "frequency",main = paste(name, freq_vec[j]), xlim = c(0,55), ylim=c(0,1))
  
    
  }
}

#stops writing to the PDF 
dev.off()

alleles[isthere]
#####################################New loop to plot loci together################################################
library(readr)

#vector of loci we will iterate through
loci_vec <- c("FGA","TPOX","CSF1PO","TH01","D8S1179","D3S1358","D7S820","D16S539",
             "D21S11","D2S1338","D19S433",
              "vWA","D18S51")


#contains all the files we want to plot 
freq_vec <- c("Afri_paper.csv","Apache.csv","Bahamian.csv","Caucasian.csv","Chamorro.csv","Filipino.csv",
              "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trini.csv")

#starts the graphics device driver for PDF graphics 
pdf("11_loci_D18S51", width=20, height=20)

#This is how many plots per row,column 
par(mfrow=c(5,5))

#position loop- loci 
for (k in 1:length(loci_vec)){
  current_loci <- loci_vec[k]

#nested for look will turn all the csv's into a matrix
##will remove NA's and plot the alleles vs. Frequency & write them to a PDF 
for( j in 1:length(freq_vec)){
  file_name <- read_csv(freq_vec[j])
  #file_name <- read_csv("afri.gene.csv")
  file_mat <- as.matrix(file_name)
  #dev.off()
  
  #start of the nested for loop
  for (i in 2:dim(file_mat)[2]){
    name =  colnames(file_mat)[i]
    
    if(name == current_loci){
      freq <- file_mat[,i]
      isthere = !is.na(freq)
      alleles <- file_mat[,1]
      name =  colnames(file_mat)[i]
      plot(x= alleles[isthere], y = freq[isthere], type = "h",lwd = 4,col="red", xlab = "Alleles",
           ylab = "frequency",main = paste(name, freq_vec[j]), xlim = c(0,55), ylim=c(0,1))
    }
  }
}
}
#stops writing to the PDF 
dev.off()


############################Expected Heterozygosity Code############################################################
### This loop goes calculates expected heterozygosity ############


sigma_eq <- c()
for (i in 2:dim(file_mat)[2]){
  loci_vector <- file_mat[,i]
  print(loci_vector)
  isthere = !is.na(loci_vector)

  isthere_loci <- loci_vector[isthere]
  print(isthere_loci)
  loci_sqrd <- ((as.numeric(isthere_loci))^2)
  print(loci_sqrd)
  loci_sum <- sum(loci_sqrd)
  sigma_eq <- c(loci_sum,sigma_eq)
  exp_homo <- (1/length(sigma_eq)*sigma_eq)
  exp_heteroZ <- (1-exp_homo)
  
}







#################################### Matrix building (w/ Carina) 4/21/19 ################################################################

#vector of the 13 loci that we want to calculate exp heterozygosity for 
loci_vec <- c("FGA","TPOX","CSF1PO","TH01","D8S1179","D3S1358","D7S820","D16S539",
              "D21S11","D2S1338","D19S433",
              "vWA","D18S51")

#vector that holds the 11 populations in it 
freq_vec <- c("Afri_paper.csv","Apache.csv","Bahamian.csv","Caucasian.csv","Chamorro.csv","Filipino.csv",
              "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trini.csv")

# empty matrix that we will use to 
hetero_mat <- matrix(nrow = 13, ncol = length(freq_vec))

#assgns the col name to the matrix from the freq vector (changes)
colnames(hetero_mat) <- freq_vec


#asigns the row name to the matrix (changes)
rownames(hetero_mat) <- loci_vec



#initializes sigma_eq
sigma_eq <- c()

#loop that iterates through the 13 loci in Loci_vec
for (k in 1:length(loci_vec)){
  current_loci <- loci_vec[k]

  #loop iterates through each CSV file   
  for( j in 1:length(freq_vec)){
    file_name <- read_csv(freq_vec[j])
    file_mat <- as.matrix(file_name)
 
    #loop iterates through colums of file and checks if it matches current loci 
for (i in 2:dim(file_mat)[2]){
  name =  colnames(file_mat)[i]
  
  if(name == current_loci){
  loci_vector <- file_mat[,i]
  print(loci_vector)
  isthere = !is.na(loci_vector)
  isthere_loci <- loci_vector[isthere]
  loci_sqrd <- ((as.numeric(isthere_loci))^2)
  loci_sum <- sum(loci_sqrd)
  ####opposite orer of te file --switch those two ---ask Rori 
  sigma_eq <- c(sigma_eq,loci_sum)
  exp_homo <- (1/length(sigma_eq[i])*sigma_eq)
  exp_heteroZ <- (1-exp_homo)
  #print(exp_heteroZ)
  
  }
}
}
}
hetero_mat2 <- matrix(exp_heteroZ,ncol=11,nrow = 1, byrow = T)
hetero_mat2

hetero_mat2 <- matrix(exp_heteroZ,ncol=11,nrow = 13, byrow = T)
hetero_mat2

colnames(hetero_mat2) <- freq_vec


#asigns the row name to the matrix (changes)
rownames(hetero_mat2) <- loci_vec


for(l in 1:length(feq_vec)){
  mat[l,] <- exp_heteroZ
}


#turns exp hetero into a matrix 
exp_mat <-as.matrix(exp_heteroZ)
#################################### Working Matrix (Friday w/ Maria) 2/22/19 ################################################################


#vector that holds the 11 populations in it 
freq_vec <- c("Afri_paper.csv","Apache.csv","Bahamian.csv","Caucasian.csv","Chamorro.csv","Filipino.csv",
              "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trini.csv")



loci_vec <- c("FGA","TPOX","CSF1PO","TH01","D8S1179","D3S1358","D7S820","D16S539",
                  "D21S11","D2S1338","D19S433",
                  "vWA","D18S51")

#intiates the vector that holds all exp heterozygosity values 
exp_hvec <- c()

#loop that iterates through the 13 loci in Loci_vec
for (k in 1:length(loci_vec)){
  current_loci <- loci_vec[k]
  
  #loop iterates through each CSV file   
  for( j in 1:length(freq_vec)){
    file_name <- read_csv(freq_vec[j])
    file_mat <- as.matrix(file_name)
    
    #loop iterates through colums of file and checks if it matches current loci 
    for (i in 2:dim(file_mat)[2]){
      name =  colnames(file_mat)[i]
      
      if(name == current_loci){
        loci_vector <- file_mat[,i]
        isthere = !is.na(loci_vector)
        isthere_loci <- loci_vector[isthere]
        print(isthere_loci)
        #work on the math
        loci_sqrd <- as.numeric(isthere_loci)^2
        #print(loci_sqrd)
        loci_sum <- sum(loci_sqrd)
        #print(loci_sum)
        #print(sigma_eq)
        exp_hetero <- (1 -(loci_sum))
        exp_hvec <- c(exp_hvec,exp_hetero)
        print(paste("heterozygosity for locus  ",loci_vec[k],"for pop ",freq_vec[j], "is ",exp_hetero))
      }
    }
  }
}


# empty matrix that we will use to 
hetero_mat <- matrix(exp_hvec,nrow = 11, ncol = length(loci_vec))

#assgns the col name to the matrix from the freq vector (changes)
colnames(hetero_mat) <- loci_vec


#asigns the row name to the matrix (changes)
rownames(hetero_mat) <- freq_vec


#Calculate average expected heterozygosity for each population 
hetero_mat
mean(hetero_mat[1,])