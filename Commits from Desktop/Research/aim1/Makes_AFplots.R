library(readr)

#vector of loci we will iterate through
loci_vec <- c("FGA","TPOX","CSF1PO","TH01","D8S1179","D3S1358","D7S820","D16S539",
              "D21S11","D2S1338","D19S433",
              "vWA","D18S51")


#contains all the files we want to plot 
freq_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
              "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")

#starts the graphics device driver for PDF graphics 
pdf("plots_11popnozero", width=20, height=20)

#This is how many plots per row,column 
par(mfrow=c(3,3))

#position loop- loci 
for (k in 1:length(loci_vec)){
  current_loci <- loci_vec[k]
  
  #nested for look will turn all the csv's into a matrix
  ##will remove NA's and plot the alleles vs. Frequency & write them to a PDF 
  for( j in 1:length(freq_vec)){
    file_name <- read_csv(freq_vec[j], na = "0")
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
