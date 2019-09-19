library(readr)
Chamorrofinal <- read_csv("Chamorrofinal.csv", 
                          na = "empty", skip = 1)
View(Chamorrofinal)

#turns file into a matrix
pop_mat <- as.matrix(Chamorro)

head(Chamorro)

#removes unwanted columns
matlesscol = pop_mat[ ,-18:-43]

#removes unwanted rows
finalfreqtbl_chamorro<-matlesscol[-100:-172,]








#histogram code for 

hist(h3, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,0.9),
     main="False Positive (2 Person Mixture)", xlab="Likelihood Ratio (Log10)")










