library(forensim)
data("strveneto")
?data("strusa")
data("Tu")
tabfreq("Hisp")
################hisp 
hist(h1, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (2P Mix Hispanic)", xlab="Likelihood Ratio (Log10)", ylab = "Frequency")
hist(h2, breaks=seq(from=-1,to=25,by=1), col=rgb(0,0,1,0.5), ylim=c(0,1.0), prob = TRUE, add=T)
box()
?hist

hist(h1, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (2P Mix Hispanic)", xlab="Likelihood Ratio (Log10)", ylab = "Frequency")
hist(h2, breaks=seq(from=-1,to=25,by=1), col=rgb(0,0,1,0.5), ylim=c(0,1.0), prob = TRUE, add=T)
box()