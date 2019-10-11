h1 = total_LR_plot_vector1
h2 =total_LR_plot_vector_non1

plot.new()
par(mfrow=c(2,2))
hist(h1, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (2P Mix Hispanic)", xlab="Likelihood Ratio (Log10)", ylab = "Frequency")
hist(h2, breaks=seq(from=-1,to=25,by=1), col=rgb(0,0,1,0.5), ylim=c(0,1.0), prob = TRUE, add=T)
box()

h3 = total_LR_plot_vector1
h4 = total_LR_plot_vector_non1


hist(h3, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (3P Mix Hispanic)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(h4, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
######################################### 4 person mix #####################################################
h5 = total_LR_plot_vector3
h6 = total_LR_plot_vector_non1


hist(h5, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (4P Mix Hispanic)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(h6, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
######################################### 5 person mix #####################################################
h7 = total_LR_plot_vector3
h8 = total_LR_plot_vector_non1


hist(h7, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (5P Mix Hispanic)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(h8, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
######################################### 6 person mix #####################################################
sixperson = total_LR_plot_vector3
sixperson_non = total_LR_plot_vector_non3


hist(sixperson, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (6P Mix Hisp)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(sixperson_non, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
######################################### 7 person mix #####################################################
sevenperson = total_LR_plot_vector3
sevenperson_non = total_LR_plot_vector_non3


hist(sevenperson, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (7P Mix Hisp)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(sevenperson_non, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()

######################################### 8 person mix #####################################################
eiperson = total_LR_plot_vector3
eiperson_non = total_LR_plot_vector_non3


hist(eiperson, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (8P Mix Hisp)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(eiperson_non, breaks=seq(from=-1,to=300,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
######################################### 1 person  #####################################################
onep = total_LR_plot_vector3
onep = total_LR_plot_vector_non3


hist(onep, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (7P Mix His)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(h6, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()

