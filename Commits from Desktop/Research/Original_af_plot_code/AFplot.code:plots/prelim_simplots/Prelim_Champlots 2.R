onecham = total_LR_plot_vector3
onecham_non = total_LR_plot_vector_non3

plot.new()
par(mfrow=c(2,2))
hist(onecham, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (7P Mix His)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(onecham_non, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
#####################################2 person mix##############################################
twocham = total_LR_plot_vector3
twocham_non = total_LR_plot_vector_non3


hist(twocham, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (2P Mix Chamorro)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(twocham_non, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
#####################################3 person mix##############################################
threecham = total_LR_plot_vector3
threecham_non = total_LR_plot_vector_non3


hist(threecham, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (3P Mix Chamorro)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(threecham_non, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
#####################################4 person mix##############################################
fourcham = total_LR_plot_vector3
fourcham_n = total_LR_plot_vector_non3


hist(fourcham, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (4P Mix Chamorro)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(fourcham_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
#####################################5 person mix##############################################
fivecham = total_LR_plot_vector3
fivecham_n = total_LR_plot_vector_non3


hist(fivecham, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (4P Mix Chamorro)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(fivecham_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
#####################################6 person mix##############################################
sixcham = total_LR_plot_vector3
sixcham_n = total_LR_plot_vector_non3


hist(sixcham, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (6P Mix Chamorro)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(sixcham_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
#####################################y person mix##############################################
sevcham = total_LR_plot_vector3
sevcham_n = total_LR_plot_vector_non3


hist(sevcham, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (7P Mix Chamorro)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(sevcham_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
##################################### 8 person mix##############################################
eicham = total_LR_plot_vector3
eicham_n = total_LR_plot_vector_non3


hist(eicham, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (8P Mix Chamorro)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(eicham_n, breaks=seq(from=-1,to=300,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
