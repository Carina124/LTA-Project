twoapa = total_LR_plot_vector3
twoapa_non = total_LR_plot_vector_non3

plot.new()
par(mfrow=c(2,2))
hist(twoapa, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (2P Mix Apache)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(twoapa_non, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################3 person mix############################################
threeapa = total_LR_plot_vector3
threeapa_n = total_LR_plot_vector_non3

hist(threeapa, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (3P Mix Apache)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(threeapa_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################4 person mix############################################
fourapa = total_LR_plot_vector3
fourapa_n = total_LR_plot_vector_non3

hist(fourapa, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (4P Mix Apache)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(fourapa_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################5 person mix############################################
fiveapa = total_LR_plot_vector3
fiveapa_n = total_LR_plot_vector_non3

hist(fiveapa, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (5P Mix Apache)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(fiveapa_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################ 6 person mix############################################
sixapa = total_LR_plot_vector3
sixapa_n = total_LR_plot_vector_non3

hist(sixapa, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (6P Mix Apache)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(sixapa_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################ 7 person mix############################################
sevapa = total_LR_plot_vector3
sevapa_n = total_LR_plot_vector_non3

hist(sevapa, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (7P Mix Apache)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(sevapa_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################ 8 person mix############################################
eiapa = total_LR_plot_vector3
eiapa_n = total_LR_plot_vector_non3

hist(eiapa, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (8P Mix Apache)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(eiapa_n, breaks=seq(from=-1,to=300,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
