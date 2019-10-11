twotrini = total_LR_plot_vector3
twotrini_n = total_LR_plot_vector_non3

plot.new()
par(mfrow=c(2,2))
hist(twotrini, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (2P Mix Trinidad)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(twotrini_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################3 person mix############################################
threeafri = total_LR_plot_vector3
threeafri_n = total_LR_plot_vector_non3

hist(threeafri, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (3P Mix Trinidad)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(threeafri_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################4 person mix############################################
fourafri = total_LR_plot_vector3
fourafri_n = total_LR_plot_vector_non3

hist(fourafri, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (4P Mix Trinidad)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(fourafri_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################5 person mix############################################
fiveafri = total_LR_plot_vector3
fiveafri_n = total_LR_plot_vector_non3

hist(fiveafri, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (5P Mix Trinidad)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(fiveafri_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################ 6 person mix############################################
sixafri = total_LR_plot_vector3
sixafri_n = total_LR_plot_vector_non3

hist(sixafri, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (6P Mix Trinidad)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(sixafri_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################ 7 person mix############################################
sevafri = total_LR_plot_vector3
sevafri_n = total_LR_plot_vector_non3

hist(sevafri, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (7P Mix Trinidad)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(sevafri_n, breaks=seq(from=-1,to=20,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
################################ 8 person mix############################################
eiapa = total_LR_plot_vector3
eiapa_n = total_LR_plot_vector_non3

hist(eiapa, breaks=seq(from=-1,to=25,by=1), col=rgb(1,0,0,0.5), prob = TRUE, ylim=c(0,1.0),
     main="False Postive (8P Mix Trinidad)", xlab="Likelihood Ratio (Log10)",ylab = "Frequency")
hist(eiapa_n, breaks=seq(from=-1,to=300,by=1), col=rgb(0,0,1,0.5), ylim=c(0,0.9), prob = TRUE, add=T)
box()
