# install limma from bioconductor to extract desired columns from file
source("http://bioconductor.org/biocLite.R")
biocLite("limma")

# load limma
library(limma)

#code below reads desired columns into dataframe (read.columns() is function in limma package)
expHetFiledf = read.columns("Wrld_database.csv", c("Average EH", "Identifier"), sep=",")
View(expHetFiledf)

#code below creates scatterplot for expHetFiledf
library(devtools)

expHetFiledf$Identifier = as.factor(expHetFiledf$Identifier) #This line converts 'Identifier' column from numeric to a factor variable
head(expHetFiledf)
ggplot(expHetFiledf, aes(x=`Average EH`, y=Identifier, shape=Identifier, color=Identifier, size=Identifier)) + 
  geom_point() + scale_shape_manual(values=seq(0,12))
##################################################################
#code below creates barplot for expHetFiledf
library(reshape2) #package where melt function is
expHetFiledf.m = melt(expHetFiledf, id.vars = 'Identifier')
expHetFiledf.m 
ggplot(expHetFiledf.m, aes(x=Identifier, y=value, fill=variable)) + geom_bar(stat = "identity")

#Code below creates histogram
p = ggplot(expHetFiledf, aes(x=`Average EH`, group = Identifier, fill = Identifier)) + geom_histogram(color="black")
p + theme_bw() + theme_classic() + theme(legend.position = "top")




