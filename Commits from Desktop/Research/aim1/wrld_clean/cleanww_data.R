library(tidyverse)
library(readr)
library(ggplot2)




################################### Raw World Data Read In #########################################
#reads in the file with 60,000 + data points 
world_data <- read_csv("world_data_clean.csv", col_names = FALSE)

view(world_data)

#turns data type to a matrix or DF first then creates a file 
write.table(exphetero_mat_clean, file = paste("X1", "new.csv", sep ="_"), sep =",", row.names=FALSE)

# changes the generic names given to the columns when it is imported
col_names <- c("Population","Loci","Allele","Frequency")

for (i in 1:4){
  print(i)
  colnames(world_data)[i] <- col_names[i]
  print(names(world_data))
  
}



#subset_unique<- subset(, world_data[,1:3] )

#tells you the number in each group (informative)
table(world_data$Population)

#summarizes all the data that is found in world_data (less informative)
summary(world_data)

qplot(data = world_data, x = population)

################################ Loop that formats populations for Forensim ######################################

#creates a list of all the unique populations in the dataset 
uni_pop <- unique(world_data$Population)


#saves data points that are duplicated to a vector, dup
dup <- world_data[which(duplicated(world_data[,1:3])),]

#saves an unlisted version of unique populations from the dup into a vector, unlisted_duplicates 
unlisted_dup_pops <- unlist(unique(dup[,1]))

# The loop removes the duplicated populations in uni_pop by looking through unlisted_dup_pops for mathcing pops
for (i in length(uni_pop):1){
 if (uni_pop[i] %in% unlisted_dup_pops){
   uni_pop = uni_pop[-i]
   print(paste(i,"is a match"))
 }
  
}

#I have to re-remove these matches 



#initiates a vector that will catch all the file names passed through the loop 
file_names<- c()


for (pop in 1:length(uni_pop)) {
  #saves the data for a particular population into pop_lines
  pop_lines <- which(world_data$Population == uni_pop[pop])
  
  #spreads the data from pop_lines to the format that we want
  pop_data_spread <-spread(world_data[pop_lines,], key=Loci,value=Frequency)  %>% 
    mutate_at(vars(-Population), funs( ifelse(is.na(.),0,.)))
  
  # file naming format: population_new.csv
  file_name <- paste(uni_pop[pop], "new.csv", sep ="_")
  #saves each file name to this vector that grows at each iteration; for later use
  file_names <- c(file_names,file_name)
  #creates files 
  if( sum(colnames(pop_data_spread)==  "vWA")){
    location <-which(colnames(pop_data_spread)== "vWA")
    colnames(pop_data_spread)[location]<-"VWA"
  }
 # write.table(pop_data_spread[,2:dim(pop_data_spread)[2]], file = file_name, sep =",", row.names=FALSE)
  
  
}

#trying to save file_names as an object in r
save(list = file_names, file = "/Users/niquoceberio/Desktop/Research/aim1/wrld_clean")
write.table(file_names, file = file_name_list, sep =",", row.names=FALSE)

###Matrix that sorts populations####
pop_matrix <- matrix(data = NA, nrow = 13, ncol = 13)
write.table(exphetero_mat_clean, file = exp_hetero_clean, sep = ",")
