library(tidyverse)
library(readr)
library(ggplot2)




################################### Raw World Data Read In #########################################
#reads in the file with 60,000 + data points 
world_data <- read_csv("world_data_clean.csv", col_names = FALSE)

view(world_data)

#turns data type to a matrix or DF first then creates a file 
write.table(world_data, file = paste("X1", "new.csv", sep ="_"), sep =",", row.names=FALSE)

# changes the generic names given to the columns when it is imported
col_names <- c("Population","Loci","Allele","Frequency")

for (i in 1:4){
  print(i)
  colnames(world_data)[i] <- col_names[i]
  print(names(world_data))
  
}

###practice subsetting
#subset_unique<- subset(, world_data[,1:3] )




#tells you the number in each group (informative)
table(world_data$Population)

#summarizes all the data that is found in world_data (less informative)
summary(world_data)

qplot(data = world_data, x = population)

################################ Loop that formats populations for Forensim ######################################

uni_pop <- (unique(world_data$Population))

#filter any pop that has less than 13 loci #keep track of the data sets we dropped #whats the minimum amount of people that 
#these allele freqs are here

dup <- world_data[which(duplicated(world_data[,1:3])),]
new_dup <- as.vector(dup[,1])
new_dup_unique <- unique(new_dup)


rmv_list <- c()

#because the i's im looking for are right next to each other
for (i in length(uni_pop):1){
 # print(i)
  if(uni_pop[i] %in% new_dup_unique[[1,]]){
    print(paste(i,"is a match ",uni_pop[i]))
    uni_pop <- uni_pop[-i]
    
  }
  else{
    #print(paste(i,"is false"))
  }
}
  




#initiates a vector that will catch all the file names passed through the loop 
file_names <- rep(NA,length(uni_pop))
file_names <-c()

for (pop in 1:length(uni_pop)) {
  #print(paste("pop is",pop))
  pop_lines <- which(world_data$Population == uni_pop[pop])
  #print(paste("pop lines",pop_lines))
  test <-spread(world_data[pop_lines,], key=Loci,value=Frequency)  %>% 
    mutate_at(vars(-Population), funs( ifelse(is.na(.),0,.)))
  #print(paste("test is",test))
  file_name <- paste(uni_pop[pop], "new.csv", sep ="_")
  #print(paste("file name is",file_name))
  #define file_names as an empty vector
  file_names <- c(file_names,file_name)
  #write.table(test, file = paste(uni_pop[pop], "new.csv", sep ="_"), sep =",", row.names=FALSE)
  
  
}
write.csv(file_names,"Original_427_Filenames.csv")
save(list = character(), file = "Original_427_Filenames.RData")
######DATA: For the other team D2441 #######
for (pop in 1:length(uni_pop)) {
  #print(paste("pop is",pop))
  pop_lines <- which(world_data$Population == uni_pop[pop])
  #print(paste("pop lines",pop_lines))
  test <-spread(world_data[pop_lines,], key=Loci,value=Frequency)  %>% 
    mutate_at(vars(-Population), funs( ifelse(is.na(.),0,.)))
  #print(paste("test is",test))
  file_name <- paste(uni_pop[pop], "new.csv", sep ="_")
  #print(paste("file name is",file_name))
  #define file_names as an empty vector
  file_names <- c(file_names,file_name)
  write.table(test, file = paste(uni_pop[pop], "new.csv", sep ="_"), sep =",", row.names=FALSE)
  
  
}



    
      
  


