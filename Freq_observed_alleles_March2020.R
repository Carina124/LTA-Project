######## Checking Assigned Alleles vs. Actual allele Frequency ###########

num_contrib = 100

file.name <-  "USAf.1_new.csv"


loci.names.inorder <- c("CSF1PO",
                        "D3S1358",
                        "D5S818",
                        "D7S820",
                        "D8S1179",
                        "D13S317",
                        "D16S539",
                        "D18S51",
                        "D21S11",
                        "FGA",
                        "TH01",
                        "TPOX",
                        "VWA"
)


###############################Total contrib function MATT's CODE ##############################################################
library(forensim)
library(readr)
set.seed (123560) 
  #reads in the files we are working with 
  afs_csv <- read.csv(file.name)   
  
  pop.afs <- tabfreq(tab = afs_csv, 
                     pop.names = as.factor(file.name)
  )
  
  genos.in.mix <- simugeno(tab = pop.afs, 
                           n = num_contrib, 
                           which.loc = c("CSF1PO",
                                         "D3S1358",
                                         "D5S818",
                                         "D7S820",
                                         "D8S1179",
                                         "D13S317",
                                         "D16S539",
                                         "D18S51",
                                         "D21S11",
                                         "FGA",
                                         "TH01",
                                         "TPOX",
                                         "VWA"
                           )
  )

####comparing genos.in.mix to imported allele frequncey table 
# the population name is the only thing not hard coded ###
genos.in.mix$tab.freq$USAf.1_new.csv[k]
pop.afs$tab$USAf.1_new.csv[k]
#### Make a plan####
# Loop 1: 


k =13
known.contrib.all.atk.pros <- c()
for(i in 1:(num_contrib)){
  known.contrib.all.atk.pros = c( known.contrib.all.atk.pros, 
                                  as.numeric(strsplit(sim.mix$mix.prof[i,k], 
                                                      "/")[[1]]))}

#### loop that zeros out all the positions in the tab freq objest's list ##############
#double square bracket allows you to go into a list. 
#you can only go deeper into a list if the thing before it has double square brackets

#INDEXING genos.in.mix$tab.freq[[allele frequency all loci]][[freq of specific STR]][value at each allele]

# genos.in.mix needs to be saved to a diff name because we still need the original data from the forensim list
empty_list = genos.in.mix$tab.freq$USAf.1_new.csv

for (locus in 1:13){
  # loop that goes through each STR in the list 
 
  empty_list[[locus]]
  print(locus)
  
  for (position in 1:length(empty_list[[locus]])){
    
    #This loop sets each position at a particular STR to o 0
    
    empty_list[[locus]][[position]] = 0
  }
}

print(empty_list)

###### Next steps: loop that writes ratios of alleles assigned to empty_list
observed.allele.count[[1]][1] = 0
observed.allele.count[[1]]["9"]


freq_assigned_vec <- c()
j = 1
for (k in 1:13){
for (j in 1:length(genos.in.mix@mix.all[[k]])){
  #looks in
  print(paste("allele",sim.mix@mix.all$VWA[j]))
  current_allele = sim.mix@mix.all$VWA[j]
  count = 0
  #observed.allele.count[[k]][current_allele] = observed.allele.count[[k]][current_allele] + 1
  for(i in 1:length(known.contrib.all.atk.pros)){
    if (known.contrib.all.atk.pros[i] == sim.mix@mix.all$v){
      count = count + 1
      observed.allele.count[[k]][current_allele] = observed.allele.count[[k]][current_allele] + 1
      
      
      
    }
    
  }
  #print(count)
  freq_assigned = count/(length(known.contrib.all.atk.pros))
  #print(freq_assigned)
  freq_assigned_vec <- c(freq_assigned_vec,freq_assigned)
  #ass.freq.matrix <- ass.freq.matrix[j,1]
}
  
}