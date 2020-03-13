######## Checking Assigned Alleles vs. Actual allele Frequency ###########

#### Niquo 03/12/2020 Turning this code into a function######

Freq_observed <- function(file.name, num_contrib){



######## These need to be intiatied first 

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



  sim.mix <- simumix(genos.in.mix,
                     ncontri = num_contrib
  )
  

#### loop that zeros out all the positions in the tab freq objest's list ##############
#double square bracket allows you to go into a list. 
#you can only go deeper into a list if the thing before it has double square brackets

#INDEXING genos.in.mix$tab.freq[[allele frequency all loci]][[freq of specific STR]][value at each allele]

# genos.in.mix needs to be saved to a diff name because we still need the original data from the forensim list



empty_list = genos.in.mix$tab.freq[[file.name]]

for (locus in 1:13){
  # loop that goes through each STR in the list 
 
  empty_list[[locus]]
  
  
  for (position in 1:length(empty_list[[locus]])){
    
    #This loop sets each position at a particular STR to o 0
    
    empty_list[[locus]][[position]] = 0
  }
}



for (k in 1:13){
  # this looks at each STR in the list 1-13 
  
  #generate the allel
  known.contrib.all.atk.pros <- c()
  for(i in 1:(num_contrib)){
    known.contrib.all.atk.pros = c( known.contrib.all.atk.pros, 
                                    as.numeric(strsplit(sim.mix$mix.prof[i,k], 
                                                        "/")[[1]]))}
  #goes into the list at each STR 1=13
  current_loci = unique(known.contrib.all.atk.pros)

for (j in 1:length(current_loci)){
  ## this is how many alleles there are in a particular population 
  #print(paste("allele",
  current_allele =  current_loci[j]
  count = 0


  #observed.allele.count[[k]][current_allele] = observed.allele.count[[k]][current_allele] + 1
  for(l in 1:length(known.contrib.all.atk.pros)){
    if (known.contrib.all.atk.pros[l] == current_allele){
     count = (count + 1) 
      #empty_list[[k]][current_allele] =  empty_list[[k]][current_allele] + 1
      #this line goes into a specific allele and frequency @ STR K 
      #genos.in.mix$tab.freq[[file.name]][[k]][[j]]))
    # empty_list[[k]][[as.character(current_allele)]] = count
      
    }
    
  }
  count = (count+1)/length(known.contrib.all.atk.pros)
  empty_list[[k]][[as.character(current_allele)]] <- count
}
  
}
  Freq_actually_obs_list <<- empty_list
}

# Freq_observed (file.name, num_contrib)
#num_contrib = 100
#file.name <-  "USAf.1_new.csv"







