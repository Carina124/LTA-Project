

########################## VARIABLES THAT NEED TO BE INITIATED####################################################

# for the frequency part of the LR calculation 
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

# "_" == of 
# "." == a space between words
# "



###########################
#### START OF FUNCTION ####
###########################
LR_calculator <- function(file.name, num_contrib, num_sims, is.a.truecontrib){
  
  
  # We have to use this particular seed for reproduceability 
  set.seed (123560) 
  
  #reads in the files we are working with 
  afs_csv <- read.csv(file.name)   
  
  
  #####tabfreq makes an onject that holds the following information###
  # @TAB function reads in all the allele frequencies per loci 
  # @which.loc reports the loci that are taken into consderation
  # @pop.names are right now just reads "population" but will try to change it to the actual population (Changed)
  
  pop.afs <- tabfreq(tab = afs_csv, 
                     pop.names = as.factor(file.name)
  )
  
  
  # This is where the LRs from each iteration is stored
  truecontrib.LR.vec <- c()
  noncontrib.LR.vector <- c()
  
  
  
  ##############################
  ### true contributor code ####
  ##############################
  
  
  if (is.a.truecontrib == 0){
    j = 1 
    
    #I think this should be <=
    while (j < num_sims){
      
      #######simugeno objects store genotypes from the tabfreq ###########
      #popgen$tab.geno gives the genotyprs of all individuals (n) 
      sim.genotypes <- simugeno(tab = pop.afs, 
                                n = num_contrib, 
                                which.loc = loci.names.inorder
                                # APRIL 29, 20; NIQUO C
                                # CHANIGING THIS SO THAT THE CODE READS STRS FROM THE SAME VECTOR THROUGHOUT 
                                #c("CSF1PO",
                                  #c("CSF1PO",
                                              #"D3S1358, "D5S818" "D7S820","D8S1179 "D13S317","D16S539","D18S51","D21S11","FGA", "TH01",
                                             # "TPOX",
                                             # "VWA"
                                
      )
      
      #simulate a mixture using the simulated genotypes  
      sim.mix <- simumix(sim.genotypes,
                         ncontri = num_contrib
      ) 
      
      
      singleLR_vector <- c()
      
      k = 1 
      while (k < 14){  
        known.contrib.all.atk.pros = c()
        
        for(i in 1:num_contrib){
          known.contrib.all.atk.pros = c(known.contrib.all.atk.pros, 
                                         as.numeric(strsplit(sim.mix$mix.prof[i,k], "/")[[1]])
          )
        } 
        if   (num_contrib == 1){
          known.contrib.all.atk.def = 0
        } else {known.contrib.all.atk.def = c()
        for (i in 1:(num_contrib - 1)){
          known.contrib.all.atk.def = c(known.contrib.all.atk.def,
                                        as.numeric(strsplit(sim.mix$mix.prof[i, k], "/")[[1]])
          )
        }
        
        
        }
        known.noncontrib.all.atk.def <- c()
        known.noncontrib.all.atk.def <-  as.numeric(strsplit(sim.mix$mix.prof[num_contrib, k], "/")[[1]])
        
        ##############################
        ### True contrib Single_LR ####
        ##############################
        single_LR <<- LR( Repliste = c(sim.mix@mix.all[[k]]),
                          Tp = c(known.contrib.all.atk.pros),
                          ## Vd = Non contrib under Hd - Suspect ##
                          Td = c(known.contrib.all.atk.def),
                          Vp = 0,
                          ## Vp = Non contrib under Hp - 0 ##
                          Vd = known.noncontrib.all.atk.def, 
                          xd = 1,
                          xp = 0,
                          theta = 0,
                          prDHet = c(0.2,0.2),
                          prDHom = c(0.04,0.04),
                          prC = 0,
                          freq = pop.afs@tab[[1]][[loci.names.inorder[k]]]
        )
        
        
        
        
        
        
        #######################################
        ##### True Contrib Infinity loop #######        
        ######################################       
        while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0){
          single_LR <<- LR( Repliste = c(sim.mix@mix.all[[k]]),
                            Tp = c(known.contrib.all.atk.pros),
                            ## Vd = Non contrib under Hd - Suspect ##
                            Td = c(known.contrib.all.atk.def),
                            #tells Dorothy to return the last two digits of vector (dont think we need the strsplit commands)
                            Vp = 0,
                            ## Vp = Non contrib under Hp - 0 ##
                            Vd = known.noncontrib.all.atk.def, 
                            xd = 1,
                            xp = 0,
                            theta = 0,
                            prDHet = c(0.2,0.2),
                            prDHom = c(0.04,0.04),
                            prC = 0,
                            freq = pop.afs@tab[[1]][[loci.names.inorder[k]]]
          )
          
          
          #end of while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0)  
        }
        
        
        singleLR_vector <- c(singleLR_vector,single_LR$LR)
        
        k = k + 1
        
        #end of while (k < 14), for LR calculation   
      }
      
      # Adds the LR for this sim to the log10_LR vector 
      truecontrib.LR.vec[j] <- log10(prod(singleLR_vector))
      j = j + 1
      
      #end of while (j < num_sims)
    }
    
    #end of if (is.a.truecontrib == 0)
  }
  
  ##########################
  ### non contributors #####
  ##########################
  
  else if (is.a.truecontrib == 1){   
    #set.seed (458) 
    j = 1
    
    while (j < num_sims){ 
      
      genos.in.mix <- simugeno(tab = pop.afs, 
                               n = num_contrib, 
                               which.loc = loci.names.inorder
                               # APRIL 29, 20; NIQUO C
                               # CHANIGING THIS SO THAT THE CODE READS STRS FROM THE SAME VECTOR THROUGHOUT 
                               #c("CSF1PO",
                                # c("CSF1PO","D3S135,D5S818","D7S820","D8S1179","D13S317","D16S539","D18S51","D21S11", "FGA","TH01","TPOX", "VWA"
                               
      )
      
      sim.mix <- simumix(genos.in.mix,
                         ncontri = num_contrib
      )
      
      
      
      
      noncon.sus <- simugeno(tab = pop.afs, 
                             n =  1, 
                             which.loc = loci.names.inorder
                              # APRIL 29, 20; NIQUO C
                              # CHANIGING THIS SO THAT THE CODE READS STRS FROM THE SAME VECTOR THROUGHOUT 
                              #c("CSF1PO",
                                           #"D3S1358",
                                          # "D5S818",
                                           #"D7S820",
                                           #"D8S1179",
                                          # "D13S317",
                                          # "D16S539",
                                           #"D18S51",
                                           #"D21S11",
                                           #"FGA",
                                           #"TH01",
                                           #"TPOX",
                                           #"VWA"
                             
      )
      
      singleLR_vector <- c()
      #############################
      ### non contributors == 1+ ##
      #############################
      
      k = 1
      
      while (k < 14){  
        
        known.contrib.all.atk.pros <- c()
        
        if(num_contrib == 1){
          known.contrib.all.atk.pros = as.numeric(strsplit(noncon.sus$tab.geno[1,k], 
                                                           "/")[[1]])
        } else {
          
          for(i in 1:(num_contrib)){
            known.contrib.all.atk.pros = c( known.contrib.all.atk.pros, 
                                            as.numeric(strsplit(sim.mix$mix.prof[i,k], 
                                                                "/")[[1]])
            )
          }
          known.contrib.all.atk.pros = c(known.contrib.all.atk.pros,as.numeric(strsplit(noncon.sus$tab.geno[1,k], 
                                                                                        "/")[[1]]))
        }
        
        known.contrib.all.atk.def = c()
        
        if(num_contrib == 1){
          known.contrib.all.atk.def = 0
        } else { 
          # debugging: Niquo 2/19/20 (i in 1:(num_contribs - 1))
          for (i in 1:(num_contrib)){
            known.contrib.all.atk.def = c(known.contrib.all.atk.def,
                                          as.numeric(strsplit(sim.mix$mix.prof[i, k], "/")[[1]])
            )
          }
        }
        
        #this is the known non contributor according to the defense in this case (i.e the suspect)
        known.noncontrib.all.atk.def <- c()
        known.noncontrib.all.atk.def = c(as.numeric(strsplit(noncon.sus$tab.geno[1,k], 
                                                             "/")[[1]]))
        
        #############################
        ### non contrib Single_LR ##
        #############################
        
        single_LR <- LR( Repliste = c(sim.mix$mix.all[[k]]),
                         Tp = known.contrib.all.atk.pros,
                         Td = known.contrib.all.atk.def,
                         Vp = 0,
                         #Vd = known.noncontrib is noncon.sus
                         Vd = known.noncontrib.all.atk.def,
                         xd = 1,
                         xp = 0,
                         theta = 0,
                         prDHet = c(0.2,0.2),
                         prDHom = c(0.04,0.04),
                         prC = 0,
                         freq = pop.afs@tab[[1]][[loci.names.inorder[k]]]
        )
        
        while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0){
          print(paste("in while loop", k))
          single_LR <- LR( Repliste = c(sim.mix$mix.all[[k]]),
                           Tp = known.contrib.all.atk.pros,
                           Td = known.contrib.all.atk.def,
                           Vp = 0,
                           #Vd = known.noncontrib is noncon.sus
                           Vd = known.noncontrib.all.atk.def,
                           xd = 1,
                           xp = 0,
                           theta = 0,
                           prDHet = c(0.2,0.2),
                           prDHom = c(0.04,0.04),
                           prC = 0,
                           freq = pop.afs@tab[[1]][[loci.names.inorder[k]]]
                           
          )
        }
        
        
        singleLR_vector <- c(singleLR_vector,single_LR$LR)
        k = k + 1
        
        
        #end of  while (k < 14)  
      }
      # debugging 2/20 Niquo - I added log10
      #noncontrib.LR.vector[j] <- prod(singleLR_vector)
      
      noncontrib.LR.vector[j] <- log10(prod(singleLR_vector))
      j = j + 1
      #end of else, for more than 1 contributor 
      
    }
    
    #end of while (j < num_sims), which loops over num_sims
    #this is where the return statement goes for this loop 
    return(noncontrib.LR.vector)
  }
  
  return(truecontrib.LR.vec)
  
}



#### END OF FUNCTION #############################
###################debugging 2/18/20 - 2/19Niquo C######################################

LR_calculator <- function(file.name, num_contribs, num_sims, is.a.truecontrib)
  
  result_4contrib <- LR_calculator(pops.for.sims,4,100,1)
  result_5contrib <- LR_calculator(pops.for.sims,5,100,1)
  result_6contrib <- LR_calculator(pops.for.sims,6,100,1)
  result_7contrib <- LR_calculator(pops.for.sims,7,100,1)
  result_8contrib <- LR_calculator(pops.for.sims,8,100,1)
  result_9contrib <- LR_calculator(pops.for.sims,9,100,1)
  result_10contrib <- LR_calculator(pops.for.sims,10,100,1)
  ####################################################################################
  


########################## VARIABLES THAT NEED TO BE INITIATED####################################################

#I think this will require that we load the file "273_populations_sims.csv" to the server

pops.for.sims <- readRDS("pop.for.sims.csv")

file.exists("pop.for.sims.csv")

#used in the 3d array 
contributer = 1:10

# Niquo 2/19 changed to num_contribs throughtout  
#used as a parameter in the function and throughout
num_contrib = 10

#amount of iterations
#used as a parameter in the function and throughout

num_sims = 100


#########################################################################################################################
#################m STEP:1  Makes a 3d array for TRUE CONTRIBUTORS ################################################
#STEP:1 Is the code that makes an empty 3d arratthen passes each population to the total contrib function 




#make a variable that is num of sims 
#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
true.array_3d = array(rep(1,num_sims*num_contrib*length(pops.for.sims)), dim =c(num_sims,num_contrib,length(pops.for.sims)))



#This loop uses and x,y,z variable scheme to represent locations in 3d space 
# x is the number of contributors [,x,]
# y are the populations [,,y]
# z are the slices of 1000 num_sims [z,,]
#we need n+1 simulations
for (y in 1:length(pops.for.sims)) {
  #print(y)
  for (x in 1:length(contributer)) {
                        # file.name,num_contrib,num_sims,is.a.truecontrib)  
    result <-LR_calculator(pops.for.sims[y], x, num_sims + 1, 0)
    
    true.array_3d[,x,y] = result
  }
}

write.csv(true.array_3d,"A_truecontrib.array.csv")
saveRDS(true.array_3d, file = "A_truecontrib.array.csv")





######################################STEP:1 3d array for NON CONTRIBUTOR all populations ########################################
# STEP 1: same as above but with NON CONTRIBUTOR code



#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
non.contrib.fpr.array = array(rep(1,num_sims*num_contrib*length(pops.for.sims)), dim =c(num_sims,num_contrib,length(pops.for.sims)))




# x is the number of contributors [,x,]
# y are the populations [,,y]
# z are the slices of 1000 num_sims [z,,]
#we need n+1 simulations



for (y in 1:length(pops.for.sims)) {
  print(y)
  for (x in 1:length(contributer)) {
    result = LR_calculator(pops.for.sims[y], x, num_sims +1 , 1)
    non.contrib.fpr.array[,x,y] = result
    
  }
}
write.csv(non.contrib.fpr.array,"A_non.contrib.fpr.array.csv")
saveRDS(non.contrib.fpr.array, file = "A_non.contrib.fpr.array.csv")


#non.contrib.fpr.array

########################################## STEP:2 A Loop that counts FPR greater than -1 for the entire 3d array NON CONTRIB #######################################################################
## This loop counts the number of FPRs found in each population per contributor and stores it in the matrix FPRcount


#fprcount_mat = matrix(1:88,nrow = length(pops.for.sim),ncol = 8)

#pops.for.sims = pops.for.sims[1:57]

fprcount_mat = matrix(1:(length(pops.for.sims)*num_contrib),nrow = length(pops.for.sims),ncol = num_contrib)
for (i in 1:length(pops.for.sims)){
  # i makes the connection from population 
  #population = pops.for.sims[i]
  
  
  for(j in 1:num_contrib){
    contrib = j
    count=1
    #sticks is num_sims 
    #sticks = array_3d[,j,i]
    for(k in 1:num_sims){
      sticks = non.contrib.fpr.array[k,j,i]
      
      
      if(sticks > 0){
        
        count = count + 1
        print(count)
      }
      
      
    }
    
    
    result = (count - 1)
    # print(result)
    fprcount_mat[i,j] = result
    
  }
}


total.fpr.non.contrib.mat <- fprcount_mat
total.fpr.non.contrib.mat

#This is to name the matrix row and columns 

iterations = 1:num_contrib
colnames(total.fpr.non.contrib.mat) <- iterations
rownames(total.fpr.non.contrib.mat) <- pops.for.sims
total.fpr.non.contrib.mat

#Write the Matrix to a CSV
write.csv(total.fpr.non.contrib.mat, "A_noncontrib_mat2.csv")



########################################STEP:4 LOOP that calculates the AVERAGE FPR NON CONTRIB ############################################################
#STEP 2: This is done after LRs have been calculated for all populations and have been put into a 3darray.
#this is an empty matrix that looks through the 3d array and counts the instances where the LR is greater than -1
# the loop then divides the # of LR by the amount of iterations and it is stored in a Matrix that will be used for plotting 

fpr.average.mat = matrix(1:(length(pops.for.sims)*num_contrib),nrow = length(pops.for.sims),ncol = num_contrib)

for (i in 1:length(pops.for.sims)){
  population = pops.for.sims[i]
  
  for(j in 1:num_contrib){
    contrib = j
    count=1
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:num_sims){
      #this needs to be changed depending on which array is being used 
      sticks = non.contrib.fpr.array[k,j,i]
      
      
      if(sticks > 0){
        
        count = count + 1
        print(count)
      }
      
      
    }
    
    #this needs to be changed depending on how many iterations are used 
    result = (count - 1)/num_sims
    print(result)
    fpr.average.mat[i,j] = result
    
  }
}



non.contrib.ave.fpr_mat <- fpr.average.mat
non.contrib.ave.fpr_mat

iterations = 1:num_contrib
colnames(non.contrib.ave.fpr_mat) <- iterations
rownames(non.contrib.ave.fpr_mat) <- pops.for.sims
non.contrib.ave.fpr_mat

write.csv(non.contrib.ave.fpr_mat, "A_non.contrib.ave.fpr_mat_feb2020.csv")
