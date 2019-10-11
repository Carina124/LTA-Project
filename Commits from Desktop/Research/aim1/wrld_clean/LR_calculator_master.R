

########################## VARIABLES THAT NEED TO BE INITIATED####################################################
 
#I think this will require that we load the file "273_populations_sims.csv" to the server

pops.for.sims <- readRDS("273_populations_sims.csv")



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


#num of contribs
contributer = 1:8

###############################Total contrib function MATT's CODE ##############################################################
library(forensim)
library(readr)

# "_" == of 
# "." == a space between words
# "



###########################
#### START OF FUNCTION ####
###########################
LR_calculator <- function(file.name, num_contribs, num_sims, is.a.truecontrib){
  
  # set.seed(657) 
  # We have to use this particular seed for reproduceability 
  #set.seed (123560) 
  
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
      #set.seed (123560) 
      #######simugeno objects store genotypes from the tabfreq ###########
      #popgen$tab.geno gives the genotyprs of all individuals (n) 
      sim.genotypes <- simugeno(tab = pop.afs, 
                                n = num_contribs, 
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
      
     #simulate a mixture using the simulated genotypes  
     sim.mix <- simumix(sim.genotypes,
                        ncontri = num_contribs
                        ) 
     
     
     singleLR_vector <- c()
     
       k = 1 
       while (k < 14){  
         known.contrib.all.atk.pros = c()
         
         for(i in 1:num_contribs){
           known.contrib.all.atk.pros = c(known.contrib.all.atk.pros, 
                                 as.numeric(strsplit(sim.mix$mix.prof[i,k], "/")[[1]])
                                 )
         } 
         if   (num_contribs == 1){
           known.contrib.all.atk.def = 0
         } else {known.contrib.all.atk.def = c()
         for (i in 1:(num_contribs - 1)){
           known.contrib.all.atk.def = c(known.contrib.all.atk.def,
                               as.numeric(strsplit(sim.mix$mix.prof[i, k], "/")[[1]])
                               )
         }
        
           
         }
         known.noncontrib.all.atk.def <- c()
         known.noncontrib.all.atk.def <-  as.numeric(strsplit(sim.mix$mix.prof[num_contribs, k], "/")[[1]])
        
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
    set.seed (458) 
    j = 1
  
    while (j < num_sims){ 
      
      genos.in.mix <- simugeno(tab = pop.afs, 
                                n = num_contribs, 
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
                         ncontri = num_contribs
      )
      
    
      
      
      noncon.sus <- simugeno(tab = pop.afs, 
                                n =  1, 
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
      
      singleLR_vector <- c()
      #############################
      ### non contributors == 1+ ##
      #############################
      
        k = 1
        set.seed(657)
        while (k < 14){  
         
         known.contrib.all.atk.pros <- c()
         
         if(num_contribs == 1){
           known.contrib.all.atk.pros = as.numeric(strsplit(noncon.sus$tab.geno[1,k], 
                                                            "/")[[1]])
         } else {
          
            for(i in 1:(num_contribs - 1)){
           known.contrib.all.atk.pros = c( known.contrib.all.atk.pros, 
                                   as.numeric(strsplit(sim.mix$mix.prof[i,k], 
                                                       "/")[[1]])
                                   )
            }
           known.contrib.all.atk.pros = c(known.contrib.all.atk.pros,as.numeric(strsplit(noncon.sus$tab.geno[1,k], 
                                                                                         "/")[[1]]))
         }
        
            known.contrib.all.atk.def = c()
            
            if(num_contribs == 1){
              known.contrib.all.atk.def = 0
            } else { 
            for (i in 1:(num_contribs - 1)){
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

      noncontrib.LR.vector[j] <- prod(singleLR_vector)
      j = j + 1
      #end of else, for more than 1 contributor 
      
      }
      
    #end of while (j < num_sims), which loops over num_sims
    #this is where the return statement goes for this loop 
    return(noncontrib.LR.vector)
    }
    
    return(truecontrib.LR.vec)
   
  }


#########################
#### END OF FUNCTION ####

#TEST AREA 
#file.name,num_contribs,num_sims,is.a.truecontrib)  
LRs <- LR_calculator("Africa_new.csv", 2, 100, 0)


#######################################################
##Clean -- Niquo 10/10/19##
#######################################################


#########################################################################################################################
#################m STEP:1 updated 9.19.19 solo Makes a 3d array for TRUE CONTRIBUTORS ################################################
#STEP:1 Is the code that makes an empty 3d arratthen passes each population to the total contrib function 



#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
true.array_3d = array(rep(1,10*8*length(pops.for.sims)), dim =c(10,8,length(pops.for.sims)))

#view empty array
true.array_3d 


#This loop uses and x,y,z variable scheme to represent locations in 3d space 
# x is the number of contributors [,x,]
# y are the populations [,,y]
# z are the slices of 1000 num_sims [z,,]
#we need n+1 simulations
for (y in 1:length(pops.for.sims)) {
  
  for (x in 1:length(contributer)) {
    result <-LR_calculator(pops.for.sims[y], x, 11, 0)
    print(y)
    true.array_3d[,x,y] = result
  }
}


######################################STEP:1 3d array for NON CONTRIBUTOR all populations 3/16/19########################################
# STEP 1: same as above but with NON CONTRIBUTOR code



#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
non.contrib.fpr.array = array(rep(1,10*8*length(pops.for.sims)), dim =c(10,8,length(pops.for.sims)))

#view empty array
non.contrib.fpr.array 

# x is the number of contributors [,x,]
# y are the populations [,,y]
# z are the slices of 1000 num_sims [z,,]
#we need n+1 simulations


ptm <- proc.time()
for (y in 1:length(pops.for.sims)) {
  
  for (x in 1:length(contributer)) {
    print(paste("population", pops.for.sims[y]))
    LR_calculator(pops.for.sims[y], x, 11, 1)
    non.contrib.fpr.array[,x,y] = result
    
  }
}

proc.time() - ptm
non.contrib.fpr.array

########################################## STEP:2 Loop that counts FPR greater than -1 for the entire 3d array #######################################################################
## This loop counts the number of FPRs found in each population per contributor and stores it in the matrix FPRcount

fprcount_mat = matrix(1:88,nrow = length(pops.for.sims),ncol = 8)

for (i in 1:11){
  population = pops.for.sims[i]
  
  for(j in 1:8){
    contrib = j
    count=1
    #sticks is num_sims 
    #sticks = array_3d[,j,i]
    for(k in 1:1000){
      sticks = non.contrib.fpr.array[k,j,i]
      
      
      if(sticks > -1){
        
        count = count + 1
        print(count)
      }
      
      
    }
    
    
    result = (count - 1)
    print(result)
    fprcount_mat[i,j] = result
    
  }
}

total.fpr.true.contrib.mat <-fprcount_mat
total.fpr.non.contrib.mat <- fprcount_mat


#This is all naming schemes for the matrix row and columns 
iterations = 1:8
colnames(total.fpr.true.contrib.mat) <- iterations
rownames(total.fpr.true.contrib.mat) <- pops.for.sims
total.fpr.true.contrib.mat

#Write the Matrix to a CSV
write.csv(total.fpr.true.contrib.mat, "truecontrib_mat2.csv")

iterations = 1:8
colnames(total.fpr.non.contrib.mat) <- iterations
rownames(total.fpr.non.contrib.mat) <- pops.for.sims
total.fpr.non.contrib.mat

write.csv(non_fpr_mat, "noncontrib_mat2.csv")

################################# STEP:3 PLOT for  FPR  counts of TRUE & NON CONTRIBUTORS#####################################################################
#STEP 3 takes the non contrib and true contrib matrixes and puts them on one plot 


#This is the vector used to make the different colored lines for all of the populations 
colorvec = c( "blue","red","green","purple","violet","yellow","pink","black")

plot.new()

#This intiates the plot 
plot(x=c(),y=c(),xlim = c(1,8),ylim = c(0,1000),main =" False Positives  (1000 iterations)",ylab = "false postive", xlab = 
       "number of contributors")


#this loop adds points onto the intiated plot with i reffering to the location in the matrix 
for(i in 1:8){
  points(non_fprcount_mat[i,], col= colorvec[i], lty=1)
  points(true_fprcount_mat[i,], col=colorvec[i],lty=1)
}

legend("center", legend=c("AA","Apache","Bahamian","White","Chamorro","Filipino",
                          "Jamaican","Navajo","SEH","SWH","Trinidadian"),
       col=c( "blue","red","green","purple","violet","yellow","pink","black"), lty=1:2, cex=0.8)



########################################STEP:4 LOOP that calculates the AVERAGE FPR carina 3/20/19############################################################
#STEP 2: This is done after LRs have been calculated for all populations and have been put into a 3darray.
#this is an empty matrix that looks through the 3d array and counts the instances where the LR is greater than -1
# the loop then divides the # of LR by the amount of iterations and it is stored in a Matrix that will be used for plotting 

fpr.average.mat = matrix(1:88,nrow = length(pops.for.sims),ncol = 8)

for (i in 1:11){
  population = pops.for.sims[i]
  
  for(j in 1:8){
    contrib = j
    count=1
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:1000){
      #this needs to be changed depending on which array is being used 
      sticks = non.contrib.fpr.array[k,j,i]
      
      
      if(sticks > -1){
        
        count = count + 1
        print(count)
      }
      
      
    }
    
    #this needs to be changed depending on how many iterations are used 
    result = (count - 1)/10000
    print(result)
    fpr.average.mat[i,j] = result
    
  }
}

true.contrib.ave.fpr_mat <-fpraverage_mat
non.contrib.ave.fpr_mat <- fpraverage_mat
################################# STEP:5  PLOT for AVerage FPR  counts of TRUE & NON CONTRIBUTORS#####################################################################
#STEP 3 takes the non contrib and true contrib matrixes and puts them on one plot 


#This is the vector used to make the different colored lines for all of the populations 
colorvec = c( "blue","red","green","purple","violet","yellow","pink","black")



#This intiates the plot 
plot.new()
plot(x=c(),y=c(),xlim = c(1,8),ylim = c(0,1000),main =" False Positives  (1000 iterations)",ylab = "false postive", xlab = 
       "number of contributors")


#this loop adds points onto the intiated plot with i reffering to the location in the matrix 
for(i in 1:8){
  points(non.contrib.ave.fpr_mat[i,], col= colorvec[i], lty=1)
  points(true.contrib.ave.fpr_mat[i,], col=colorvec[i],lty=1)
}

#creates a lengend the colors need to be checked to make sure there are no duplicate colors!
legend("center", legend=c("AA","Apache","Bahamian","White","Chamorro","Filipino",
                          "Jamaican","Navajo","SEH","SWH","Trinidadian"),
       col=c( "blue","red","green","purple","violet","yellow","pink","black"), lty=1:2, cex=0.8)


################### STEP 6: PLtos FPR greater than -1 vs Expected heterozygozity of the 7th contrib###################
# STEP 6 creates an empty matrix then adds the average expected heterozygosity for each population to column 1 
# And adds oncontributor count of values greater than -1 for the 7th contributor in all populations to column 2 
# Next we plot the matrix

FPR7_matrix <-matrix(1:11,nrow = 11,ncol = 2)
FPR7_matrix[,2] <-non.contrib.ave.fpr_mat[,7]

col_vec <- c("Expected Heterozygosity","FPR")
colnames(FPR7_matrix) <- col_vec

rownames(FPR7_matrix) <- pop_vec2
FPR7_matrix[,1] <- exphetero_mat

FPR7_matrix



#This is the vector used to make the different colored lines for all of the populations 
colorvec = c("blue","red","green","purple","violet","yellow","pink","black","chartreuse",
             "coral","aquamarine4")


#This intiates the plot 
plot.new()
plot(FPR7_matrix,main ="FPR vs Expected Heterozygosity (7 contributors)",ylab = "FPR", xlab = 
       "Expected Heterozygosity")


col_vec <- c("Expected Heterozygosity","FPR")
colnames(fprcount_mat) <- 1:8

rownames(fprcount_mat) <- pop_vec2
FPR7_matrix[,1] <- exphetero_mat

FPR7_matrix

write.csv(fprcount_mat,"FPR_allpop.csv")
