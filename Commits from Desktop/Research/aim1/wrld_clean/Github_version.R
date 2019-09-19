#notes for next time: 

# do we need a pop size? 


########################## VARIABLES THAT NEED TO BE INITIATED####################################################
#The files that will be used 

#pop_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv", "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")

pop_vec <- clean_popnames

#The names of populations to add to the rows of the matrix 
pop_vec2 <- c("African American",
              "Apache",
              "Bahamian",
              "Caucasian",
              "Chamorro",
              "Filipino",
              "Jamaican",
              "Navajo",
              "SE Hispanic",
              "SW Hispanic",
              "Trinidadian"
              )

loci_vec <- c("CSF1PO",
              "D3S1358",
              "D5S818",
              "D7S820","
              D8S1179",
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

#########################################
## Ignore, For Debugging Function Only ##
#########################################

#file.name = "USAf.1_new.csv"
#pop.size = 100
#num_contribs = 4
#num_sims = 1000
#non.or.truecontrib = 0


###########################
#### START OF FUNCTION ####
###########################
total_contrib_function <- function(file.name, pop.size, num_contribs, num_sims, non.or.truecontrib){
  # set.seed(657) 
  # We have to use this seed for reproduceability 
  set.seed (123560) 
  
  #reads in the files we are working with 
  afs_csv <- read.csv(file.name)    
  
  #####tabfreq makes an onject that holds the following information###
  # @TAB function reads in all the allele frequencies per loci 
  # @which.loc reports the loci that are taken into consderation
  # @pop.names are right now just reads "population" but will try to change it to the actual population (Changed)
  #you can access indidual STR freq data like this: pop.afs$tab$USAf.1_new.csv$VWA
  
  pop.afs <- tabfreq(tab = afs_csv, 
                     pop.names = as.factor('population')
                     )
  
  #######simugeno objects store genotypes from the tabfreq ###########
  #popgen$tab.geno gives the genotyprs of all individuals (n) 
  #look at a particular indivduals genotype: sim.genotypes$tab.geno[1:10,]
  

  sim.genotypes <- simugeno(tab = pop.afs, 
                            n = pop.size, 
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
 
  ### QUESTION:figure out the population situation they are diffrent depending on the population size 
  #i.e I think the population size matters because when it only generates 3 people it uses just those 3
  #when it geenrates 100 people other indivudals that are simulated can be used but I don't know if this matters 
 
  ###this does nothing right now
  #noncontrib <- simugeno(tab = pop.afs, 
                         #n = pop.size, 
                         #which.loc = c("CSF1PO",
                                       #"D3S1358",
                                       #"D5S818",
                                       #"D7S820",
                                       #"D8S1179",
                                       #"D13S317",
                                       #"D16S539",
                                       #"D18S51",
                                       #"D21S11",
                                       #"FGA",
                                       #"TH01",
                                       #"TPOX",
                                       #"VWA"
                                       #)
                         #)
  
  #SIMUMIX objects store DNA mixtures#########
  ####noncontrib is the simulated genotype data
  #### ncontri is the number of contributors to the mixture 
  #this use to use noncontrib but I didn't see how it was diff from sim. genotypes
  noncon.sus <- simumix(sim.genotypes, 
                        ncontri = num_contribs
                        )
  
  ###### QUESTION: Why do we need this here if we make the mixutres inside the loop? probably two diff populations 
  #commenting out for now 
  
  ## Inititialize variables for vector ##
  csf.2 <- pop.afs$tab$population$CSF1PO
  d3s.2 <- pop.afs$tab$population$D3S1358
  #D2S1338 *changed out for a diff STR
  d2s.2 <- pop.afs$tab$population$D5S818         
  d7s.2 <- pop.afs$tab$population$D7S820
  d8s.2 <- pop.afs$tab$population$D8S1179
  d16.2 <- pop.afs$tab$population$D16S539
  d18.2 <- pop.afs$tab$population$D18S51
  d21.2 <- pop.afs$tab$population$D21S11
  fga.2 <- pop.afs$tab$population$FGA          
  th01.2 <- pop.afs$tab$population$TH01
  tp0x.2 <- pop.afs$tab$population$TPOX
  vwa.2 <- pop.afs$tab$population$VWA
  #D19S433 *
  d19.2 <- pop.afs$tab$population$D13S317       
  
  ## Creating Matrix ## Numeric, 45? Not sure this is working 
  pop.afs.matrix <- matrix(list(), 
                           nrow = 13, 
                           ncol = 1
                           ) 
  
  pop.afs.matrix[[1,1]] <- c(d8s.2)
  pop.afs.matrix[[2,1]] <- c(th01.2)
  pop.afs.matrix[[3,1]] <- c(fga.2)
  pop.afs.matrix[[4,1]] <- c(csf.2)
  pop.afs.matrix[[5,1]] <- c(tp0x.2)
  pop.afs.matrix[[6,1]] <- c(d3s.2)
  pop.afs.matrix[[7,1]] <- c(d7s.2)
  pop.afs.matrix[[8,1]] <- c(d16.2)
  pop.afs.matrix[[9,1]] <- c(d18.2)
  pop.afs.matrix[[10,1]] <- c(d21.2)
  pop.afs.matrix[[11,1]] <- c(d2s.2)
  pop.afs.matrix[[12,1]] <- c(d19.2)
  pop.afs.matrix[[13,1]] <- c(vwa.2)
  

  #total_LR <- c()
  #LR_vector <- c()
  log10_LR <- c(1:num_sims)
  
  
  #  If non.or.truecontrib == 0 then we are simulationing a mixture where the 
  # suspect did not contribute DNA to the sample
  
  ##############################
  ### true contributor code ####
  ##############################
  
  #is a true contributor mix   #loop 1:
  if (non.or.truecontrib == 0){
    j = 1 
    
    #I think this should be <=
    while (j < num_sims){
      
     #simulate a mixture using the simulated genotypes that consider this populations 
     #allele frequencies for the appropriate population size 
     sim.mix <- simumix(sim.genotypes,
                        ncontri = num_contribs
                        ) 
     singleLR_vector <- c() 
    
     #loop 2:
     if (num_contribs == 1){                                    
       
       k = 1
       #K stops at 14 because we are considering the 13 core STR loci
       while (k < 14){
         
          single_LR <- LR( Repliste = c(sim.mix$mix.all[[k]]),
                           Tp = c(as.numeric(strsplit(sim.mix$mix.prof[1,k], "/")[[1]])),
                           Td = 0,
                           Vp = 0,
                           Vd = c(as.numeric(strsplit(sim.mix$mix.prof[1,k], "/")[[1]])),
                           xd = 1,
                           xp = 0,
                           theta = 0,
                           prDHet = c(0.2,0.2),
                           prDHom = c(0.04,0.04),
                           prC = 0,
                           freq = pop.afs.matrix[[k,1]]
                          )
          
          #print(single_LR$LR)
          singleLR_vector <- c(singleLR_vector,single_LR$LR)
          k = k + 1
       }
     } 
     
     #if more than one contributor it goes here 
     else {  
       
       #loop 3## LR Calculation
       k = 1 
       while (k < 14){  
         prosecutor.all.atk = c()
         
         for(i in 1:num_contribs){
           
           ### QUESTION: why do we need to use strsplit ? Rename these alleles so they make sense 
           prosecutor.all.atk = c(prosecutor.all.atk, 
                                 as.numeric(strsplit(sim.mix$mix.prof[i,k], "/")[[1]])
                                 )
         } 
         defense.all.atk = c()
         #i think num_contribs should be the same as above (i in 1:(num_contribs - 1)
         for (i in 1:(num_contribs - 1)){
           defense.all.atk = c(defense.all.atk,
                               as.numeric(strsplit(sim.mix$mix.prof[i, k], "/")[[1]])
                               )
         }
         
         #Loop 4 # detail what goes in each var
         single_LR <<- LR( Repliste = c(sim.mix@mix.all[[k]]),
                           Tp = c(prosecutor.all.atk),
                           ## Vd = Non contrib under Hd - Suspect ##
                           Td = c(defense.all.atk),
                           #tells Dorothy to return the last two digits of vector (dont think we need the strsplit commands)
                           Vp = 0,
                           ## Vp = Non contrib under Hp - 0 ##
                           Vd = c(as.numeric(tail(prosecutor.all.atk,2), "/")[[1]]), 
                           xd = 1,
                           xp = 0,
                           theta = 0,
                           prDHet = c(0.2,0.2),
                           prDHom = c(0.04,0.04),
                           prC = 0,
                           freq = pop.afs.matrix[[k,1]]
                           )
          
        
#######################################
##### Not using this while loop #######        
######################################       
          #while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0){
           #single_LR <- LR( Repliste = c(sim.mix@mix.all[[k]]),
                             #Tp = c(prosecutor.all.atk),
                             #Td = c(defense.all.atk), 
                              ##Vd = Non contrib under Hd - Suspect ##
                            # Vp = 0,
                            ## Vp = Non contrib under Hp - 0 ##
                             #Vd = c(as.numeric(tail( prosecutor.all.atk,2), "/")[[1]]),
                             #xd = 1,
                            # xp = 0,
                           # theta = 0,
                            #prDHet = c(0.2,0.2),
                             #prDHom = c(0.04,0.04),
                             #prC = 0,
                            #freq = pop.afs.matrix[[k,1]]
                           # )
                        #print(paste("in stupid while loop", k))     
            
          
          #end of while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0)  
         # }
          
         print(single_LR$LR)
         if ((single_LR$LR) == "Inf"){
           single_LR$LR <- 1
           singleLR_vector <<- c(singleLR_vector,single_LR$LR)
            } else {
            singleLR_vector <<- c(singleLR_vector,single_LR$LR)
            }
         print(singleLR_vector)
          #singleLR_vector <<- c(singleLR_vector,single_LR$LR)
           k = k + 1
           
       #end of while (k < 14), for LR calculation   
       }
       
     #end of else, for more than 1 contributor   
     }
     
     print(prod(singleLR_vector))
     print(log10(prod(singleLR_vector)))
     
    ### Q:annot find log10 is not found ????
     print(log10_LR[j] <<- log10(prod(singleLR_vector)))
      j = j + 1
      
    #end of while (j < num_sims)
    }
    
  #end of if (non.or.truecontrib == 0)
  }
  
  ##########################
  ### non contributors #####
  ##########################
  
  else if (non.or.truecontrib == 1){     
    j = 1
    
    #set.seed(657)
    
    #loop over num_sims
    while (j < num_sims){ 
      
      sim.mix <- simumix(sim.genotypes,
                         ncontri = num_contribs
                         )
      #Q344#I think we should make this person here, but need clarity on what it should be
      #noncon.sus <- simumix(noncontrib, 
                            #ncontri = (num_contribs) 
     # )
      
      
      singleLR_vector <- c()  
      #log10_LRvector <- c()
      
      if (num_contribs == 1){
        k = 1
        while (k < 14){
          single_LR <- LR(Repliste = c(sim.mix$mix.all[[k]]),                         
                          Tp = c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),  
                          Td = 0,
                          Vp = 0,
                          Vd = c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                          xd = 1,
                          xp = 0,
                          theta = 0,
                          prDHet = c(0.2,0.2),
                          prDHom = c(0.04,0.04),
                          prC = 0,
                          freq = pop.afs.matrix[[k,1]]
                          )
          
          #print(list_Alleles)
          singleLR_vector <- c(singleLR_vector,single_LR$LR)
          k = k + 1
         
        }
        
      } 
      
      else {
        k = 1
        
        while (k < 14){  
         #9/12 c = num_contribs - 1
          prosecutor.all.atk = c()
          
          
          ################################################################
          ## I changed the { } for these loops so that for(i in 1:c) is ## 
          ## not included in loop for(i in 1:num_contribs). this was a  ##
          ## issue rori mentioned at meeting 8/30 - CK :)               ##
          ################################################################
          for(i in 1:num_contribs){
            prosecutor.all.atk = c(prosecutor.all.atk, 
                                   as.numeric(strsplit(sim.mix$mix.prof[i,k], 
                                                       "/")[[1]])
                                   )
          }
            defense.all.atk = c()
            
            for (i in 1:(num_contribs - 1)){
              defense.all.atk = c(defense.all.atk,
                                  as.numeric(strsplit(sim.mix$mix.prof[i,k], 
                                                      "/")[[1]])
                                  )
            }
            
            single_LR <- LR( Repliste = c(sim.mix$mix.all[[k]]),
                             Tp = c(prosecutor.all.atk,
                                    as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                             Td = c(defense.all.atk),
                             Vp = 0,
                             Vd = c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                             xd = 1,
                             xp = 0,
                             theta = 0,
                             prDHet = c(0.2,0.2),
                             prDHom = c(0.04,0.04),
                             prC = 0,
                             freq = pop.afs.matrix[[k,1]]
                             )

            #while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0){
                #single_LR <- LR( Repliste = c(sim.mix$mix.all[[k]]),                         
                                #Tp = c(prosecutor.all.atk, 
                                #as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                                #Td = c(defense.all.atk),
                               # Vp = 0,
                                #Vd = c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                                #xd = 1,
                                 #xp = 0,
                                 #theta = 0,
                                # prDHet = c(0.2,0.2),
                                #prDHom = c(0.04,0.04),
                                #prC = 0,
                               # freq = pop.afs.matrix[[k,1]]
                                #)
           # }
            print(single_LR$LR)
            ####### replces inf with 1 ######
            if ((single_LR$LR) == "Inf"){
              single_LR$LR <- 1
              singleLR_vector <<- c(singleLR_vector,single_LR$LR)
            } else {
              singleLR_vector <<- c(singleLR_vector,single_LR$LR)
            }
            
               # singleLR_vector <- c(singleLR_vector,single_LR$LR)
                k = k + 1
                
          
        #end of  while (k < 14)  
        }
        
      #end of else, for more than 1 contributor  
      }
      
        log10_LR[j] <- log10(prod(singleLR_vector))
        j = j + 1
    
    #end of while (j < num_sims), which loops over num_sims
    }
    
  #end of else if (non.or.truecontrib == 1)  
  }
    
    
    return(log10_LR) 
  
#end of total_contrib_function <- function(file.name,pop.size,num_contribs,num_sims,non.or.truecontrib)  
}


#########################
#### END OF FUNCTION ####
#########################
lrs <- c()
lrs <- total_contrib_function("USAf.1_new.csv", 15, 4, 11, 0)



#######################################################
## I did not clean anything below this line - CK 9/6 ##
#######################################################


#########################################################################################################################
#################m STEP:1 Monday 3/11 solo Makes a 3d array for TRUE CONTRIBUTORS ################################################
#STEP:1 Is the code that makes an empty 3d arratthen passes each population to the total contrib function 






#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
array_3d = array(rep(1,1000*8*11), dim =c(1000,8,11))

#view empty array
array_3d 


#This loop uses and x,y,z variable scheme to represent locations in 3d space 
# x is the number of contributors [,x,]
# y are the populations [,,y]
# z are the slices of 1000 num_sims [z,,]
for (y in 1:length(pop_vec)) {
  
  for (x in 1:length(contributer)) {
    total_contrib_function(pop_vec[y], 1000, x, 1001, 0)
    
    # total_LR_plot_vectornav6 is the gobal vaiable that is created from the forensim function that has the vector of LR's 
    # total_LR_plot_vectornav6 is only for true contributers
    result = total_LR_plot_vectornav6
    #print(result)
    array_3d[,x,y] = result
  }
}

array_3d
######################################STEP:1 3d array for NON CONTRIBUTOR all populations 3/16/19########################################
# STEP 1: same as above but with NON CONTRIBUTOR code 


# need to initalize for the forensim function to work
log10_totalLRs_vec <- c()
total_LR_plot_vector_nonnav8<- c()

#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
array_3dnon_big10 = array(rep(1,10000*8*11), dim =c(10000,8,11))
#view empty array

array_3dnon_big 



ptm <- proc.time()
for (y in 1:length(pop_vec)) {
  
  for (x in 1:length(contributer)) {
    total_contrib_function(pop_vec[y], 100000, x, 10001, 1)
    
    # total_LR_plot_vectornav6 is the gobal vaiable that is created from the forensim function that has the vector of LR's 
    # total_LR_plot_vectornav6 is only for true contributers
    result = total_LR_plot_vector_nonnav8
    #print(result)
    array_3dnon_big10[,x,y] = result
    
  }
}

proc.time() - ptm
array_3dnon

########################################## STEP:2 Loop that counts FPR greater than -1 for the entire 3d array #######################################################################
## This loop counts the number of FPRs found in each population per contributor and stores it in the matrix FPRcount

fprcount_mat = matrix(1:88,nrow = 11,ncol = 8)

for (i in 1:11){
  population = pop_vec[i]
  
  for(j in 1:8){
    contrib = j
    count=1
    #sticks is num_sims 
    #sticks = array_3d[,j,i]
    for(k in 1:10000){
      sticks = array_3dnon_big10[k,j,i]
      
      
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

true_fprcount_mat <-fprcount_mat
non_fprcount_mat <- fprcount_mat


#This is all naming schemes for the matrix row and columns 
iterations = 1:8
colnames(true_fprcount_mat) <- iterations
rownames(true_fprcount_mat) <- pop_vec2
true_fprcount_mat

#Write the Matrix to a CSV
write.csv(true_fprcount_mat, "truecontrib_mat2.csv")

iterations = 1:8
colnames(non_fprcount_mat) <- iterations
rownames(non_fprcount_mat) <- pop_vec2
non_fprcount_mat

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

fpraverage_mat = matrix(1:88,nrow = 11,ncol = 8)

for (i in 1:11){
  population = pop_vec[i]
  
  for(j in 1:8){
    contrib = j
    count=1
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:10000){
      #this needs to be changed depending on which array is being used 
      sticks = array_3dnon_big[k,j,i]
      
      
      if(sticks > -1){
        
        count = count + 1
        print(count)
      }
      
      
    }
    
    #this needs to be changed depending on how many iterations are used 
    result = (count - 1)/10000
    print(result)
    fpraverage_mat[i,j] = result
    
  }
}

true_fpraverage_mat <-fpraverage_mat
non_fpraverage_mat <- fpraverage_mat
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
  points(non_fpraverage_mat[i,], col= colorvec[i], lty=1)
  points(true_fpraverage_mat[i,], col=colorvec[i],lty=1)
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
FPR7_matrix[,2] <-non_fprcount_mat[,7]

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
