########################## VARIABLES THAT NEED TO BE INITIATED####################################################
#The files that will be used 

#pop_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv", "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")

pop_vec <- clean_popnames

#The names of populations to add to the rows of the matrix 
pop_vec2 <- c("African American","Apache","Bahamian","Caucasian","Chamorro","Filipino",
              "Jamaican","Navajo","SE Hispanic","SW Hispanic","Trinidadian")

loci_vec <- c("CSF1PO","D3S1358","D5S818","D7S820","D8S1179","D13S317","D16S539","D18S51","D21S11",
              "FGA","TH01","TPOX","VWA")
#num of contribs
contributer = 1:8



###############################Total contrib function MATT's CODE ##############################################################
library(forensim)

total_contrib_function <- function(file_name,pop_size,num_contribs,iterations,non_or_true_contrib) 
{
  set.seed(657)
  allelefreq_tbl <- read.csv(file_name)    
  pop_freq <- tabfreq(tab=allelefreq_tbl, pop.names=as.factor('population'))
  popgeno <- simugeno(tab=pop_freq, n=pop_size, which.loc = c("CSF1PO","D3S1358","D5S818","D7S820","D8S1179",
                                                              "D13S317","D16S539","D18S51","D21S11",
                                                              "FGA","TH01","TPOX","VWA"))
  
  
  noncontrib <- simugeno(tab=pop_freq, n=1, which.loc = c("CSF1PO","D3S1358","D5S818","D7S820","D8S1179",
                                                          "D13S317","D16S539","D18S51","D21S11",
                                                          "FGA","TH01","TPOX","VWA"))
  #print(afrigeno@tab.geno)
  noncon.sus <- simumix(noncontrib, ncontri = 1)
  
  
                                                 ## Inititialize variables for vector ##
  csf.2 <- pop_freq$tab$population$CSF1PO
  d3s.2 <- pop_freq$tab$population$D3S1358
  d2s.2 <- pop_freq$tab$population$D5S818         #D2S1338 *changed out for a diff STR
  d7s.2 <- pop_freq$tab$population$D7S820
  d8s.2 <- pop_freq$tab$population$D8S1179
  d16.2 <- pop_freq$tab$population$D16S539
  d18.2 <- pop_freq$tab$population$D18S51
  d21.2 <- pop_freq$tab$population$D21S11
  fga.2 <- pop_freq$tab$population$FGA          
  th01.2 <- pop_freq$tab$population$TH01
  tp0x.2 <- pop_freq$tab$population$TPOX
  vwa.2 <- pop_freq$tab$population$VWA
  d19.2 <- pop_freq$tab$population$D13S317       #D19S433 *
  
  freq_vector <- matrix(list(), nrow=13, ncol=1)  ## Creating Vector ##
  freq_vector[[1,1]] <- c(d8s.2)
  freq_vector[[2,1]] <- c(th01.2)
  freq_vector[[3,1]] <- c(fga.2)
  freq_vector[[4,1]] <- c(csf.2)
  freq_vector[[5,1]] <- c(tp0x.2)
  freq_vector[[6,1]] <- c(d3s.2)
  freq_vector[[7,1]] <- c(d7s.2)
  freq_vector[[8,1]] <- c(d16.2)
  freq_vector[[9,1]] <- c(d18.2)
  freq_vector[[10,1]] <- c(d21.2)
  freq_vector[[11,1]] <- c(d2s.2)
  freq_vector[[12,1]] <- c(d19.2)
  freq_vector[[13,1]] <- c(vwa.2)
  
  LR_vector_true <- c()
  total_LR <- c()
  plot_vector <- c()
  
  if (non_or_true_contrib == 0){                                      #loop 1: 
    
    j = 1         
    while (j < iterations){
      mix.func <- simumix(popgeno,ncontri = num_contribs)
      
      LR_vector <- c()
      total_LR_vector <- c()
      print((paste("in loop 1")))
      
      if (num_contribs == 1){                                     #loop 2:
        print(paste("top of loop 2"))
        k = 1
        while (k < 14){
          single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                          Tp=c(as.numeric(strsplit(mix.func$mix.prof[1,k], "/")[[1]])),  
                          Td=0,Vp=0,Vd=c(as.numeric(strsplit(mix.func$mix.prof[1,k], "/")[[1]])),
                          xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                          prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
          #print(single_LR$LR)
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
          print(paste("in loop 2"))
        }
        
      } else {  #if more than one contributor 
        k = 1                             ## LR Calculation. ##
        while (k < 14){  
          ## Calculated 13 times for 13 loci ##
          c = num_contribs - 1  
          listOfAlleles = c()
          for(i in 1:num_contribs){
            listOfAlleles = c(listOfAlleles, as.numeric(strsplit(mix.func$mix.prof[i,k], "/")[[1]]))
            listofAllels.Defense = c()
            for (i in 1:c){
              listofAllels.Defense = c(listofAllels.Defense, as.numeric(strsplit(mix.func$mix.prof[i,k], "/")[[1]]))   
              
            }
          }
          single_LR <- LR(Repliste=c(mix.func@mix.all[[k]]),
                          Tp=c(listOfAlleles),
                          Td=c(listofAllels.Defense), ## Vd = Non contrib under Hd - Suspect ##
                          Vp=0,Vd=c(as.numeric(tail(listOfAlleles,2), "/")[[1]]), ## Vp = Non contrib under Hp - 0 ##
                          xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                          prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
          while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0) {
            single_LR <- LR(Repliste=c(mix.func@mix.all[[k]]),
                            Tp=c(listOfAlleles),
                            Td=c(listofAllels.Defense), ## Vd = Non contrib under Hd - Suspect ##
                            Vp=0,Vd=c(as.numeric(tail(listOfAlleles,2), "/")[[1]]), ## Vp = Non contrib under Hp - 0 ##
                            xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                            prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
            
            
          
            
          } 
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      }
      LR_vector_true[j] <- prod(LR_vector)
      
      j = j + 1
      total_LR_vector <- c(log10(LR_vector_true))
      total_LR_plot_vectornav6 <<- total_LR_vector 
      g = 1
      while (g < 3){
        allelefreq_tbl <- read.csv(file_name)    
        pop_freq <- tabfreq(tab=allelefreq_tbl, pop.names=as.factor('Afri'))
        popgeno <- simugeno(tab=pop_freq, n=pop_size, which.loc = c("CSF1PO","D3S1358","D5S818","D7S820",
                            "D8S1179","D13S317","D16S539","D18S51","D21S11",
                            "FGA","TH01","TPOX","VWA"))
                            
        g = g + 1
      }
      
    }
  } else if (non_or_true_contrib == 1) {
    j = 1
    #set.seed(657)
    while (j < iterations){ #loop over simulations
      mix.func <- simumix(popgeno,ncontri = num_contribs)
      
      LR_vector <- c()  
      total_LR_vector <- c()
      
      if (num_contribs == 1){
        k = 1
        while (k < 14){
          single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                          Tp=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),  
                          Td=0,Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                          xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                          prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
          #print(listOfAlleles)
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      } else {
        k = 1
        while (k < 14){  
          c = num_contribs - 1
          listOfAlleles = c()
          for(i in 1:num_contribs){
            listOfAlleles = c(listOfAlleles, as.numeric(strsplit(mix.func$mix.prof[i,k], "/")[[1]]))
            listofAllels.Defense = c()
            for (i in 1:c){
              listofAllels.Defense = c(listofAllels.Defense, as.numeric(strsplit(mix.func$mix.prof[i,k], "/")[[1]]))
              
              
            }
          }  
          
            single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                            Tp=c(listOfAlleles,  
                                 as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])), 
                            Td=c(listofAllels.Defense),   
                            Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                            xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                            prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
            while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0) {
              single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                              Tp=c(listOfAlleles,  
                                   as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])), 
                              Td=c(listofAllels.Defense),   
                              Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                              xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                              prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
            }
            
            
            
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      }     
      total_LR[j] <- (LR_vector[1]*LR_vector[2]*LR_vector[3]*
                        LR_vector[4]*LR_vector[5]*LR_vector[6]*
                        LR_vector[7]*LR_vector[8]*LR_vector[9]*
                        LR_vector[10]*LR_vector[11]*LR_vector[12]*
                        LR_vector[13])
      
      total_LR_vector <- c(log10(total_LR))
      
      total_LR_vector[total_LR_vector < 0.0] <- -1
      
      total_LR_plot_vector_nonnav8 <<- total_LR_vector
      j = j + 1
    }
  }
  return(total_LR_vector)
}


total_contrib_function <- function(file_name,pop_size,num_contribs,iterations,non_or_true_contrib)   
  
LRs = total_contrib_function("USAf.1_new.csv", 100, 1, 20, 0)


#################m STEP:1 Monday 3/11 solo Makes a 3d array for TRUE CONTRIBUTORS ################################################
#STEP:1 Is the code that makes an empty 3d arratthen passes each population to the total contrib function 




# need to initalize for the forensim function to work
total_LR_vector <- c()
total_LR_plot_vectornav6 <- c()

#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
array_3d = array(rep(1,1000*8*11), dim =c(1000,8,11))

#view empty array
array_3d 


#This loop uses and x,y,z variable scheme to represent locations in 3d space 
# x is the number of contributors [,x,]
# y are the populations [,,y]
# z are the slices of 1000 iterations [z,,]
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
total_LR_vector <- c()
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
    #sticks is iterations 
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


plot(exphetero_mat_clean)
