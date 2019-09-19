########################## VARIABLES THAT NEED TO BE INITIATED####################################################
#The files that will be used 
pop_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
             "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")

#The names of populations to add to the rows of the matrix 
pop_vec2 <- c("African American","Apache","Bahamian","Caucasian","Chamorro","Filipino",
              "Jamaican","Navajo","SE Hispanic","SW Hispanic","Trinidadian")

#num of contribs
contributer = 1:8



###############################Total contrib function MATT's CODE ##############################################################
library(forensim)

total_contrib_function_eval_pop <- function(file_name,eval_pop,num_contribs,iterations,non_or_true_contrib) 
{
  set.seed(657)
  
  LR_vector_true <- c()
  total_LR <- c()
  
  allelefreq_tbl <- read.csv(file_name)    
  pop_freq <- tabfreq(tab=allelefreq_tbl, pop.names=as.factor('population'))
  
  allelefreq_tbl_eval <- read.csv(eval_pop)    
  pop_freq_eval <- tabfreq(tab=allelefreq_tbl_eval, pop.names=as.factor('population'))
  
  popgeno_eval <- simugeno(tab=pop_freq_eval, n=1, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                                 "TPOX", "D3S1358", "D7S820", "D16S539",
                                                                 "D18S51", "D21S11", "D2S1338", "D19S433","vWA"))
  
  if (non_or_true_contrib == 0){
    j = 1         
    while (j < iterations){
      
      popgeno <- simugeno(tab=pop_freq, n=num_contribs, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                                      "TPOX", "D3S1358", "D7S820", "D16S539",
                                                                      "D18S51", "D21S11", "D2S1338", "D19S433","vWA"))
      
      mix.func <- simumix(popgeno,ncontri = num_contribs)
      
      LR_vector <- c()
      total_LR_vector <- c()
      
      if (num_contribs == 1){
        
        k = 1
        while (k < 14){
          single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                          Tp=c(as.numeric(strsplit(mix.func$mix.prof[1,k], "/")[[1]])),  
                          Td=0,Vp=0,Vd=c(as.numeric(strsplit(mix.func$mix.prof[1,k], "/")[[1]])),
                          xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                          prDHom=c(0.04,0.04),prC=0,freq=popgeno_eval@tab.freq[[1]][[k]])
          
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
        
      } else {  
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
                          prDHom=c(0.04,0.04),prC=0,freq=popgeno_eval@tab.freq[[1]][[k]])
          while(is.na(single_LR$LR) #|| #is.infinite(single_LR$LR) 
                || is.nan(single_LR$LR) || single_LR$LR<0) {
            single_LR <- LR(Repliste=c(mix.func@mix.all[[k]]),
                            Tp=c(listOfAlleles),
                            Td=c(listofAllels.Defense), ## Vd = Non contrib under Hd - Suspect ##
                            Vp=0,Vd=c(as.numeric(tail(listOfAlleles,2), "/")[[1]]), ## Vp = Non contrib under Hp - 0 ##
                            xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                            prDHom=c(0.04,0.04),prC=0,freq=popgeno_eval@tab.freq[[1]][[k]])
          } 
          #print(popgeno_eval@tab.freq[[1]][[k]])
          #print(single_LR)
          #print(listOfAlleles)
          #print(listofAllels.Defense)
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      }
      LR_vector_true[j] <- prod(LR_vector)
      
      j = j + 1
      total_LR_vector <- c(log10(LR_vector_true))
      
    }
  } else if (non_or_true_contrib == 1) {
    j = 1
    while (j < iterations){ #loop over simulations
     # print(paste("on iteration number", j))
      popgeno <- simugeno(tab=pop_freq, n=num_contribs, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                                      "TPOX", "D3S1358", "D7S820", "D16S539",
                                                                      "D18S51", "D21S11", "D2S1338", "D19S433","vWA"))
      
      noncontrib <- simugeno(tab=pop_freq, n=1, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                              "TPOX", "D3S1358", "D7S820", "D16S539",
                                                              "D18S51", "D21S11", "D2S1338", "D19S433",
                                                              "vWA"))
      noncon.sus <- simumix(noncontrib, ncontri = 1)
      
      
      mix.func <- simumix(popgeno,ncontri = num_contribs)
      #print("made mixture")
      LR_vector <- c()  
      total_LR_vector <- c()
      
      if (num_contribs == 1){
        k = 1
        while (k < 14){
          single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                          Tp=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),  
                          Td=0,Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                          xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                          prDHom=c(0.04,0.04),prC=0,freq=popgeno_eval@tab.freq[[1]][[k]])
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
                          prDHom=c(0.04,0.04),prC=0,freq=popgeno_eval@tab.freq[[1]][[k]])
          while(is.na(single_LR$LR) #|| #is.infinite(single_LR$LR) 
                || is.nan(single_LR$LR) || single_LR$LR<0) {
            single_LR <- LR(Repliste=c(mix.func$mix.all[[k]]),                         
                            Tp=c(listOfAlleles,  
                                 as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])), 
                            Td=c(listofAllels.Defense),   
                            Vp=0,Vd=c(as.numeric(strsplit(noncon.sus$mix.prof[1,k], "/")[[1]])),
                            xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                            prDHom=c(0.04,0.04),prC=0,freq=popgeno_eval@tab.freq[[1]][[k]])
            print("in the while crappy forensim result check")
          }
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      }
      LR_vector_true[j] <- prod(LR_vector)
      
      total_LR_vector <- c(log10(LR_vector_true))
      #total_LR_vector[total_LR_vector < 0.0] <- -1
      j = j + 1
    }
  }
  return(total_LR_vector)
}
LRs = total_contrib_function_eval_pop('Filipino.csv','swh.csv', 6, 10, 0)





#################m STEP:1 Monday 3/11 solo Makes a 3d array for TRUE CONTRIBUTORS ################################################
#STEP:1 Is the code that makes an empty 3d arratthen passes each population to the total contrib function 
#(redone on May21, 2019 using updated code)



# need to initalize for the forensim function to work
total_LR_vector <- c()
 

#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
array_3d = array(rep(1:1000,20*8*11), dim =c(1000,8,11))

#view empty array
array_3d 

popsize_vec = 2:9


#This loop uses and x,y,z variable scheme to represent locations in 3d space 
# x is the number of contributors [,x,]
# y are the populations [,,y]
# z are the slices of 1000 iterations [z,,]

for (y in 1:length(pop_vec)) {
  print(pop_vec[y])
  for (x in 1:length(contributer)) {
    result = total_contrib_function_eval_pop(pop_vec[y], pop_vec[y], x, 1001, 0)
    
      
      # total_LR_plot_vectornav6 is the gobal vaiable that is created from the forensim function that has the vector of LR's 
      # total_LR_plot_vectornav6 is only for true contributers
      #result = total_LR_vector
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

array_3dnoncontrib = array(rep(1,1000*8*11), dim =c(1000,8,11))
#view empty array

array_3dnoncontrib



ptm <- proc.time()
for (y in 1:length(pop_vec)) {
  print(pop_vec[y])
  
  for (x in 1:length(contributer)) {
   result = total_contrib_function_eval_pop (pop_vec[y], pop_vec[y], x, 1001, 1)
    
    # total_LR_plot_vectornav6 is the gobal vaiable that is created from the forensim function that has the vector of LR's 
    # total_LR_plot_vectornav6 is only for true contributers
    #result = total_LR_plot_vector_nonnav8
    #print(result)
    array_3dnoncontrib[,x,y] = result
   
  }
}

proc.time() - ptm
array_3dnoncontrib

########################################## STEP:2 Loop that counts FPR greater than -1 for the entire 3d array #######################################################################
## This loop counts the number of FPRs found in each population per contributor and stores it in the matrix FPRcount

fprcount_mat = matrix(1:77,nrow = 11,ncol = 7)

for (i in 1:11){
  population = pop_vec[i]
  print(pop_vec[i])
  
  for(j in 1:7){
    contrib = j
    count=1
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:1000){
      sticks <<- array_3dnoncontrib[k,j,i]
      
      
      if(is.nan(sticks)){
        count = count
        print(paste("nan found at", k,j,i))
      
      }
      else if(sticks > -1){
        count = count + 1
       
      }
      
    }
    
    result = (count - 1)/1000
    #result = (count - 1)
    #print(result)
    fprcount_mat[i,j] = result
    
  }
}

non_fprcount_mat <- fprcount_mat



#This is all naming schemes for the matrix row and columns 

iterations = 1:7
colnames(non_fprcount_mat) <- iterations
rownames(non_fprcount_mat) <- pop_vec2
non_fprcount_mat

write.csv(non_fprcount_mat, "noncontrib_mat2.csv")

#################################### STEP 2.1 Loop that counts FPR greater than -1 for the entire 3d array #####################################################

#### SAME LOOP FOR TRUE CONTRIBUTORS 

fprcount_mat = matrix(1:77,nrow = 11,ncol = 7)

for (i in 1:11){
  population = pop_vec[i]
  print(pop_vec[i])
  
  for(j in 1:7){
    contrib = j
    count=1
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:1000){
      sticks = array_3d[k,j,i]
      
      
      if(is.nan(sticks)){
        count = count
        print(paste("nan found at", k,j,i))
        
      }
      else if(sticks > -1){
        count = count + 1
        
      }
      
    }
    
    result = (count - 1)/1000
    #result = (count - 1)
    #print(result)
    fprcount_mat[i,j] = result
    
  }
}


true_fprcount_mat <- fprcount_mat


#This is all naming schemes for the matrix row and columns 
iterations = 1:7
colnames(true_fprcount_mat) <- iterations
rownames(true_fprcount_mat) <- pop_vec2
true_fprcount_mat

#Write the Matrix to a CSV
write.csv(true_fprcount_mat, "truecontrib_mat2.csv")


################################# STEP:3 PLOT for  FPR  counts of TRUE & NON CONTRIBUTORS#####################################################################
#STEP 3 takes the non contrib and true contrib matrixes and puts them on one plot 


#This is the vector used to make the different colored lines for all of the populations 
colorvec = c( "blue","red","green","purple","violet","yellow","pink","black")

plot.new()

#This intiates the plot 
plot(x=c(),y=c(),xlim = c(1,7),ylim = c(0,.06),main =" False Positives  (1000 iterations)",ylab = "false postive", xlab = 
       "number of contributors")


#this loop adds points onto the intiated plot with i reffering to the location in the matrix 
for(i in 1:7){
  points(non_fprcount_mat[i,], col= colorvec[i], lty=1)
  #points(true_fprcount_mat[i,], col=colorvec[i],lty=1)
}

legend("center", legend=c("AA","Apache","Bahamian","White","Chamorro","Filipino",
                          "Jamaican","Navajo","SEH","SWH","Trinidadian"),
       col=c( "blue","red","green","purple","violet","yellow","pink","black"), lty=1:2, cex=0.8)
####
qplot(data = test_fpr, xlab = "Contributors", ylab = "Number of False Postives") + scale_color_brewer(palette = "Accent", ylab("Population"))

col = viridis(10))
########################################STEP:4 LOOP that calculates the AVERAGE FPR carina 3/20/19############################################################
#STEP 2: This is done after LRs have been calculated for all populations and have been put into a 3darray.
#this is an empty matrix that looks through the 3d array and counts the instances where the LR is greater than -1
# the loop then divides the # of LR by the amount of iterations and it is stored in a Matrix that will be used for plotting 

fpraverage_mat = matrix(1:88,nrow = 11,ncol = 8)

for (i in 1:11){
  population = pop_vec[i]
  print(pop_vec[i])
  
  for(j in 1:8){
    contrib = j
    count=1
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:1000){
      #this needs to be changed depending on which array is being used 
      sticks = array_3dnon_big[k,j,i]
      
      
      if(sticks > -1){
        
        count = count + 1
        print(count)
        
      }
      
      
    }
    
    #this needs to be changed depending on how many iterations are used 
    result = (count - 1)/1000
    print(result)
    fpraverage_mat[i,j] = result
    
  }
}

true_fpraverage_mat <-fpraverage_mat
non_fpraverage_mat <- fpraverage_mat
################################# STEP:5  PLOT for AVerage FPR  counts of TRUE & NON CONTRIBUTORS#####################################################################
#STEP 3 takes the non contrib and true contrib matrixes and puts them on one plot 


#This is the vector used to make the different colored lines for all of the populations 
colorvec = c("blue","red","green","purple","violet","yellow","pink","black")



#This intiates the plot 
plot.new()
plot(x=c(),y=c(),xlim = c(1,8),ylim = c(0,1000),main =" False Positives  (1000 iterations)",ylab = "false postive", xlab = 
       "number of contributors")


#this loop adds points onto the intiated plot with i reffering to the location in the matrix 
for(i in 1:8){
  points(non_fpraverage_mat[i,], col= colorvec[i], lty=1)
 # points(true_fpraverage_mat[i,], col=colorvec[i],lty=1)
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
FPR7_matrix[,2] <-non_fprcount_mat[,6]

col_vec <- c("Expected Heterozygosity","FPR")
colnames(FPR7_matrix) <- col_vec

rownames(FPR7_matrix) <- pop_vec2
FPR7_matrix[,1] <- exphetero_mat

FPR7_matrix

FPR7_df <- as.data.frame.table(FPR7_matrix)

FPR7_df

#This is the vector used to make the different colored lines for all of the populations 
colorvec = c("blue","red","green","purple","violet","yellow","pink","black","chartreuse",
             "coral","aquamarine4")


#This intiates the plot 
plot.new()
plot(FPR7_matrix,main ="FPR vs Expected Heterozygosity",ylab = "False Postive Rate", xlab = 
       "Expected Heterozygosity")


FPR7_plot <- qplot(data = FPR7_, x = Var2, y = Freq, color = Var1, xlab = "Number of Contributors", ylab = "False Postive Rate", lwd = 5) + scale_color_brewer(palette = "Spectral", ylab("Population")) 


col_vec <- c("Expected Heterozygosity","FPR")
colnames(fprcount_mat) <- 1:8

rownames(fprcount_mat) <- pop_vec2
FPR7_matrix[,1] <- exphetero_mat

FPR7_matrix

write.csv(fprcount_mat,"FPR_allpop6contributors.csv")
#################################### STEP 7 FALSE # of contributors vs false postive rates#########################
#This graph is plot all 8 contribs vs. FPRs for Apache, African American, and SW Hispanic 
popvec_3 = c("African American","Apache", "SE Hispanic")

FPR_3pop_matrix <-matrix(1:8,nrow = 3,ncol = 8)
FPR_3pop_matrix[1,] <-non_fprcount_mat[1,]
FPR_3pop_matrix[2,] <-non_fprcount_mat[2,]
FPR_3pop_matrix[3,] <-non_fprcount_mat[9,]

for(i in FPR_3pop_matrix[])

col_vec <- c(1:8)
colnames(FPR_3pop_matrix) <- col_vec

rownames(FPR_3pop_matrix) <- popvec_3
FPR_3pop_matrix[,1] <- exphetero_mat

FPR_3pop_matrix


###### poster all pop table that looks at 7 contributors ##################################

FPR_allpop_tbl <- as.data.frame.table(non_fprcount_mat)

FPR_allpop_tbl_plot <- qplot(data = FPR_allpop_tbl, x = Var2, y = Freq, color = Var1, lwd = 3,xlab = "Number of Contributors", ylab = "False Postive Rate") + scale_color_brewer(palette = "Spectral", ylab("Population")) 
#### Removes the grey in the background 
FPR_allpop_plot_2 <- FPR_allpop_tbl_plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

###############################################################################################################