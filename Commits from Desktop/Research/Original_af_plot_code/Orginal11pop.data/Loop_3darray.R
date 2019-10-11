########################## VARIABLES THAT NEED TO BE INITIATED####################################################
#The files that will be used 
pop_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
             "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")

#The names of populations to add to the rows of the matrix 
pop_vec2 <- c("African American","Apache","Bahamian","Caucasian","Chamorro","Filipino",
              "Jamaican","Navajo","SE Hispanic","SW Hispanic","Trinidadian")




###############################Total contrib function##############################################################
library(forensim)

total_contrib_function <- function(file_name,pop_size,num_contribs,iterations,non_or_true_contrib) 
{
  set.seed(657)
  afri.gene <- read.csv(file_name)    
  afripop <- tabfreq(tab=afri.gene, pop.names=as.factor('Afri'))
  afrigeno <- simugeno(tab=afripop, n=pop_size, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                              "TPOX", "D3S1358", "D7S820", "D16S539",
                                                              "D18S51", "D21S11", "D2S1338", "D19S433",
                                                              "vWA"))
  
  noncontrib_Afri <- simugeno(tab=afripop, n=1, which.loc = c("D8S1179", "TH01" ,"FGA", "CSF1PO",
                                                              "TPOX", "D3S1358", "D7S820", "D16S539",
                                                              "D18S51", "D21S11", "D2S1338", "D19S433",
                                                              "vWA"))
  
  noncon.sus <- simumix(noncontrib_Afri, ncontri = 1)
  
  
  
  fga.2 <- afripop$tab$Afri$FGA           ## Inititialize variables for vector ##
  csf.2 <- afripop$tab$Afri$CSF1PO
  th01.2 <- afripop$tab$Afri$TH01
  tp0x.2 <- afripop$tab$Afri$TPOX
  vwa.2 <- afripop$tab$Afri$vWA
  d3s.2 <- afripop$tab$Afri$D3S1358
  d7s.2 <- afripop$tab$Afri$D7S820
  d8s.2 <- afripop$tab$Afri$D8S1179
  d16.2 <- afripop$tab$Afri$D16S539
  d18.2 <- afripop$tab$Afri$D18S51
  d21.2 <- afripop$tab$Afri$D21S11
  d2s.2 <- afripop$tab$Afri$D2S1338
  d19.2 <- afripop$tab$Afri$D19S433
  
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
  
  if (non_or_true_contrib == 0){
    
    j = 1         
    while (j < iterations){
      mix.func <- simumix(afrigeno,ncontri = num_contribs)
      
      LR_vector <- c()
      total_LR_vector <- c()
      
      if (num_contribs == 1){
        
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
                          prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
          while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0) {
            single_LR <- LR(Repliste=c(mix.func@mix.all[[k]]),
                            Tp=c(listOfAlleles),
                            Td=c(listofAllels.Defense), ## Vd = Non contrib under Hd - Suspect ##
                            Vp=0,Vd=c(as.numeric(tail(listOfAlleles,2), "/")[[1]]), ## Vp = Non contrib under Hp - 0 ##
                            xd=1,xp=0,theta=0,prDHet=c(0.2,0.2),
                            prDHom=c(0.04,0.04),prC=0,freq=freq_vector[[k,1]])
            
            
            #LR_vector <- c(LR_vector,single_LR$LR)
            #if(single_LR$LR>0)
            #if(is.nan(single_LR$LR)==FALSE){
            #  print(paste("Lr for iteration ", j, "allele", k, "is not nan it is ", single_LR$LR))
            #break
            #LR_vector <- c(LR_vector,single_LR$LR)
            
          } 
          LR_vector <- c(LR_vector,single_LR$LR)
          k = k + 1
        }
      }
      LR_vector_true[j] <- prod(LR_vector)#((LR_vector[1])*(LR_vector[2])*(LR_vector[3])*
      #(LR_vector[4])*(LR_vector[5])*(LR_vector[6])*
      #(LR_vector[7])*(LR_vector[8])*(LR_vector[9])*
      #(LR_vector[10])*(LR_vector[11])*(LR_vector[12])*
      #(LR_vector[13]))
      
      j = j + 1
      total_LR_vector <- c(log10(LR_vector_true))
      total_LR_plot_vectornav6 <<- total_LR_vector 
    }
  } else if (non_or_true_contrib == 1) {
    j = 1
    #set.seed(657)
    while (j < iterations){ #loop over simulations
      mix.func <- simumix(afrigeno,ncontri = num_contribs)
      
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



LRs = total_contrib_function('Trinidadian.csv', 100, 3, 101, 1)
#################m STEP:1 Monday 3/11 solo Makes a 3d array for TRUE CONTRIBUTORS ################################################
#STEP:1 Is the code that makes an empty 3d arratthen passes each population to the total contrib function 

pop_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
                 "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")
#num of contribs
contributer = 1:8


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

pop_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
                 "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")
#num of contribs
contributer = 1:8


# need to initalize for the forensim function to work
total_LR_vector <- c()
total_LR_plot_vector_nonnav8<- c()

#setting up the empty 3d array, rows are pop, columns are # of contrib, iteate-1 is the sheets of the cube
array_3dnon_big10 = array(rep(1,10000*8*11), dim =c(10000,8,11))
#view empty array

array_3dnon_big 
??byrow

#goes th
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



########################################STEP:4 Final code that counts FPR carina 3/20/19############################################################
#STEP 2: This is done after LRs have been calculated for all populations and have been put into a 3darray.
#this is an empty matrix that looks through the 3d array and counts the instances where the LR is greater than -1
# the loop then divides the # of LR by the amount of iterations and it is stored in a Matrix that will be used for plotting 

falsepositive_mat = matrix(1:88,nrow = 11,ncol = 8)

for (i in 1:11){
  population = pop_vec[i]
  
  for(j in 1:8){
    contrib = j
    count=1
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:10000){
      sticks = array_3dnon_big[k,j,i]
      
      
      if(sticks > -1){
        
        count = count + 1
        print(count)
      }
      
      
    }
    
    
    result = (count - 1)/10000
    print(result)
    falsepositive_mat[i,j] = result
    
  }
}

true_fpr_mat <-falsepositive_mat
non_fpr_mat <- falsepositive_mat
