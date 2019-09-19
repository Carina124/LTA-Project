########################################STEP:2 Final code that counts FPR carina 3/20/19############################################################
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
############################### STEP:3 Plot for False Postive Rates after producing the Matrix above ######################################
# STEP:3 this code is used to make plots after  the falsepostive matrix has been produced.
# this code will create plots of the FPR all populations 


#The files that will be used 
pop_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
                 "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")

#The names of populations to add to the rows of the matrix 
pop_vec2 <- c("African American","Apache","Bahamian","Caucasian","Chamorro","Filipino",
                  "Jamaican","Navajo","SE Hispanic","SW Hispanic","Trinidadian")

#Starts a PDF so the plots can be stored 
pdf("FPR_noncontrib_plot", width=20, height=20)

#This is how many plots per row,column 
par(mfrow=c(1,1))

######################################################################################################################

#This is all naming schemes for the matrix row and columns 
iterations = 1:8
colnames(true_fpr_mat) <- iterations
rownames(true_fpr_mat) <- pop_vec2
true_fpr_mat

#Write the Matrix to a CSV
write.csv(true_fpr_mat, "truecontrib_mat2.csv")

iterations = 1:8
colnames(non_fpr_mat) <- iterations
rownames(non_fpr_mat) <- pop_vec2
non_fpr_mat

write.csv(non_fpr_mat, "noncontrib_mat2.csv")
############################################## PLOT for AVerage FPR #####################################################################

#This is the vector used to make the different colored lines for all of the populations 
colorvec = c( "blue","red","green","purple","violet","yellow","pink","black")
symbol_vec = c('o',"+","*",".","-","s","n",",")
length(symbol_vec)

plot.new()

#This intiates the plot 
plot(x=c(),y=c(),xlim = c(1,8),ylim = c(0,1),main ="Average False Positive Rate (1000 iterations)",ylab = "average false postive rate", xlab = 
       "number of contributors")


#this loop adds points onto the intiated plot with i reffering to the location in the matrix 
for(i in 1:8){
  points(non_fpr_mat[i,], col= colorvec[i], lty=1)
  points(true_fpr_mat[i,], col=colorvec[i],lty=1)
}

legend("center", legend=c("AA","Apache","Bahamian","White","Chamorro","Filipino",
                           "Jamaican","Navajo","SEH","SWH","Trinidadian"),
       col=c( "blue","red","green","purple","violet","yellow","pink","black"), lty=1:2, cex=0.8)

?lty


################################### Plots for FPR counts ############################################################
######################################################################################################################

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
############################################## PLOT for  FPR  counts #####################################################################

#This is the vector used to make the different colored lines for all of the populations 
colorvec = c( "blue","red","green","purple","violet","yellow","pink","black")
symbol_vec = c('o',"+","*",".","-","s","n",",")
length(symbol_vec)

plot.new()

#This intiates the plot 
plot(x=c(),y=c(),xlim = c(1,8),ylim = c(0,1000),main =" False Positives  (1000 iterations)",ylab = "false postive", xlab = 
       "number of contributors")


#this loop adds points onto the intiated plot with i reffering to the location in the matrix 
for(i in 1:8){
  points(non_fprcount_mat[i,], col= colorvec[i], lty=1)
  #points(true_fprcount_mat[i,], col=colorvec[i],lty=1)
}

legend("center", legend=c("AA","Apache","Bahamian","White","Chamorro","Filipino",
                          "Jamaican","Navajo","SEH","SWH","Trinidadian"),
       col=c( "blue","red","green","purple","violet","yellow","pink","black"), lty=1:2, cex=0.8)

?lty
################### Rori's requested plot of FPR greater than -1 vs Expected heterozygozity of the 7th contrib

#This is the vector used to make the different colored lines for all of the populations 
colorvec = c("blue","red","green","purple","violet","yellow","pink","black","chartreuse",
             "coral","aquamarine4")
symbol_vec = c('o',"+","*",".","-","s","n",",")
length(symbol_vec)

plot.new()
?plot
#This intiates the plot 
plot(FPR7_matrix,main ="FPR vs Expected Heterozygosity (7 contributors)",ylab = "FPR", xlab = 
       "Expected Heterozygosity")
points(x=exphetero_mat, col= "blue", lty=5)
points(y=exphetero_mat)

#new matrix 
FPR7_matrix <-matrix(1:11,nrow = 11,ncol = 2)
FPR7_matrix[,2] <-non_fprcount_mat[,7]

col_vec <- c("Expected Heterozygosity","FPR")
colnames(FPR7_matrix) <- col_vec

rownames(FPR7_matrix) <- pop_vec2
FPR7_matrix[,1] <- exphetero_mat

FPR7_matrix

?color.scale


#This intiates the plot 
plot(x=c(),y=c(),xlim = c(1,15),ylim = c(.6,1),main =" FPR vs Expected Heterozygosity (7 contributors)",ylab = "Expected Heterozygosity", xlab = 
       "False Positives")
points(x=FPR7_matrix[1,],y=FPR7_matrix[,1],col="red")

#this loop adds points onto the intiated plot with i reffering to the location in the matrix 
for(i in ){
  points(FPR7_matrix[i,], col= colorvec[i], lty=1)
  #points(true_fprcount_mat[i,], col=colorvec[i],lty=1)
}

legend("center", legend=c("AA","Apache","Bahamian","White","Chamorro","Filipino",
                          "Jamaican","Navajo","SEH","SWH","Trinidadian"),
       col=c( "blue","red","green","purple","violet","yellow","pink","black"), lty=1:2, cex=0.8)


fpr_df <- as.data.frame(FPR7_matrix)

library(ggplot2)
g <- ggplot(fpr_df, aes(x =1:15, y = 0:1 , color = fp0))
g <- g + geom_point(shape = 15, aes(fill = fp0))

################### Rori's requested plot of FPR greater than -1 vs Expected heterozygozity of the 7th contrib
heterozygosity 
true_fprcount_mat[]



######################################### STEP:4 ###################################################################

# STEP 4: Creates plots of the difference between FPR 


fpr_diff_mat = matrix(1:88,nrow = 11,ncol = 8)

for (i in 1:11){
  population = pop_vec[i]
  
  for(j in 1:8){
    contrib = j
    count=1
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:100){
      sticks = array_3dnon[k,j,i]
      
      
      if(sticks > -1){
        
        count = count + 1
        print(count)
      }
      
      
    }
    
    
    result = (count - 1)/1000
    print(result)
    falsepositive_mat[i,j] = result
    
  }
}

fpr_diff_mat


########################################STEP:0  Friday 3/15/19 with Carina making plots#######################################
#STEP: 0 makes plots that we may not use but wull keep the code for
###this makes plots like the ones MAtt made that shows power and FPR getting closer to one another 
popNist_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
                 "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")

pdf("True_noncontrib_plots", width=20, height=20)

#This is how many plots per row,column 
par(mfrow=c(5,5))


# i 1:8 because There are 8 contributors 
for (i in 1:11){
  population = popNist_vec[i]
  #print(population)
  
  
  for(j in 1:8){
    contrib= j
    trueIteration_data = array_3d[, j, i]
    nonIteration_data = array_3dnon[, j,i]
    hist( trueIteration_data,col=rgb(1,0,0,0.5),prob = TRUE, ylim=c(0,1.0),
          main = paste(contrib, population), xlab="Likelihood Ratio (Log10)", ylab = "Frequency")
    hist(nonIteration_data,col=rgb(0,0,1,0.5), ylim=c(0,1.0), prob = TRUE, add=T)
    box()
    
    
  }
  
}

dev.off()



