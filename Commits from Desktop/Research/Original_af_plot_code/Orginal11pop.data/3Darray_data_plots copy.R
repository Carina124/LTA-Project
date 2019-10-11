######################################## Friday 3/15/19 with Carina making plots#######################################

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

######################### 3/20/19 Maria,Carina, Plot for Rori#############################

count = 1
for (i in 1:11){
  population = popNist_vec[i]

  for(j in 1:8){
  contrib = j
  #sticks is iterations 
  sticks = array_3d[,j,i]
  
  if(sticks >= 1){
    count = count + 1 
    
  }
  

}}

######################
count = 1

issues =  c("NaN","NA","Inf")


for (i in 1:2){
  population = popNist_vec[i]
  
  for(j in 1:7){
    contrib = j
    #sticks is iterations 
    #sticks = array_3d[,j,i]
    for(k in 1:100){
      sticks = array_3dnon[k,j,i]
     

       if(sticks >= 1){
         
        count = count + 1 
        #print(count)
       }
        #ifelse(sticks %in% issues,k,k)
        
    }
    
  }}
testmat = matrix(data=count)

############################




########################################Final code that counts FPR carina 3/20/19############################################################

testmat = matrix(1:88,nrow = 11,ncol = 8)

for (i in 1:11){
  population = popNist_vec[i]
  
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
     
    
    result = (count - 1)/100
    print(result)
    testmat[i,j] = result
  
    }
}

testmat
############################### Plot for False Postive Rates after producing the Matrix above ######################################

#The files that will be used 
popNist_vec <- c("afri_paper.csv","apache.csv","Bahamian.csv","white.csv","Chamorro.csv","Filipino.csv",
                 "Jamaican.csv","Navajo.csv","seh.csv","swh.csv","Trinidadian.csv")

#The names of populations to add to the rows of the matrix 
popNist_vec2 <- c("afri","apache","Bahamian","white","Chamorro","Filipino",
                  "Jamaican","Navajo","seh","swh","Trinidadian")

#Starts a PDF so the plots can be stored 
pdf("FPR_noncontrib_plot", width=20, height=20)

#This is how many plots per row,column 
par(mfrow=c(1,1))



#This is all naming schemes for the matrix row and columns 
iterations = 1:8
colnames(testmat) <- iterations
rownames(testmat) <- popNist_vec2
testmat




#This is the vector used to make the different colored lines for all of the populations 
colorvec = c( "blue","red","green","purple","violet","yellow","pink","black")
symbol_vec = c('o',"+","*","^","#","$","-","a")
length(symbol_vec)

#This intiates the plot 
plot(x=c(),y=c(),xlim = c(1,8),ylim = c(1,88))


#this loop adds points onto the intiated plot with i reffering to the location in the matrix 
for(i in 1:8){
  points(testmat[i,], type="o", pch=symbol_vec[i], lty=1)
}


##################################ggplots alone Monday 3/18/19#############################################################
install.packages("ggplot2")
library(ggplot2)

#cannot use gg plot on an array? 
array_3d_vec = c(array_3d[,1,1])
array_3d_vec
typeof(array_3d_vec)

#ggplot code 
ggplot(data=array_3d_vec)
ggplot(data=array_3d_vec,aes(x="contrib"))+geom_histogram()


?lims

















