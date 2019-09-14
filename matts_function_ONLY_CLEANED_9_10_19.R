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

################################
## Ignore, For Debugging Only ##
################################

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
  
  #######simugeno onjects store genotypes from the tabfreq ###########
  #popgen$tab.geno gives the genotyprs of all individuals (n) 
  #look at a particular indivduals genotype: popgeno$tab.geno[1:10,]
  
  
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
 
  
  noncontrib <- simugeno(tab = pop.afs, 
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
  
  #SIMUMIX objects store DNA mixtures#########
  ####noncontrib is the simulated genotype data
  #### ncontri is the number of contributors to the mixture 
  noncon.sus <- simumix(noncontrib, 
                        ncontri = 3
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
  
  singleLR_vector <- c()
  total_LR <- c()
  LR_vector <- c()
  log10_LR <- c()
  
  
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
          
          print(single_LR$LR)
          print(LR_vector <- c(LR_vector,single_LR$LR))
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
          
         print(single_LR$LR)
         print(single_LR)
          
          while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0){
            single_LR <- LR( Repliste = c(sim.mix@mix.all[[k]]),
                             Tp = c(prosecutor.all.atk),
                             Td = c(defense.all.atk), 
                             ## Vd = Non contrib under Hd - Suspect ##
                             Vp = 0,
                             ## Vp = Non contrib under Hp - 0 ##
                             Vd = c(as.numeric(tail( prosecutor.all.atk,2), "/")[[1]]),
                             xd = 1,
                             xp = 0,
                             theta = 0,
                             prDHet = c(0.2,0.2),
                             prDHom = c(0.04,0.04),
                             prC = 0,
                             freq = pop.afs.matrix[[k,1]]
                             )
                             
            #singleLR_vector starts collecting LRs for all 13 loci
            singleLR_vector <- c(LR_vector,single_LR$LR)
          
          #end of while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0)  
          }
         
           singleLR_vector <<- c(singleLR_vector,single_LR$LR)
           k = k + 1
           
       #end of while (k < 14), for LR calculation   
       }
     
     #end of else, for more than 1 contributor   
     }
     
      #when K == 14 the entire vector is multiplied here
      log10_LR[j] <<- log10(prod(singleLR_vector))
      #j moves on to the next simulation until there are none left to do
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
      
      singleLR_vector <- c()  
      log10_LRvector <- c()
      
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
          c = num_contribs - 1
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
            
            for (i in 1:c){
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
            
              while(is.na(single_LR$LR) || is.infinite(single_LR$LR) || is.nan(single_LR$LR) || single_LR$LR<0){
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
              }
                singleLR_vector <- c(singleLR_vector,single_LR$LR)
                k = k + 1
          
        #end of  while (k < 14)  
        }
        
      #end of else, for more than 1 contributor  
      }
      
        log10_LR[j] <- log10(prod(singleLR_vector))
        
        #log10_totalLRs_vec <- c(log10(total_LR))
        log10_LR[log10_LR < 0.0] <- -1
        
        #total_LR_plot_vector_nonnav8 <<- log10_totalLRs_vec
        j = j + 1
    
    #end of while (j < num_sims), which loops over simulations
    }
    
  #end of else if (non.or.truecontrib == 1)  
  }
    
    #singleLRs_vec_product[j] - this is what I want reterned 
    return(log10_LR) 
  
#end of total_contrib_function <- function(file.name,pop.size,num_contribs,num_sims,non.or.truecontrib)  
}


#########################
#### END OF FUNCTION ####
#########################

#file.name, pop.size, num_contribs, num_sims, non.or.truecontrib

LRs = total_contrib_function("USAf.1_new.csv", 1000, 3, 100, 0)

