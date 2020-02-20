# READ IN NIQUOS AVERAGE EXP HET FILE 
NiqAveExpHet = read.csv('Average_Expected_Hetero_9.19.csv', header = T)

#READ IN MY AVERAGE EXP HET FILE 
MariaAveExpHet = read.csv('430PopAverageExpHetero.csv', header = T)

# MERGE BOTH FILES TO BE ABLE TO DO COMPARISON OF AVE EXP HET VALUES 
MergedAveExpHET = merge(NiqAveExpHet, MariaAveExpHet, by=c('X'), all.x = T)

# DIFFERENCE IS CALCULATED BETWEEN MY AVE EXP HET VALUES AND NIQUOS AVE EXP HET VALUES
MergedAveExpHET$Difference = (MergedAveExpHET$Average.Expected.Heterozygosity - MergedAveExpHET$AverageExpHet)
MergedAveExpHET
