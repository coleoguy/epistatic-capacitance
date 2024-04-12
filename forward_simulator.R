# heading

N <- 10 #population
loci <- 100 #positions on the genome 
mu <- 10^-5 #human mutation rate 10^-9 for an individual nucleotide

GetPopulation <- function(N,loci){
  
  # The numbers present in the matrix represent the genotype at a locus. 
  # First number is maternal and second is paternal 
  # 0,0 = 1
  # 0,1 = 2
  # 1,0 = 3
  # 1,1 = 4
  
  pop <- matrix(1, N, loci)
  pop[(N/2 + 1):N,1] <- 2
  # first column is sex determining locus with 1=X 2=Y 
  
  return(pop)
  
}


pop <- GetPopulation(N=N,loci=loci)

MutatePop <- function(pop, mu){
  mut.prob <- 1-(1-mu)^1547 # 1547 is the number of base pairs in the coding region of our genome (treating it as the length of a gene)
  mutants <- rbinom(n=nrow(pop), size=ncol(pop), mut.prob)
  for(i in 1:nrow(pop)){
    if(mutants[i] > 0){
      hit <- sample(2:ncol(pop), size=mutants[i], replace=F)
      for(j in 1:length(hit)){
        switch(hit[j],
               sample(c(2,3), 1),  #for 1
               sample(c(1,4), 1),  #for 2
               sample(c(1,4), 1),  #for 3
               sample(c(2,3), 1))  #for 4
      
              
        return(hit[j])
               
        
        #select columns where the change of drawing a column is equal to the mut.prob
        #change genotype at selected location
        #repeat proccess for each individual               
               
      }
    }
  }
  

}
  
PickParents <- function(pop){

  parents <- sample(1:nrow(pop), size = N, replace = TRUE, prob = )
  return (parents)



  

gametes 
}

