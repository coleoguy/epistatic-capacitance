

N <- 10 
loci <- 20 
mu <- 10^-9
GetPopulation <- function(N,loci){
  

  # The numbers present in the matrix represent the genotype at a locus. 
  # First number is maternal and second is paternal 
  # 0,0 = 1
  # 0,1 = 2
  # 1,0 = 3
  # 1,1 = 4
  
  pop <- matrix(1, N, loci)
  pop[(N/2 + 1):N,1] <- 2
  # first column is sex determining locus with 0=X 1=Y 
  
  return(pop)
  
}


foop <- GetPopulation(N=20,loci=50)

MutatePop <- function(pop, mu){
  mutation.probability 1-(1-mu)^1547 #NOTHING MUTATED
  #select columns where the change of drawing a column is equal to the mut.prob
  #change genotype at selected location
  #repeat proccess for each individual
}
  
