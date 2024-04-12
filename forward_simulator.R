# heading



#package to look at 
https://search.r-project.org/CRAN/refmans/simcross/html/sim_meiosis.html



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


# A: Let's keep pick parents as a random sampling and not worry about fitness yet
PickParents <- function(pop){

  parents <- sample(1:nrow(pop), size = N, replace = TRUE, prob = )
  return (parents)



  

gametes 
}

#### TODO #A:
## Finish MutatePop
# Let's get it to a point were we can start a population and change
# some of their genomes
#
## Start GetGametes
# This steps comes after we have picked parents to mate, since we only make
# one "successfull" gamete for each time a parent mates.
# There's a some moving pieces that need to work in this step:
# 1) We need to select at which point recombination will happen, preferably
#    having a different random recombination spot for each parent
# 2) We need to select one strand to start making a gamete and a method to 
#    switch strand after the recombination spot
# 3) We need to actually create the gametes, which should contain 0 and 1 as
#    as their genotype, with values taken from the parent genotype dependent
#    on the strand and recombination spot
#
## Fertilization
# Since we should have twice as many gametes, we need to pick two random gametes
# and have them go trough a compressing step that will take the gametes from 0,1
# genotypes with two row to a 1-4 genotype on a single row. This step results
# in the next generation's initial population
#
#
#### This are the steps required to get a population that can complete
#### it's life cycle and change will happen due to drift and mutations.
#### The steps bellow will allow us to apply selection to our population
#### But if you finish all the steps above we should meet before we start
#### working on them.
#
#
## GetPhenotype Function
# This requires mapping the genotype to some phenotype (let's not worry about this one yet)
#
## GetFitness Function
# This requires a phenotype that will be compared to an optimum. Should also record population's mean fitness
#
# After this we will add the steps specific to our study (measuring epi and add
# variance, making an optimum shift, creating a bottleneck, etc.)

