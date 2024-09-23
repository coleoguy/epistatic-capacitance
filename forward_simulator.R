## Zoya Wani zoyawani1@tamu.edu
## Andres Barboza P. andresdbp@tamu.edu
## April 16, 2024

###### Starting Conditions #########
N <- 1000 #population
loci <- 100 #positions on the genome 
mu <- 10^-5 #human mutation rate 10^-9 for an individual nucleotide
baseval <- 10 # this is a base minimum value for our phenotype
loci.imp <- sort(sample(2:loci, 10))
opt <- 15
###### End Starting Conditions #########


###### FUNCTIONS #########

# Helper Functions #
# TODO Try commenting out this function and then go to line 94 and change
# using the apply_mutation function instead
apply_switch <- function(pop, coordinates) {
  # Create a function to apply the switch statement to a single element
  apply_switch_single <- function(row_index, col_index) {
    value <- pop[row_index, col_index]
    pop[row_index, col_index] <- switch(value,
                                        sample(c(2,3), 1),  #for 1
                                        sample(c(1,4), 1),  #for 2
                                        sample(c(1,4), 1),  #for 3
                                        sample(c(2,3), 1))  #for 4
    return(pop[row_index, col_index])  # return the updated value
  }
  # Apply the function to each row of coordinates and collect results
  updated_values <- mapply(apply_switch_single, coordinates[,1], coordinates[,2])
  # Update the pop matrix with the updated values
  pop[coordinates] <- updated_values
  return(pop)
}

#This function is how we are carrying out mutations and updating them in the sample 
# TODO add helper functin that converts diploid vector to haplotype
# End of Helper Functions #


GetPopulation <- function(N,loci){
  
  # The numbers present in the matrix represent the genotype at a locus. 
  # First number is maternal and second is paternal 
  # 0,0 = 1
  # 0,1 = 2
  # 1,0 = 3
  # 1,1 = 4
  tot.loci <- N*loci
  pop <- matrix(sample(1:4, tot.loci, replace=T), N, loci)
  #set as female
  pop[1:(N/2),1] <- 1
  # set half as male
  pop[(N/2 + 1):N,1] <- 2
  # first column is sex determining locus with 1=XX 2=XY 
  # remaining columns are the rest of the genome
  # each row in pop is one individual
  return(pop)
  
}

# TODO We need to determine if this is giving us what we want with regard to the 
# number of mutations being entered into the population.
MutatePop <- function(pop, mu){
  # TODO add a comment that explains what this is the probability of a 
  # mutation at any given locus the probability of an individual having one 
  # mutation the average number of mutations per an individual?
  mut.prob <- 1-(1-mu)^1547 # 1547 is the number of base pairs in the coding region of our genome (treating it as the length of a gene)
  # TODO add a comment that explains what this does.
  mut.coord <- which(matrix(runif(nrow(pop) * ncol(pop)), nrow = nrow(pop), ncol = ncol(pop)) < mut.prob, arr.ind = TRUE)
  
  # This switch function applies the mutations as they are produced
  apply_mutations <- function(pop, coordinates) {
    # Create a function to apply the switch statement to a single element
    apply_switch_single <- function(row_index, col_index) {
      value <- pop[row_index, col_index]
      pop[row_index, col_index] <- switch(value,
                                          sample(c(2,3), 1),  #for 1
                                          sample(c(1,4), 1),  #for 2
                                          sample(c(1,4), 1),  #for 3
                                          sample(c(2,3), 1))  #for 4
      return(pop[row_index, col_index])  # return the updated value
    }
    # Apply the switch function to each row of coordinates and collect results
    updated_values <- mapply(apply_switch_single, coordinates[,1], coordinates[,2])
    # Update the pop matrix with the updated values
    pop[coordinates] <- updated_values
    
    return(pop)
  }

  mutants <- apply_switch(pop, mut.coord)
  return(mutants)
}

GetPheno <- function(pop, loci.imp, baseval){
  temppop <- pop[,loci.imp]
  temppop[temppop==1] <- 0
  temppop[temppop %in% c(2, 3)] <- 1
  temppop[temppop == 4] <- 2
  phenos <- rowSums(temppop) + baseval
  return(phenos)
}

GetFit <- function(obs, opt, sigma){
  numer <- (obs - opt)^2
  denom <- (2 * sigma)^2
  w <- exp(-(numer / denom))
  return(w)
  }

PickParents <- function(pop, w){

}


# 1) function simplification repettitve functions
# 2) picking parents
 
###### END FUNCTIONS #########


###### Running Sims ##########
pop <- GetPopulation(N = N, loci = loci)
pop <- MutatePop(pop = pop, mu = mu)
phenos <- GetPheno(pop = pop, loci.imp = loci.imp, baseval = baseval)
w <- GetFit(obs=phenos, opt=opt, sigma=2)

# TODO put this into a function up above and combine the two vectors
# perhaps into a list of length 2
mothers <- sample(1:(N/2), size = N, replace = T, prob = w[1:(N/2)])
fathers <- sample((N/2 +1):N, size = N, replace = T, prob = w[(N/2 +1):N])

# TODO start wokring on the get gametes function

###### Running Sims ##########




#### TODO #A:
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

