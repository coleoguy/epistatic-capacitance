## Zoya Wani zoyawani1@tamu.edu
## Andres Barboza P. andresdbp@tamu.edu
## April 16, 2024


###### Starting Conditions #########
N <- 10 #population
loci <- 100 #positions on the genome 
mu <- 10^-5 #human mutation rate 10^-9 for an individual nucleotide
baseval <- 0 # this is a base minimum value for our phenotype
loci.imp <- sort(sample(2:loci, loci/10))
opt <- 100
sigma <- 10
gen <- 50
arch <- "add" # add, addxadd
###### End Starting Conditions #########

###### FUNCTIONS #########

GetPopulation <- function(N,loci){
  
  # The numbers present in the matrix represent the genotype at a locus. 
  # First number is maternal and second is paternal 
  # 0,0 = 1
  # 0,1 = 2
  # 1,0 = 3
  # 1,1 = 4
  pop <- matrix(sample(1:4, N*loci, replace=T), N, loci)
  
  # We are moving away into a model without sex, so the bottom part is commented out
  
  # set as female
  # pop[1:(N/2),1] <- 1
  # set half as male
  # pop[(N/2 + 1):N,1] <- 2
  # first column is sex determining locus with 1=XX 2=XY 
  # remaining columns are the rest of the genome
  # each row in pop is one individual
  
  
  return(pop)
  
}

MutatePop <- function(pop, mu) {
  # Probability of an individual having at least one mutation
  mut.prob <- 1 - (1 - mu)^1547  # 1547 is the number of base pairs in the coding region
  
  # Determine mutation coordinates where random values are less than mut.prob
  mut.coord <- which(matrix(runif(nrow(pop) * ncol(pop)), 
                            nrow = nrow(pop), ncol = ncol(pop)) < mut.prob, 
                     arr.ind = TRUE)
  
  apply_mutations <- function(pop, coordinates) {
    for (i in seq_len(nrow(coordinates))) {
      row_index <- coordinates[i, 1]
      col_index <- coordinates[i, 2]
      value <- pop[row_index, col_index]
      
      # switches values 
      pop[row_index, col_index] <- switch(value,
                                          sample(c(2, 3), 1),  # for 1
                                          sample(c(1, 4), 1),  # for 2
                                          sample(c(1, 4), 1),  # for 3
                                          sample(c(2, 3), 1))  # for 4
    }
    return(pop)
  }
  
  # Apply mutations to the population matrix
  mutants <- apply_mutations(pop, mut.coord)
  return(mutants)
}

GetPheno <- function(pop, loci.imp, baseval, arch){
  temppop <- pop[,loci.imp]
  if(arch == "add") {
    temppop[temppop==1] <- 0
    temppop[temppop %in% c(2, 3)] <- 1
    temppop[temppop == 4] <- 2
  }
  if(arch == "addxadd") {
    # Epistasis: additive by additive
  }
  phenos <- rowSums(temppop) + baseval
  return(phenos)
}

GetFit <- function(obs, opt, sigma){
  numer <- (obs - opt)^2
  denom <- (2 * sigma)^2
  w <- exp(-(numer / denom))
  return(w)
}

Reproduction <- function(pop, N, w, loci) {
  # Sample parents based on fitness probabilities
  parents <- sample(N, size = 2 * N, replace = TRUE, prob = w)
  
  # Generate haplotype selection matrix
  haplotypes <- matrix(sample(1:2, 2 * N * loci, replace = TRUE), nrow = 2 * N, ncol = loci)
  
  # Gamete formation
  gametes <- ifelse(
    pop[parents, ] == 1,  # Genotype 1 -> Both alleles are 0
    0,
    ifelse(
      pop[parents, ] == 4,  # Genotype 4 -> Both alleles are 1
      1,
      ifelse(
        pop[parents, ] == 2,  # Genotype 2 -> Depends on haplotype
        ifelse(haplotypes == 1, 0, 1),
        ifelse(haplotypes == 1, 1, 0)  # Genotype 3 -> Depends on haplotype
      )
    )
  )
  
  # Fertilization: Pair gametes to form diploid genotypes
  maternal_gametes <- gametes[seq(1, 2 * N, by = 2), ]
  paternal_gametes <- gametes[seq(2, 2 * N, by = 2), ]
  new_population <- 1 + maternal_gametes + paternal_gametes * 2
  
  return(new_population)
}

SimulateGenerations <- function(N, loci, mu, baseval, loci.imp, opt, gen, sigma) {
  pop <- GetPopulation(N, loci)
  avg_phenos <- numeric(gen)
  for (generation in 1:gen) {
    pop <- MutatePop(pop, mu)
    phenos <- GetPheno(pop, loci.imp, baseval)
    avg_phenos[generation] <- mean(phenos)
    w <- GetFit(phenos, opt, sigma)
    pop <- Reproduction(pop, N, w, loci)
    print(generation)
  }
  return(list(final_population = pop, avg_phenos = avg_phenos))
}

# Run the simulation
iter <- 50
for (i in 1:iter) {
  simulation_result <- SimulateGenerations(N = 1000, loci, mu, baseval, loci.imp, opt, gen, sigma)
  lines(simulation_result$avg_phenos, col = "#314CB6")
  if(i == 1) {
    plot(simulation_result$avg_phenos, type = "l", col = "#314CB6",
         ylim = c(min(simulation_result$avg_phenos), max(simulation_result$avg_phenos)))
  }
}
