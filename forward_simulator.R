## Zoya Wani zoyawani1@tamu.edu
## Andres Barboza P. andresdbp@tamu.edu
## April 16, 2024


###### Starting Conditions #########
N <- 100 #population
loci <- 10 #positions on the genome 
mu <- 10^-5 #human mutation rate 10^-9 for an individual nucleotide
baseval <- 10 # this is a base minimum value for our phenotype
loci.imp <- sort(sample(2:loci, loci/2))
opt <- 20
sigma <- 3
iter <- 1000
###### End Starting Conditions #########


###### FUNCTIONS #########



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
  
  mothers <- sample(1:(N/2), size = N, replace = TRUE, prob = w[1:(N/2)])
  fathers <- sample((N/2 + 1):N, size = N, replace = TRUE, prob = w[(N/2 + 1):N])
  
  return(list(mothers = mothers, fathers = fathers))
}

GetGametes <- function(mothers, fathers, pop) {
  genotype_lookup <- data.frame(
    genotype = c(1, 2, 3, 4),
    allele_1 = c(0, 0, 1, 1),
    allele_2 = c(0, 1, 0, 1)
  )
  
  get_maternal_alleles <- function(pop, mothers) {
    # Call the mothers and their genotypes
    genotype <- pop[mothers, ]
    
    
    matched_rows <- match(genotype, genotype_lookup$genotype)
    # Match rows of the genotype to the genotypes in the table 
    
    # Get alleles that match rows 
    allele_1 <- genotype_lookup$allele_1[matched_rows]
    allele_2 <- genotype_lookup$allele_2[matched_rows]
    
    return(list(allele_1 = allele_1, allele_2 = allele_2))
  }
  
  #same as mothers
  get_paternal_alleles <- function(pop, fathers) {
    # Get the genotypes of the fathers
    genotype <- pop[fathers, ]
    
    
    matched_rows <- match(genotype, genotype_lookup$genotype)
    
    # Extract alleles using the matched rows
    allele_1 <- genotype_lookup$allele_1[matched_rows]
    allele_2 <- genotype_lookup$allele_2[matched_rows]
    
    return(list(allele_1 = allele_1, allele_2 = allele_2))
  }
  
  # Get maternal alleles
  maternal_alleles <- get_maternal_alleles(pop, mothers)
  
  # Get paternal alleles
  paternal_alleles <- get_paternal_alleles(pop, fathers)
  
  return(list(
    maternal_alleles = maternal_alleles,
    paternal_alleles = paternal_alleles
  ))
}

MakeFertilization <- function(gametes){
  
  zygotes <- paste0(gametes$maternal_alleles$allele_1, "-", gametes$paternal_alleles$allele_1)
  
  lookup <- data.frame(
    zygote = c("0-0", "0-1", "1-0", "1-1"),
    genotype = c(1, 2, 3, 4)
  )
  
  genotype <- lookup$genotype[match(zygotes, lookup$zygote)]
  
  # Create new population matrix with updated genotypes
  new_pop <- matrix(genotype, nrow = length(zygotes), ncol = loci, byrow = TRUE)
  return(new_pop)
}

SimulateGenerations <- function(N, loci, mu, baseval, loci.imp, opt, iter, sigma) {
  pop <- GetPopulation(N, loci)
  avg_phenos <- numeric(iter)
  for (generation in 1:iter) {
    pop <- MutatePop(pop, mu)
    phenos <- GetPheno(pop, loci.imp, baseval)
    avg_phenos[generation] <- mean(phenos)
    w <- GetFit(phenos, opt, sigma)
    parents <- PickParents(pop, w)
    gametes <- GetGametes(parents$mothers, parents$fathers, pop)
    pop <- MakeFertilization(gametes)
    print(generation)
  }
  return(list(final_population = pop, avg_phenos = avg_phenos))
}

# Run the simulation
simulation_result <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt, iter, sigma)
final_population <- simulation_result$final_population
average_phenotypes <- simulation_result$avg_phenos

plot(average_phenotypes)

