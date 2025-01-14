## Zoya Wani zoyawani1@tamu.edu
## Andres Barboza P. andresdbp@tamu.edu
## April 16, 2024


###### Starting Conditions #########
N <- 500 #population
loci <- 100 # positions on the genome 
mu <- 10^-5 # human mutation rate 10^-9 for an individual nucleotide
baseval <- 0 # this is a base minimum value for our phenotype
loci.imp <- sort(sample(1:loci, loci/10))
opt <- 10
sigma <- 10
gen <- 100
arch <- "sign" # add, sign, inc, dec
sag <- 1.1
sign_flag <- "alter" # half, alter
###### End Starting Conditions #########

###### FUNCTIONS #########

GetPopulation <- function(N,loci){

  #pop <- matrix(sample(1:4, N*loci, replace=T), N, loci)
  pop <- matrix(rep(4, N*loci), N, loci) 

  # 0,0 = 1
  # 0,1 = 2
  # 1,0 = 3
  # 1,1 = 4
  
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

GetPheno <- function(pop, loci.imp, baseval, arch, sign_flag){
  temppop <- pop[,loci.imp]
  temppop[temppop==1] <- 0
  temppop[temppop %in% c(2, 3)] <- 1
  temppop[temppop == 4] <- 2
  
  if(arch == "add") {
    phenos <- rowSums(temppop) + baseval
  }
  
  if(arch == "sign") {
    if (sign_flag == "half") {
      dir_loci <- temppop[, 1:(length(loci.imp)/2)]
      inv_loci <- temppop[, (length(loci.imp)/2)+1:(length(loci.imp)/2)]
    } else if (sign_flag == "alter") {
      dir_loci <- temppop[, seq(1, length(loci.imp), by = 2)]
      inv_loci <- temppop[, seq(2, length(loci.imp), by = 2)]
    }
    
    cond1 <- (dir_loci == 2 & inv_loci == 2) | (dir_loci == 0 & inv_loci == 0)
    cond2 <- (dir_loci == 2 & inv_loci == 0) | (dir_loci == 0 & inv_loci == 2)
    cond3 <- (dir_loci == 1 | inv_loci == 1)

    phenos <- baseval + rowSums(4 * cond1 + 0 * cond2 + 2 * cond3)
  }
  
  if(arch == "inc") {
    phenos <- (length(loci.imp)*2)*((rowSums(temppop) + baseval)/(length(loci.imp)*2))^sag
  }
  
  if(arch == "dec") {
    phenos <- baseval + (length(loci.imp) * 2 - baseval) * ((rowSums(temppop) + baseval - baseval) / (length(loci.imp) * 2 - baseval))^(1/sag)
  }
  
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
  recomb <- pmax(pmin(rpois(n = 2 * N, lambda = 10), loci), 1)
  haplotypes <- matrix(unlist(lapply(recomb, function(r) rep(sample(1:2, r, replace = TRUE), each = ceiling(loci / r), length.out = loci))),
                       nrow = 2 * N, ncol = loci, byrow = TRUE)

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

GetArch <- function(pop, loci.imp, phenos) {
  
  pop <- pop[, loci.imp]
  
  # Additive
  a_mat <- matrix(0, nrow = N, ncol = length(loci.imp))
  a_mat[pop == 1] <- 1
  a_mat[pop == 4] <- -1
  
  fit <- lm(phenos ~ a_mat)
  R2a <- summary(fit)$adj.r.squared
  
  # Dominance
  d_mat <- matrix(1, nrow = N, ncol = length(loci.imp))
  d_mat[pop == 1 | pop == 4] <- 0
  
  fit <- lm(phenos ~ a_mat + d_mat)
  R2ad <- summary(fit)$adj.r.squared
  
  # Epistasis
  # Precompute pairs of indices for interaction terms
  if (sign_flag == "alter") {
    indices <- seq(1, ncol(pop), by = 2)
  } else {
    half_cols <- ncol(pop) / 2
    indices <- 1:half_cols
  }
  
  # Compute epistatic interaction matrices
  e_mat_aa <- a_mat[, indices] * a_mat[, indices + ifelse(sign_flag == "alter", 1, half_cols)]
  e_mat_ad <- a_mat[, indices] * d_mat[, indices + ifelse(sign_flag == "alter", 1, half_cols)]
  e_mat_dd <- d_mat[, indices] * d_mat[, indices + ifelse(sign_flag == "alter", 1, half_cols)]
  
  fit <- lm(phenos ~ a_mat + d_mat + e_mat_aa + e_mat_ad + e_mat_dd)
  R2ade <- summary(fit)$adj.r.squared
  
  return(c(R2a, R2ad, R2ade))
}

SimulateGenerations <- function(N, loci, mu, baseval, loci.imp, opt, gen, sigma, arch, sign_flag, verbose) {
  pop <- GetPopulation(N, loci)
  avg_phenos <- numeric(gen)
  lm_arch <- matrix(NA, nrow = gen, ncol = 3)
  for (generation in 1:gen) {
    pop <- MutatePop(pop, mu)
    phenos <- GetPheno(pop, loci.imp, baseval, arch, sign_flag)
    avg_phenos[generation] <- mean(phenos)
    lm_arch[generation, ] <- GetArch(pop, loci.imp, phenos)
    w <- GetFit(phenos, opt, sigma)
    pop <- Reproduction(pop, N, w, loci)
    if(verbose) print(generation)
  }
  return(list(final_population = pop, avg_phenos = avg_phenos, lm_arch = lm_arch))
}

# Run the simulation
iter <- 5
pheno <- add <- dom <- epi <- matrix(NA, nrow = iter, ncol = gen)
for (i in 1:iter) {
  print(i)
  simulation_result <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt, gen, sigma, arch, sign_flag, verbose=F)
  add[i,] <- simulation_result$lm_arch[,1]
  dom[i,] <- (simulation_result$lm_arch[,2]-simulation_result$lm_arch[,1])
  epi[i,] <- (simulation_result$lm_arch[,3]-simulation_result$lm_arch[,2])
  pheno[i,] <- simulation_result$avg_phenos
}

plot(colMeans(epi), type = "l", col = rgb(0.06,0.24,0.49,0.8), ylim = c(0, 1), lwd = 3)
lines(colMeans(add), col = rgb(0.89,0.34,0.18,0.8), lwd = 3)
lines(colMeans(pheno)/20, col = rgb(0.45,0.66,0.46,0.8), lwd = 3)



# Compare selection to epistasis eq
num_sigma <- 100
iter <- 100
sigma_values <- seq(1, 10, length.out = num_sigma)
epi_final_values <- numeric(num_sigma)

for (s in seq_along(sigma_values)) {
  
  sigma_now <- sigma_values[s]
  final_epi_reps <- numeric(iter)
  
  for (i in 1:iter) {
    sim_result <- SimulateGenerations(N, 
                                      loci, 
                                      mu, 
                                      baseval, 
                                      loci.imp, 
                                      opt, 
                                      gen, 
                                      sigma_now, 
                                      arch, 
                                      sign_flag, 
                                      verbose = FALSE)
    epi_values <- sim_result$lm_arch[,3] - sim_result$lm_arch[,2]
    final_epi_reps[i] <- epi_values[gen]
  }
  epi_final_values[s] <- mean(final_epi_reps)
}

cor_value <- cor(sigma_values, epi_final_values)
plot(sigma_values, epi_final_values, 
     type = "l", 
     xlab = expression(sigma), 
     ylab = "Equilibrium Epistasis",
     main = "Equilibrium Epistasis vs. Sigma")




# Compare selection to epistasis eq
library(viridis)
iter <- 100
sigma_values <- seq(1, 10, length.out = 10)
cols <- viridis(length(sigma_values))
plot(1:gen, 
     rep(NA, gen), 
     type = "n", 
     xlab = "Generation", 
     ylab = "Mean Epistasis (colMeans of 100 iter)", 
     ylim = c(0, 1))

for (s in seq_along(sigma_values)) {
  print(s)
  sigma_now <- sigma_values[s]
  epi_data <- matrix(NA, nrow = iter, ncol = gen)
  for (i in 1:iter) {
    sim_result <- SimulateGenerations(N, 
                                      loci, 
                                      mu, 
                                      baseval, 
                                      loci.imp, 
                                      opt, 
                                      gen, 
                                      sigma_now, 
                                      arch, 
                                      sign_flag, 
                                      verbose = FALSE)
    epi_data[i, ] <- sim_result$lm_arch[,3] - sim_result$lm_arch[,2]
  }
  lines(1:gen, colMeans(epi_data), 
        col = cols[s], 
        lwd = 2)
}
legend("bottomright", 
       legend = paste("Sigma =", round(sigma_values, 2)), 
       col = cols, 
       lty = 1, 
       lwd = 2)


# snippet
# plot(GetFit(obs=seq(from=0,to=20, length.out=100), 10, sigma=1)~seq(from=0,to=20, length.out=100))
