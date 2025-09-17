## Zoya Wani zoyawani1@tamu.edu
## Andres Barboza P. andresdbp@tamu.edu
## April 16, 2024

library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)

ncores <- 4

###### Starting Conditions #########
N <- 1000 #population
loci <- 100 # positions on the genome 
mu <- 10^-5 # human mutation rate 10^-9 for an individual nucleotide
baseval <- 0 # this is a base minimum value for our phenotype
num.imp.loci <- 10
loci.imp <- sort(sample(1:loci, num.imp.loci))
opt <- num.imp.loci
sigma <- num.imp.loci/2
gen <- 2000
arch <- "axd" # add, axa, axd, dxa, dxd, inc, dec
sag <- 1.1
epi_flag <- "alter" # half, alter
###### End Starting Conditions #########


###### FUNCTIONS #########

GetPopulation <- function(N,loci){

  #pop <- matrix(sample(1:4, N*loci, replace=T), N, loci)
  pop <- matrix(rep(4, N*loci), N, loci)
  #pop <- matrix(rep(4, N*loci), N, loci)

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

GetPheno <- function(pop, loci.imp, baseval, arch, epi_flag){
  temppop <- pop[,loci.imp]
  temppop[temppop==1] <- 0
  temppop[temppop %in% c(2, 3)] <- 1
  temppop[temppop == 4] <- 2
  
  if(arch == "add") {
    phenos <- rowSums(temppop) + baseval
  }
  
  if(arch == "axa") {
    if (epi_flag == "half") {
      dir_loci <- temppop[, 1:(length(loci.imp)/2)]
      inv_loci <- temppop[, (length(loci.imp)/2)+1:(length(loci.imp)/2)]
    } else if (epi_flag == "alter") {
      dir_loci <- temppop[, seq(1, length(loci.imp), by = 2)]
      inv_loci <- temppop[, seq(2, length(loci.imp), by = 2)]
    }
    
    cond1 <- (dir_loci == 2 & inv_loci == 2) | (dir_loci == 0 & inv_loci == 0)
    cond2 <- (dir_loci == 2 & inv_loci == 0) | (dir_loci == 0 & inv_loci == 2)
    cond3 <- (dir_loci == 1 | inv_loci == 1)

    phenos <- baseval + rowSums(4 * cond1 + 0 * cond2 + 2 * cond3)
  }
  
  if(arch == "dxd") {
    if (epi_flag == "half") {
      dir_loci <- temppop[, 1:(length(loci.imp)/2)]
      inv_loci <- temppop[, (length(loci.imp)/2)+1:(length(loci.imp)/2)]
    } else if (epi_flag == "alter") {
      dir_loci <- temppop[, seq(1, length(loci.imp), by = 2)]
      inv_loci <- temppop[, seq(2, length(loci.imp), by = 2)]
    }
    
    cond1 <- (dir_loci == 2 & inv_loci != 1) | (dir_loci == 0 & inv_loci != 1) | (dir_loci != 1 & inv_loci == 2) | (dir_loci != 1 & inv_loci == 0)
    cond2 <- (dir_loci == 1 & inv_loci != 1) | (dir_loci != 1 & inv_loci == 1)
    cond3 <- (dir_loci == 1 & inv_loci == 1)
    
    phenos <- baseval + rowSums(2 * cond1 + 0 * cond2 + 4 * cond3)
  }
  
  if(arch == "axd") {
    if (epi_flag == "half") {
      dir_loci <- temppop[, 1:(length(loci.imp)/2)]
      inv_loci <- temppop[, (length(loci.imp)/2)+1:(length(loci.imp)/2)]
    } else if (epi_flag == "alter") {
      dir_loci <- temppop[, seq(1, length(loci.imp), by = 2)]
      inv_loci <- temppop[, seq(2, length(loci.imp), by = 2)]
    }
    
    cond1 <- (dir_loci == 2 & inv_loci != 1)
    cond2 <- (dir_loci == 1)
    cond3 <- (dir_loci == 0 & inv_loci != 1)
    cond4 <- (dir_loci == 2 & inv_loci == 1)
    cond5 <- (dir_loci == 0 & inv_loci != 1)
      
    phenos <- baseval + rowSums(1 * cond1 + 2 * cond2 + 3 * cond3 + 4 * cond4 + 0 * cond5)
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
  if (epi_flag == "alter") {
    indices <- seq(1, ncol(pop), by = 2)
  } else {
    half_cols <- ncol(pop) / 2
    indices <- 1:half_cols
  }
  
  # Compute epistatic interaction matrices
  e_mat_aa <- a_mat[, indices] * a_mat[, indices + ifelse(epi_flag == "alter", 1, half_cols)]
  e_mat_ad <- a_mat[, indices] * d_mat[, indices + ifelse(epi_flag == "alter", 1, half_cols)]
  e_mat_dd <- d_mat[, indices] * d_mat[, indices + ifelse(epi_flag == "alter", 1, half_cols)]
  
  fit <- lm(phenos ~ a_mat + d_mat + e_mat_aa + e_mat_ad + e_mat_dd)
  R2ade <- summary(fit)$adj.r.squared
  
  return(c(R2a, R2ad, R2ade))
}

SimulateGenerations <- function(N, loci, mu, baseval, loci.imp, opt, gen, sigma, arch, epi_flag, verbose) {
  pop <- GetPopulation(N, loci)
  avg_phenos <- numeric(gen)
  lm_arch <- matrix(NA, nrow = gen, ncol = 3)
  for (generation in 1:gen) {
    pop <- MutatePop(pop, mu)
    phenos <- GetPheno(pop, loci.imp, baseval, arch, epi_flag)
    avg_phenos[generation] <- mean(phenos)
    lm_arch[generation, ] <- GetArch(pop, loci.imp, phenos)
    w <- GetFit(phenos, opt, sigma)
    pop <- Reproduction(pop, N, w, loci)
    if(verbose) print(generation)
  }
  return(list(final_population = pop, avg_phenos = avg_phenos, lm_arch = lm_arch))
}

###### ANALYSIS #########

#### NUMBER OF LOCI AND Ve ####

iter <- 200
loci_range <- seq(4, 100, by = 2)

arch <- "axa"
results_axa <- mclapply(loci_range, function(l) {
  final_epi_reps <- numeric(iter)
  for(i in 1:iter){
    loci.imp <- sort(sample(1:loci, l))
    opt <- l
    sigma <- l/2
    sim_result <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt,
                                      gen, sigma, arch, epi_flag, verbose = FALSE)
    epi_values <- sim_result$lm_arch[,3] - sim_result$lm_arch[,2]
    final_epi_reps[i] <- epi_values[gen]
  }
  mean_epi <- mean(final_epi_reps)
  sd_epi <- sd(final_epi_reps)
  data.frame(loci = l, mean_epi = mean_epi, sd_epi = sd_epi)
}, mc.cores = ncores)
df_axa <- do.call(rbind, results_axa) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1)
  )

arch <- "axd"
results_axd <- mclapply(loci_range, function(l) {
  final_epi_reps <- numeric(iter)
  for(i in 1:iter){
    loci.imp <- sort(sample(1:loci, l))
    opt <- l
    sigma <- l/2
    sim_result <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt,
                                      gen, sigma, arch, epi_flag, verbose = FALSE)
    epi_values <- sim_result$lm_arch[,3] - sim_result$lm_arch[,2]
    final_epi_reps[i] <- epi_values[gen]
  }
  mean_epi <- mean(final_epi_reps)
  sd_epi <- sd(final_epi_reps)
  data.frame(loci = l, mean_epi = mean_epi, sd_epi = sd_epi)
}, mc.cores = ncores)
df_axd <- do.call(rbind, results_axd) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1)
  )

arch <- "dxd"
results_dxd <- mclapply(loci_range, function(l) {
  final_epi_reps <- numeric(iter)
  for(i in 1:iter){
    loci.imp <- sort(sample(1:loci, l))
    opt <- l
    sigma <- l/2
    sim_result <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt,
                                      gen, sigma, arch, epi_flag, verbose = FALSE)
    epi_values <- sim_result$lm_arch[,3] - sim_result$lm_arch[,2]
    final_epi_reps[i] <- epi_values[gen]
  }
  mean_epi <- mean(final_epi_reps)
  sd_epi <- sd(final_epi_reps)
  data.frame(loci = l, mean_epi = mean_epi, sd_epi = sd_epi)
}, mc.cores = ncores)
df_dxd <- do.call(rbind, results_dxd) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1)
  )


## PLOT
# Read data
axa <- read.csv("results/df_axa.csv")
axd <- read.csv("results/df_axd.csv")
dxd <- read.csv("results/df_dxd.csv")

# Combine data into a single dataframe with an identifier
axa$group <- "axa"
axd$group <- "axd"
dxd$group <- "dxd"
data_l <- rbind(
  axa[, c("loci", "mean_epi", "ymin", "ymax", "group")], 
  axd[, c("loci", "mean_epi", "ymin", "ymax", "group")], 
  dxd[, c("loci", "mean_epi", "ymin", "ymax", "group")]
)

# Plot the data
colors <- wes_palette("FantasticFox1", 3, type = "continuous")
names(colors) <- c("axa", "axd", "dxd")

ggplot(data_l, aes(x = loci, y = mean_epi, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) + 
  geom_line(linewidth = 1.2) +
  labs(
    x = "Number of Loci",
    y = "Epistatic Variation",
    title = "Effect of Number of Loci on Epistatic Variation"
  ) +
  scale_color_manual(
    values = colors,
    name = "Epistasis Type",
    labels = c("Additive by Additive", "Additive by Dominance", "Dominance by Dominance")
  ) +
  scale_fill_manual(
    values = colors,
    name = "Epistasis Type",
    labels = c("Additive by Additive", "Additive by Dominance", "Dominance by Dominance")
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(min(data_l$loci), max(data_l$loci), length.out = 5)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )






#### SELECTION AND Ve ####

num_sigma <- 100
iter <- 200
sigma_values <- seq(1, 10, length.out = num_sigma)
colors <- wes_palette("FantasticFox1", 3, type = "continuous")
names(colors) <- c("Epistasis", "Phenotype", "Additive")

results <- mclapply(sigma_values, function(sigma_now) {
  final_epi_reps <- numeric(iter)
  for(i in 1:iter){
    sim_result <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt, gen, sigma_now, arch, epi_flag, verbose = FALSE)
    epi_values <- sim_result$lm_arch[,3] - sim_result$lm_arch[,2]
    final_epi_reps[i] <- epi_values[gen]
  }
  mean_epi <- mean(final_epi_reps)
  sd_epi <- sd(final_epi_reps)
  data.frame(sigma = sigma_now, mean_epi = mean_epi, sd_epi = sd_epi)
}, mc.cores = ncores)

results <- readRDS("axa.rds")
df_summary <- do.call(rbind, results) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1)
  )
ggplot(df_summary, aes(x = sigma, y = mean_epi)) +
  geom_line(color = colors["Epistasis"], size = 1.2) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = colors["Epistasis"], alpha = 0.2) +
  labs(
    x = "Strength of Selection",
    y = "Epistatic Variation",
    title = "Additive by Additive"
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    breaks = seq(1, 10, 2)
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "right",
    plot.title = element_text(size = 20),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
results <- readRDS("axd.rds")
df_summary <- do.call(rbind, results) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1)
  )
ggplot(df_summary, aes(x = sigma, y = mean_epi)) +
  geom_line(color = colors["Phenotype"], size = 1.2) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = colors["Phenotype"], alpha = 0.2) +
  labs(
    x = "Strength of Selection",
    y = "Epistatic Variation",
    title = "Additive by Dominance"
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    breaks = seq(1, 10, 2)
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "right",
    plot.title = element_text(size = 20),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
results <- readRDS("dxd.rds")
df_summary <- do.call(rbind, results) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1)
  )
ggplot(df_summary, aes(x = sigma, y = mean_epi)) +
  geom_line(color = colors["Additive"], size = 1.2) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = colors["Additive"], alpha = 0.2) +
  labs(
    x = "Strength of Selection",
    y = "Epistatic Variation",
    title = "Dominance by Dominance"
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    breaks = seq(1, 10, 2)
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "right",
    plot.title = element_text(size = 20),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )









# Read and process each RDS file --------------------------------------------

# For "Additive by Additive" (axa)
results_axa <- readRDS("results/axa.rds")
df_axa <- do.call(rbind, results_axa) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1),
    group = "axa"
  )

# For "Additive by Dominance" (axd)
results_axd <- readRDS("results/axd.rds")
df_axd <- do.call(rbind, results_axd) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1),
    group = "axd"
  )

# For "Dominance by Dominance" (dxd)
results_dxd <- readRDS("results/dxd.rds")
df_dxd <- do.call(rbind, results_dxd) %>%
  mutate(
    ymin = pmax(mean_epi - sd_epi, 0),
    ymax = pmin(mean_epi + sd_epi, 1),
    group = "dxd"
  )

# Combine the datasets -------------------------------------------------------

# Here we assume all data frames have the same x-variable name (here "sigma").
# If your data use "loci" instead, update the column name accordingly.
data_s <- bind_rows(
  df_axa[, c("sigma", "mean_epi", "ymin", "ymax", "group")],
  df_axd[, c("sigma", "mean_epi", "ymin", "ymax", "group")],
  df_dxd[, c("sigma", "mean_epi", "ymin", "ymax", "group")]
)

# Define the colors using a continuous palette from wesanderson.
# Note that the names here should match the group labels.
colors <- wes_palette("FantasticFox1", 3, type = "continuous")
names(colors) <- c("axa", "axd", "dxd")

# Plot the combined data ------------------------------------------------------

ggplot(data_s, aes(x = sigma, y = mean_epi, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  labs(
    x = "Strength of Selection",
    y = "Epistatic Variation",
    title = "Epistatic Variation by Strength of Selection"
  ) +
  scale_color_manual(
    values = colors,
    name = "Epistasis Type",
    labels = c("Additive by Additive", "Additive by Dominance", "Dominance by Dominance")
  ) +
  scale_fill_manual(
    values = colors,
    name = "Epistasis Type",
    labels = c("Additive by Additive", "Additive by Dominance", "Dominance by Dominance")
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(min(data_combined$sigma), max(data_combined$sigma), length.out = 5)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )


#### SELECTION AND Ve (3 types) ####
# Parameters
iter <- 5
gen <- 200
sigma <- 5
architectures <- c("axa", "axd", "dxd")

# Placeholder for results
all_results <- list()

for (arch in architectures) {
  # Initialize matrices for the current architecture
  pheno <- add <- dom <- epi <- matrix(NA, nrow = iter, ncol = gen)
  
  # Run simulations for the current architecture
  for (i in 1:iter) {
    print(paste("Architecture:", arch, "- Iteration:", i))
    simulation_result <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt, gen, sigma, arch, epi_flag, verbose = FALSE)
    add[i, ] <- simulation_result$lm_arch[, 1]
    dom[i, ] <- (simulation_result$lm_arch[, 2] - simulation_result$lm_arch[, 1])
    epi[i, ] <- (simulation_result$lm_arch[, 3] - simulation_result$lm_arch[, 2])
    pheno[i, ] <- simulation_result$avg_phenos
  }
  
  # Store results for epistasis
  all_results[[arch]] <- data.frame(
    Generation = 1:gen,
    Epistasis = colMeans(epi, na.rm = TRUE),
    Architecture = arch
  )
}

# Combine results from all architectures into a single data frame
combined_results <- do.call(rbind, all_results)

# Plot Epistasis across architectures
colors <- wes_palette("FantasticFox1", length(architectures), type = "continuous")
names(colors) <- architectures

ggplot(combined_results, aes(x = Generation, y = Epistasis, color = Architecture)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = colors,
    name = "Architecture"
  ) +
  labs(
    x = "Generation",
    y = "Mean Epistasis",
    title = "Epistasis Across Generations for Different Architectures"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "right",
    plot.title = element_text(size = 20),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 14)
  )



#### FITNESS FUNCITON ####
df <- expand.grid(
  sigma = 1:10,
  obs   = seq(0, 20, length.out = 100)
) %>%
  mutate(fitness = GetFit(obs, 10, sigma))

cols <- wes_palette("FantasticFox1", 10, type = "continuous")

ggplot(df, aes(x = obs, y = fitness, color = factor(sigma), group = sigma)) +
  geom_line(size = 1.2) +
  geom_vline(aes(linetype = "Optimal Phenotype"), linetype = "dashed", xintercept = 10, color = "black", size = 0.8) +
  scale_color_manual(
    values = cols,
    name = "Sigma",
    labels = paste0("s = ", 1:10),
    guide = guide_legend(order = 2)
  ) +
  scale_linetype_manual(
    name = "Optimal Phenotype",
    values = "dashed",
    guide = guide_legend(order = 1)
  ) +
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, 5),
    labels = seq(-10, 10, 5)
  ) +
  labs(x = "Deviation from Optimal Phenotype", y = "Fitness",
       title = "Effect of Sigma (Strength of Selection) On Fitness of Phenotypes") +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "right",
    plot.title = element_text(size = 20),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

#### SIMULATION ACROSS GEN ####
iter <- 20
pheno <- add <- dom <- epi <- matrix(NA, nrow = iter, ncol = gen)
for (i in 1:iter) {
  print(i)
  simulation_result <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt, gen, sigma, arch, epi_flag, verbose=F)
  add[i,] <- simulation_result$lm_arch[,1]
  dom[i,] <- (simulation_result$lm_arch[,2]-simulation_result$lm_arch[,1])
  epi[i,] <- (simulation_result$lm_arch[,3]-simulation_result$lm_arch[,2])
  pheno[i,] <- simulation_result$avg_phenos
}

df <- data.frame(
  Generation = 1:gen,
  Epistasis = colMeans(epi, na.rm = TRUE),
  Additive = colMeans(add, na.rm = TRUE),
  Dominance = colMeans(add, na.rm = TRUE),
  Phenotype = colMeans(pheno, na.rm = TRUE) / 20  # Scaling as per original code
)
df_long <- df %>%
  pivot_longer(
    cols = c(Epistasis, Additive, Phenotype),
    names_to = "Component",
    values_to = "Fitness"
  )

colors <- wes_palette("FantasticFox1", 3, type = "continuous")
names(colors) <- c("Epistasis", "Additive", "Phenotype")

ggplot(df_long, aes(x = Generation, y = Fitness, color = Component, linetype = Component)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = colors,
    name = "Values",
    labels = c("Additive", "Epistasis", "Phenotype")
  ) +
  scale_linetype_manual(
    values = c("Additive" = "solid", "Epistasis" = "solid", "Phenotype" = "dashed"),
    name = "Values",
    labels = c("Additive", "Epistasis", "Phenotype")
  ) +
  labs(
    x = "Generation",
    y = "Proportion of Total Variation",
    title = "Simulation Variables Across Generations"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ . * 20, name = "Phenotype")
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "right",
    plot.title = element_text(size = 20),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
