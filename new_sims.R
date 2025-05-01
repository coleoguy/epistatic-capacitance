## Andres Barboza - andresdbp@tamu.edu
## Zoya Wani - zoyawani1@tamu.edu
## April 16, 2024

library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)

ncores <- 4

###### Starting Conditions #########
N            <- 500
loci         <- 100
mu           <- 1e-5
baseval      <- 0
num.imp.loci <- 10
opt          <- num.imp.loci
sigma        <- num.imp.loci/2
gen          <- 500
arch         <- "axd"
epi_flag     <- "alter"
sag          <- 1.1
###### End Starting Conditions #########

###### FUNCTIONS #########

GetPopulation <- function(N,loci){
  
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

#### 1) NUMBER OF LOCI AND Ve ####

iter       <- 4
loci_range <- seq(4, 100, by = 2)
arches     <- c("axa", "axd", "dxd")

compute_loci_variance <- function(arch) {
  df_list <- mclapply(loci_range, function(l) {
    add_reps <- numeric(iter)
    dom_reps <- numeric(iter)
    epi_reps <- numeric(iter)
    for (i in seq_len(iter)) {
      loci.imp <- sort(sample(1:loci, l))
      sim     <- SimulateGenerations(N, loci, mu, baseval, loci.imp, opt = l,
                                     gen, sigma = l/2, arch, epi_flag, FALSE)
      vals    <- sim$lm_arch[gen, ]
      add_reps[i] <- vals[1]
      dom_reps[i] <- vals[2] - vals[1]
      epi_reps[i] <- vals[3] - vals[2]
    }
    data.frame(
      loci     = l,
      mean_add = mean(add_reps), sd_add = sd(add_reps),
      mean_dom = mean(dom_reps), sd_dom = sd(dom_reps),
      mean_epi = mean(epi_reps), sd_epi = sd(epi_reps)
    )
  }, mc.cores = ncores)
  df <- bind_rows(df_list)
  df$architecture <- arch
  return(df)
}

df_loci <- bind_rows(lapply(arches, compute_loci_variance))
write.csv(df_loci,  "df_loci.csv",  row.names = FALSE)

#### 2) SELECTION STRENGTH AND Ve ####

sigma_values <- seq(1, 10, length.out = 100)

compute_sigma_variance <- function(arch) {
  df_list <- mclapply(sigma_values, function(s) {
    add_reps <- numeric(iter)
    dom_reps <- numeric(iter)
    epi_reps <- numeric(iter)
    for (i in seq_len(iter)) {
      sim   <- SimulateGenerations(N, loci, mu, baseval, 
                                   loci.imp = sort(sample(1:loci, num.imp.loci)),
                                   opt, gen, sigma = s, arch, epi_flag, FALSE)
      vals  <- sim$lm_arch[gen, ]
      add_reps[i] <- vals[1]
      dom_reps[i] <- vals[2] - vals[1]
      epi_reps[i] <- vals[3] - vals[2]
    }
    data.frame(
      sigma    = s,
      mean_add = mean(add_reps), sd_add = sd(add_reps),
      mean_dom = mean(dom_reps), sd_dom = sd(dom_reps),
      mean_epi = mean(epi_reps), sd_epi = sd(epi_reps)
    )
  }, mc.cores = ncores)
  df <- bind_rows(df_list)
  df$architecture <- arch
  return(df)
}

df_sigma <- bind_rows(lapply(arches, compute_sigma_variance))
write.csv(df_sigma, "df_sigma.csv", row.names = FALSE)

#### 3) TIDY & PLOT ####

# helper to go long
tidy_variance <- function(df, xvar) {
  df %>%
    pivot_longer(
      cols = matches("^(mean|sd)_"),
      names_to = c(".value", "component"),
      names_pattern = "(mean|sd)_(.*)"
    ) %>%
    # clamp any negative means up to 0
    mutate(
      mean      = pmax(mean, 0.0000000001),
      ymin      = pmax(mean - sd, 0),
      ymax      = pmin(mean + sd, 1),
      component = factor(
        component,
        levels = c("add", "dom", "epi"),
        labels = c("Additive", "Dominance", "Epistasis")
      )
    ) %>%
    rename(!!xvar := all_of(xvar))
}


# tidy both datasets
loci_long  <- tidy_variance(df_loci,  "loci")
sigma_long <- tidy_variance(df_sigma, "sigma")

# shared colors by architecture
pal <- wes_palette("FantasticFox1", 3, type = "continuous")
names(pal) <- arches

# A) Number of loci plot
p1 <- ggplot(loci_long, aes(x = loci, y = mean,
                            color = architecture,
                            fill  = architecture)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax),
              color = NA, alpha = 0.2) +
  geom_line(size = 1.2) +
  facet_wrap(~component, nrow = 1) +
  labs(
    x     = "Number of Loci",
    y     = "Genetic Variance",
    title = "Effect of Number of Loci on Genetic Variance",
    color = "Genetic Architecture",
    fill  = "Genetic Architecture"
  ) +
  scale_color_manual(values = pal,
                     labels = c("Additive by Additive",
                                "Additive by Dominance",
                                "Dominance by Dominance")) +
  scale_fill_manual(values = pal,
                    labels = c("Additive by Additive",
                               "Additive by Dominance",
                               "Dominance by Dominance")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(min(loci_range),
                                  max(loci_range),
                                  length.out = 5)) +
  theme_minimal() +
  theme(
    panel.border   = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right","bottom"),
    legend.title   = element_text(size = 16),
    legend.text    = element_text(size = 14),
    plot.title     = element_text(size = 20),
    axis.title     = element_text(size = 16),
    axis.text      = element_text(size = 14)
  )

# B) Selection strength plot
p2 <- ggplot(sigma_long, aes(x = sigma, y = mean,
                             color = architecture,
                             fill  = architecture)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax),
              color = NA, alpha = 0.2) +
  geom_line(size = 1.2) +
  facet_wrap(~component, nrow = 1) +
  labs(
    x     = "Strength of Selection (σ)",
    y     = "Genetic Variance",
    title = "Effect of Selection Strength on Genetic Variance",
    color = "Genetic Architecture",
    fill  = "Genetic Architecture"
  ) +
  scale_color_manual(values = pal,
                     labels = c("Additive by Additive",
                                "Additive by Dominance",
                                "Dominance by Dominance")) +
  scale_fill_manual(values = pal,
                    labels = c("Additive by Additive",
                               "Additive by Dominance",
                               "Dominance by Dominance")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1,10, by = 2)) +
  theme_minimal() +
  theme(
    panel.border   = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "bottom",
    legend.title   = element_text(size = 14),
    legend.text    = element_text(size = 12),
    plot.title     = element_text(size = 20),
    axis.title     = element_text(size = 16),
    axis.text      = element_text(size = 14)
  )

# finally, print or save
print(p1)
print(p2)
# ggsave("loci_vs_variance.png", p1, width=12, height=4)
# ggsave("sigma_vs_variance.png", p2, width=12, height=4)
