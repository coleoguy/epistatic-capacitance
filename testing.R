# Andres Barboza
# 
# Testing










### OLD SCRIPTS

## Additional Scripts

# ---- Modified SimulateGenerations function to track parents ---- #
SimulateGenerationsTrackParents <- function(
    N, loci, mu, baseval, loci.imp,
    opt, gen, sigma, arch, epi_flag
) {
  pop             <- GetPopulation(N, loci)
  opt_current     <- opt
  parent_counts   <- numeric(gen)
  
  for (g in seq_len(gen)) {
    print(g)
    # Mutation
    pop <- MutatePop(pop, mu)
    
    # Phenotype and fitness
    phenos <- GetPheno(pop, loci.imp, baseval, arch, epi_flag)
    w      <- GetFit(phenos, opt_current, sigma)
    
    # Sample parents and count unique ones
    parents <- sample(nrow(pop), 2*N, TRUE, prob = w)
    parent_counts[g] <- length(unique(parents))
    
    # Reproduction
    reco    <- pmax(pmin(rpois(2*N,10), loci),1)
    haps    <- t(sapply(reco, function(r) rep(sample(1:2,r,TRUE), each=ceiling(loci/r), length.out=loci)))
    mat     <- pop[parents,]
    gametes <- ifelse(mat==1, 0,
                      ifelse(mat==4, 1,
                             ifelse(haps==1, 0, 1)))
    m_gametes <- gametes[seq(1,2*N,2),]
    p_gametes <- gametes[seq(2,2*N,2),]
    pop       <- 1 + m_gametes + 2*p_gametes
  }
  
  parent_counts
}

# ---- Parameters ---- #
N           <- 2000
loci        <- 100
mu          <- 1e-5
baseval     <- 0
num.imp.loci<- 50
loci.imp    <- sort(sample(1:loci, num.imp.loci))
opt         <- num.imp.loci
gen         <- 2000
epi_flag    <- "alter"
arches      <- "axa"
sigma_vals  <- c(1, 2, 4, 6, 8, 10)
ncores      <- 6

# ---- Run simulations ---- #
results_list <- mclapply(sigma_vals, function(sigma) {
  parent_counts <- SimulateGenerationsTrackParents(
    N, loci, mu, baseval, loci.imp=sort(sample(1:loci, num.imp.loci)),
    opt, gen, sigma, arch = "dxd", epi_flag
  )
  data.frame(
    Generation   = 1:gen,
    UniqueParents = parent_counts,
    Sigma         = factor(sigma)
  )
}, mc.cores = ncores)

df_parents <- bind_rows(results_list)

# ---- Plot ---- #
plotParents <- ggplot(df_parents, aes(x = Generation, y = UniqueParents, color = Sigma)) +
  geom_line(size = 1.2) +
  labs(
    title = "Number of Unique Parents per Generation",
    x = "Generation",
    y = "Unique Parents",
    color = expression(sigma)
  ) +
  scale_y_continuous(limits = c(0, 2000)) +
  scale_x_continuous(limits = c(0, 2000)) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 18),
    axis.title    = element_text(size = 14),
    axis.text     = element_text(size = 12),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12),
    panel.border  = element_rect(color = "darkgray", fill = NA, size = 1)
  )

print(plotParents)

