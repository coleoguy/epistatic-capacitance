### simulation_variance_analysis.R

# Load libraries
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(purrr)

# Number of cores for parallel loops
ncores <- 4

###### Starting Conditions #########
# Fixed parameters
N             <- 1000       # initial pop size placeholder
loci          <- 100        # total loci
mu            <- 1e-5       # mutation rate per gene
baseval       <- 0          # base phenotype value
num.imp.loci  <- 10         # # causal loci for sweeps 1 & 2
opt           <- num.imp.loci
sigma         <- num.imp.loci / 2
gen           <- 2000       # generations
epi_flag      <- "alter"  # "half" or "alter"
sag           <- 1.1        # exponent for inc/dec
arches        <- c("axa", "axd", "dxd")
###### End Starting Conditions #########

###### FUNCTION DEFINITIONS #########
GetPopulation <- function(N, loci) {
  matrix(rep(4, N * loci), nrow = N, ncol = loci)
}

MutatePop <- function(pop, mu) {
  mut.prob <- 1 - (1 - mu)^1547
  coords   <- which(matrix(runif(nrow(pop) * ncol(pop)), nrow = nrow(pop)) < mut.prob, arr.ind = TRUE)
  for(i in seq_len(nrow(coords))) {
    r <- coords[i,1]; c <- coords[i,2]
    pop[r,c] <- switch(pop[r,c], sample(c(2,3),1), sample(c(1,4),1), sample(c(1,4),1), sample(c(2,3),1))
  }
  pop
}

GetPheno <- function(pop, loci.imp, baseval, arch, epi_flag) {
  mat <- pop[, loci.imp]
  mat[mat==1]         <- 0
  mat[mat %in% c(2,3)] <- 1
  mat[mat==4]         <- 2
  if(arch == "add") return(rowSums(mat) + baseval)
  get_two_locus <- function(mat) {
    if(epi_flag == "half") {
      half <- ncol(mat)/2
      dir  <- mat[,1:half]; inv <- mat[,(half+1):(2*half)]
    } else {
      dir <- mat[,seq(1,ncol(mat),2)]
      inv <- mat[,seq(2,ncol(mat),2)]
    }
    list(dir=dir, inv=inv)
  }
  if(arch %in% c("axa","dxd","axd")) {
    tl  <- get_two_locus(mat)
    dir <- tl$dir; inv <- tl$inv
  }
  if(arch == "axa") {
    c1 <- (dir==2 & inv==2)|(dir==0 & inv==0)
    c2 <- (dir==2 & inv==0)|(dir==0 & inv==2)
    c3 <- (dir==1|inv==1)
    return(baseval + rowSums(4*c1 + 0*c2 + 2*c3))
  }
  if(arch == "dxd") {
    c1 <- (dir %in% c(2,0) & inv!=1)|(dir!=1 & inv %in% c(2,0))
    c2 <- ( (dir==1 & inv!=1)|(dir!=1 & inv==1) )
    c3 <- (dir==1 & inv==1)
    return(baseval + rowSums(2*c1 + 0*c2 + 4*c3))
  }
  if(arch == "axd") {
    c1 <- (dir==2 & inv!=1)
    c2 <- (dir==1)
    c3 <- (dir==0 & inv!=1)
    c4 <- (dir==2 & inv==1)
    c5 <- (dir==0 & inv!=1)
    return(baseval + rowSums(1*c1 + 2*c2 + 3*c3 + 4*c4 + 0*c5))
  }
  if(arch == "inc") {
    return((ncol(mat)*2)*((rowSums(mat)+baseval)/(ncol(mat)*2))^sag)
  }
  if(arch == "dec") {
    return(baseval + (ncol(mat)*2-baseval)*((rowSums(mat))/(ncol(mat)*2-baseval))^(1/sag))
  }
}

GetFit <- function(obs, opt, sigma) {
  exp(-((obs-opt)^2)/(2*sigma)^2)
}

Reproduction <- function(pop, N, w, loci) {
  parents <- sample(nrow(pop), 2*N, TRUE, prob = w)
  reco    <- pmax(pmin(rpois(2*N,10), loci),1)
  haps    <- t(sapply(reco, function(r) rep(sample(1:2,r,TRUE), each=ceiling(loci/r), length.out=loci)))
  mat     <- pop[parents,]
  gametes <- ifelse(mat==1, 0,
                    ifelse(mat==4, 1,
                           ifelse(haps==1, 0, 1)))
  m_gametes <- gametes[seq(1,2*N,2),]
  p_gametes <- gametes[seq(2,2*N,2),]
  1 + m_gametes + 2*p_gametes
}

GetArch <- function(pop, loci.imp, phenos) {
  pop2   <- pop[, loci.imp]
  a_mat  <- matrix(0, nrow=nrow(pop2), ncol=ncol(pop2))
  d_mat  <- matrix(1, nrow=nrow(pop2), ncol=ncol(pop2))
  a_mat[pop2==1] <- 1; a_mat[pop2==4] <- -1
  d_mat[pop2 %in% c(1,4)] <- 0
  R2a  <- summary(lm(phenos~a_mat))$adj.r.squared
  R2ad <- summary(lm(phenos~a_mat+d_mat))$adj.r.squared
  half <- if(epi_flag=="alter") ncol(pop2)/2 else NULL
  idx  <- if(epi_flag=="alter") seq(1,ncol(pop2),2) else 1:half
  e_aa <- a_mat[,idx]*a_mat[,idx+ifelse(epi_flag=="alter",1,half)]
  e_ad <- a_mat[,idx]*d_mat[,idx+ifelse(epi_flag=="alter",1,half)]
  e_dd <- d_mat[,idx]*d_mat[,idx+ifelse(epi_flag=="alter",1,half)]
  R2ade <- summary(lm(phenos~a_mat+d_mat+e_aa+e_ad+e_dd))$adj.r.squared
  c(R2a, R2ad, R2ade)
}

ApplyBottleneck <- function(pop, prop) {
  survivors <- pop[
    sample(nrow(pop), size = ceiling(nrow(pop) * prop), replace = FALSE),
    , drop = FALSE
  ]
  return(survivors)
}

SimulateGenerations <- function(
    N, loci, mu, baseval, loci.imp,
    opt, gen, sigma, arch, epi_flag,
    do_bottleneck   = FALSE,
    bottleneck_prop = 0.1,
    do_opt_change   = FALSE,
    gen_post        = 200,
    verbose         = FALSE
) {
  pre_gen   <- gen
  total_gen <- pre_gen + if (do_bottleneck) gen_post else 0
  event_gen <- if (do_bottleneck) pre_gen else floor(pre_gen / 2)
  
  pop         <- GetPopulation(N, loci)
  opt_current <- opt
  mean_pheno  <- numeric(total_gen)
  opt_history <- numeric(total_gen)
  lm_arch     <- matrix(NA, nrow = total_gen, ncol = 3)
  
  for (g in seq_len(total_gen)) {
    # 1) mutation
    pop <- MutatePop(pop, mu)
    
    # 2) bottleneck (only once, at generation = pre_gen)
    if (do_bottleneck && g == event_gen) {
      pop <- ApplyBottleneck(pop, bottleneck_prop)
      N   <- nrow(pop)  # recovery draws from survivors
    }
    
    # 3) optional optimal shift at the same event
    if (do_opt_change && g == event_gen) {
      opt_current <- sample(c(baseval, length(loci.imp) * 2), 1)
    }
    
    # 4) phenotype, fitness & record
    phenos        <- GetPheno(pop, loci.imp, baseval, arch, epi_flag)
    mean_pheno[g] <- mean(phenos)
    opt_history[g]<- opt_current
    lm_arch[g, ]  <- GetArch(pop, loci.imp, phenos)
    
    # 5) reproduction back up to N
    w   <- GetFit(phenos, opt_current, sigma)
    pop <- Reproduction(pop, N, w, loci)
    
    if (verbose) message("Gen ", g)
  }
  
  # 6) build the post-event deviation vector
  diff_vec     <- abs(mean_pheno - opt_history)
  mean_response <- mean(
    diff_vec[(event_gen + 1):(event_gen + gen_post)]
  )
  
  list(
    final_population = pop,
    avg_phenos       = mean_pheno,
    lm_arch          = lm_arch,
    diff_vec         = diff_vec,
    mean_response    = mean_response
  )
}



###### Helper: tidy long #########
tidy_variance <- function(df, xvar) {
  df %>%
    pivot_longer(
      cols = matches("^(mean|sd)_"),
      names_to = c(".value","component"),
      names_pattern = "(mean|sd)_(.*)"
    ) %>%
    mutate(
      mean      = pmax(mean, 0),              # clamp negatives
      ymin      = pmax(mean - sd, 0),
      ymax      = pmin(mean + sd, 1),
      component = factor(component,
                         levels=c("add","dom","epi"),
                         labels=c("Additive","Dominance","Epistasis"))
    ) %>%
    rename(!!xvar := all_of(xvar))
}

# Palettes
pal_arch <- wes_palette("FantasticFox1", length(arches), type="continuous")
names(pal_arch) <- arches
# define color palette for the three components
pal_comp <- wes_palette("FantasticFox1", 3, type = "continuous")
names(pal_comp) <- c("Additive", "Dominance", "Epistasis")


######## STEP 1: Loci sweep #########
iter       <- 200
loci_range <- seq(4, 100, by=2)
compute_loci_variance <- function(arch) {
  res <- mclapply(loci_range, function(l) {
    add_v <- dom_v <- epi_v <- numeric(iter)
    for(i in seq_len(iter)){
      loci.imp <- sort(sample(1:loci, l))
      sim      <- SimulateGenerations(N, loci, mu, baseval, loci.imp,
                                      opt=l, gen, sigma=l/2,
                                      arch, epi_flag)
      vals     <- sim$lm_arch[gen,]
      add_v[i] <- vals[1]
      dom_v[i] <- vals[2] - vals[1]
      epi_v[i] <- vals[3] - vals[2]
    }
    data.frame(
      loci     = l,
      mean_add = mean(add_v), sd_add = sd(add_v),
      mean_dom = mean(dom_v), sd_dom = sd(dom_v),
      mean_epi = mean(epi_v), sd_epi = sd(epi_v)
    )
  }, mc.cores = ncores)
  df <- bind_rows(res)
  df$architecture <- arch
  df
}
df_loci <- bind_rows(lapply(arches, compute_loci_variance))
write.csv(df_loci, "df_loci.csv", row.names=FALSE)

loci_long <- tidy_variance(df_loci, "loci")

p1 <- ggplot(loci_long, aes(x = loci, y = mean,
                            color = component, fill = component)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~architecture, nrow = 1) +
  labs(
    x     = "Number of Loci",
    y     = "Genetic Variance",
    title = "Effect of Number of Loci on Genetic Variance",
    color = "Variance\nComponent",
    fill  = "Variance\nComponent"
  ) +
  scale_color_manual(values = pal_comp) +
  scale_fill_manual(values = pal_comp) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(min(loci_range), max(loci_range), length.out = 5)) +
  theme_minimal() +
  theme(
    panel.border    = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "bottom",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(size = 20),
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 14)
  )

######## STEP 2: Selection strength sweep #########
sigma_values <- seq(1, 10, length.out=100)
compute_sigma_variance <- function(arch) {
  res <- mclapply(sigma_values, function(s) {
    add_v <- dom_v <- epi_v <- numeric(iter)
    for(i in seq_len(iter)){
      sim  <- SimulateGenerations(N, loci, mu, baseval,
                                  sort(sample(1:loci, num.imp.loci)), opt,
                                  gen, sigma=s, arch, epi_flag)
      vals <- sim$lm_arch[gen,]
      add_v[i] <- vals[1]
      dom_v[i] <- vals[2] - vals[1]
      epi_v[i] <- vals[3] - vals[2]
    }
    data.frame(
      sigma    = s,
      mean_add = mean(add_v), sd_add = sd(add_v),
      mean_dom = mean(dom_v), sd_dom = sd(dom_v),
      mean_epi = mean(epi_v), sd_epi = sd(epi_v)
    )
  }, mc.cores = ncores)
  df <- bind_rows(res)
  df$architecture <- arch
  df
}
df_sigma <- bind_rows(lapply(arches, compute_sigma_variance))
write.csv(df_sigma, "df_sigma.csv", row.names=FALSE)

sigma_long <- tidy_variance(df_sigma, "sigma")

p2 <- ggplot(sigma_long, aes(x = sigma, y = mean,
                             color = component, fill = component)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~architecture, nrow = 1) +
  labs(
    x     = "Selection Strength (σ)",
    y     = "Genetic Variance",
    title = "Effect of Selection Strength on Genetic Variance",
    color = "Variance\nComponent",
    fill  = "Variance\nComponent"
  ) +
  scale_color_manual(values = pal_comp) +
  scale_fill_manual(values = pal_comp) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(1, 10, by = 2)) +
  theme_minimal() +
  theme(
    panel.border    = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "bottom",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(size = 20),
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 14)
  )
######## STEP 3: Population size sweep #########
N_values <- c(100, 250, 500, 1000, 1500, 2000)

compute_N_variance <- function(arch) {
  # for each population size…
  res <- lapply(N_values, function(N_now) {
    # launch 'iter' replicates in parallel
    replicate_results <- mclapply(seq_len(iter), function(i) {
      sim  <- SimulateGenerations(
        N_now, loci, mu, baseval,
        sort(sample(1:loci, 20)),
        opt = 20, gen, sigma = 5,
        arch, epi_flag,
        verbose = FALSE
      )
      vals <- sim$lm_arch[gen, ]
      # return named vector of the three components
      c(
        add = vals[1],
        dom = vals[2] - vals[1],
        epi = vals[3] - vals[2]
      )
    }, mc.cores = ncores)
    
    # bind into a matrix and compute means + sds
    mat <- do.call(rbind, replicate_results)
    data.frame(
      N        = N_now,
      mean_add = mean(mat[ , "add"]), sd_add = sd(mat[ , "add"]),
      mean_dom = mean(mat[ , "dom"]), sd_dom = sd(mat[ , "dom"]),
      mean_epi = mean(mat[ , "epi"]), sd_epi = sd(mat[ , "epi"])
    )
  })
  
  # stitch together and tag the architecture
  df <- bind_rows(res)
  df$architecture <- arch
  df
}

# run it
df_N <- bind_rows(lapply(arches, compute_N_variance))
write.csv(df_N, "df_N.csv", row.names = FALSE)

# and tidy up for plotting
N_long <- tidy_variance(df_N, "N")


p3 <- ggplot(N_long, aes(x = N, y = mean,
                         color = component, fill = component)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~architecture, nrow = 1) +
  labs(
    x     = "Population Size (N)",
    y     = "Genetic Variance",
    title = "Effect of Population Size on Genetic Variance",
    color = "Variance\nComponent",
    fill  = "Variance\nComponent"
  ) +
  scale_color_manual(values = pal_comp) +
  scale_fill_manual(values = pal_comp) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = N_values) +
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

######## STEP 4: Bottleneck effect on response to selection #########

iter            <- 4
bottle_prop     <- 0.1
arches          <- c("axa","axd","dxd")
do_opt_change   <- TRUE
bottleneck_flag <- c(FALSE, TRUE)

gen      <- 200
gen_post <- 100

df_time <- bind_rows(
  lapply(arches, function(arch) {
    lapply(bottleneck_flag, function(do_bottle) {
      # run `iter` replicates in parallel
      res_list <- mclapply(seq_len(iter), function(i) {
        sim <- SimulateGenerations(
          N, loci, mu, baseval,
          sort(sample(1:loci, num.imp.loci)),
          opt, gen, sigma, arch, epi_flag,
          do_bottleneck   = do_bottle,
          bottleneck_prop = bottle_prop,
          do_opt_change   = do_opt_change,
          gen_post        = gen_post,
          verbose         = FALSE
        )
        
        # determine where the event happened
        event_gen <- if (do_bottle) gen else floor(gen / 2)
        
        data.frame(
          architecture = arch,
          bottleneck   = do_bottle,
          time         = seq_len(gen_post),
          diff         = sim$diff_vec[
            (event_gen + 1):(event_gen + gen_post)
          ]
        )
      }, mc.cores = ncores)
      
      bind_rows(res_list)
    }) %>% bind_rows()
  }) %>% bind_rows()
)

# summarise
df_time_sum <- df_time %>%
  group_by(architecture, bottleneck, time) %>%
  summarise(
    mean_diff = mean(diff),
    min_diff  = min(diff),
    max_diff  = max(diff),
    .groups   = "drop"
  )

pal_b <- c("FALSE" = "blue", "TRUE" = "red")

p4 <- ggplot(df_time_sum,
       aes(x = time, y = mean_diff,
           color = factor(bottleneck),
           fill  = factor(bottleneck))) +
  geom_ribbon(aes(ymin = min_diff, ymax = max_diff),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~architecture, nrow = 1) +
  scale_color_manual(
    name   = "Bottleneck",
    values = pal_b,
    labels = c("No", "Yes")
  ) +
  scale_fill_manual(
    name   = "Bottleneck",
    values = pal_b,
    labels = c("No", "Yes")
  ) +
  labs(
    x     = "Generations after optimum shift",
    y     = expression("|mean phenotype – optimum|"),
    title = "Response to selection with vs. without bottleneck"
  ) +
  theme_minimal() +
  theme(
    panel.border    = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "bottom",
    plot.title      = element_text(size = 18),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12)
  )

# Print plots
print(p1)
print(p2)
print(p3)
print(p4)
