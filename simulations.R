# Andres Barboza
# 
# Simulations

library(parallel)

## Starting Conditions
N             <- 1000                                # pop size
loci          <- 100                                 # total number of loci
mu            <- 1e-5                                # mutation rate per gene
baseval       <- 0                                   # base phenotype value
num.imp.loci  <- 10                                  # causal loci for sweeps 1 & 2
loci.imp      <- sort(sample(1:loci, num.imp.loci))  # assigns loci for phenotype
opt           <- num.imp.loci                        # optimal phenotype value
sigma         <- num.imp.loci / 2                    # Sigma for fitness curve
gen           <- 2000                                # generations
epi_flag      <- "alter"                             # "half" or "alter", to pick interacting loci
sag           <- 1.1                                 # exponent for inc/dec epistasis NOT USED YET
arches        <- c("axa", "axd", "dxd")              # list of genetic architectures


# Number of cores for parallel loops
ncores <- 4
iter <- 200

## Single Simulation

sim <- SimulateGenerations(N = N, loci = loci, mu = mu, baseval = baseval,
                           loci.imp = loci.imp, opt = opt, gen = gen,
                           sigma = 1, arch = arches[1], epi_flag = epi_flag)

######## STEP 1: Loci sweep #########
loci_range <- seq(4, 100, by=2)
compute_loci_variance <- function(arch) {
  res <- mclapply(loci_range, function(l) {
    add_v <- dom_v <- epi_v <- numeric(iter)
    for(i in seq_len(iter)){
      loci.imp <- sort(sample(1:loci, l))
      sim      <- SimulateGenerations(N = N, loci = loci, mu = mu, baseval = baseval,
                                      loci.imp = loci.imp,
                                      opt=l, gen = gen, sigma=l/2,
                                      arch = arch, epi_flag = epi_flag)
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




######## STEP 2: Selection strength sweep #########
iter <- 10
sigma_values <- seq(1, 10, length.out=10)
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

