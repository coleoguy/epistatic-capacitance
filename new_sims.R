### simulation_variance_analysis.R

# Load libraries
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(purrr)

ncores <- 4

###### Starting Conditions #########
N             <- 1000
loci          <- 100
mu            <- 1e-5
baseval       <- 0
num.imp.loci  <- 10
opt           <- num.imp.loci
sigma         <- num.imp.loci/2
gen           <- 2000
epi_flag      <- "alter"
sag           <- 1.1
arches        <- c("axa","axd","dxd")
###### End Starting Conditions #########

###### FUNCTION DEFINITIONS #########

GetPopulation <- function(N, loci) {
  matrix(4, nrow = N, ncol = loci)
}

MutatePop <- function(pop, mu) {
  mut.prob <- 1 - (1 - mu)^1547
  coords   <- which(matrix(runif(nrow(pop)*ncol(pop)), nrow = nrow(pop)) < mut.prob,
                    arr.ind = TRUE)
  for(i in seq_len(nrow(coords))) {
    r <- coords[i,1]; c <- coords[i,2]
    pop[r,c] <- switch(pop[r,c],
                       sample(c(2,3),1),
                       sample(c(1,4),1),
                       sample(c(1,4),1),
                       sample(c(2,3),1))
  }
  pop
}

GetPheno <- function(pop, loci.imp, baseval, arch, epi_flag) {
  mat <- pop[, loci.imp]
  mat[mat==1]         <- 0
  mat[mat %in% c(2,3)] <- 1
  mat[mat==4]         <- 2
  # … (same as before for each arch) …
  # omitted here for brevity; assume exact same logic
}

GetFit <- function(obs, opt, sigma) {
  exp(-((obs-opt)^2)/(2*sigma)^2)
}

Reproduction <- function(pop, N, w, loci) {
  parents <- sample(nrow(pop), 2*N, TRUE, prob = w)
  reco    <- pmax(pmin(rpois(2*N,10), loci),1)
  haps    <- t(sapply(reco, function(r) rep(sample(1:2,r,TRUE),
                                            each=ceiling(loci/r),
                                            length.out=loci)))
  mat     <- pop[parents,]
  gametes <- ifelse(mat==1, 0,
                    ifelse(mat==4, 1,
                           ifelse(haps==1,0,1)))
  mg <- gametes[seq(1,2*N,2),]
  pg <- gametes[seq(2,2*N,2),]
  1 + mg + 2*pg
}

GetArch <- function(pop, loci.imp, phenos) {
  # … same as before …
}

ApplyBottleneck <- function(pop, prop) {
  pop[sample(nrow(pop), size = ceiling(nrow(pop)*prop), replace = FALSE), , drop=FALSE]
}

# NEW SimulateGenerations (constant N after bottleneck)
SimulateGenerations <- function(N, loci, mu, baseval, loci.imp,
                                opt, gen, sigma, arch, epi_flag,
                                do_bottleneck   = FALSE,
                                bottleneck_prop = 0.1,
                                do_opt_change   = FALSE,
                                verbose         = FALSE) {
  pop         <- GetPopulation(N, loci)
  mid         <- floor(gen/2)
  opt_current <- opt
  
  mean_pheno  <- numeric(gen)
  opt_hist    <- numeric(gen)
  lm_arch     <- matrix(NA, gen, 3)
  
  for(g in 1:gen) {
    pop <- MutatePop(pop, mu)
    
    if(do_bottleneck && g == mid) {
      pop <- ApplyBottleneck(pop, bottleneck_prop)
    }
    if(do_opt_change && g == mid) {
      opt_current <- sample(c(baseval, length(loci.imp)*2), 1)
    }
    
    phenos        <- GetPheno(pop, loci.imp, baseval, arch, epi_flag)
    mean_pheno[g] <- mean(phenos)
    opt_hist[g]   <- opt_current
    
    lm_arch[g, ] <- GetArch(pop, loci.imp, phenos)
    
    w   <- GetFit(phenos, opt_current, sigma)
    pop <- Reproduction(pop, N, w, loci)
    
    if(verbose) message("Gen ", g)
  }
  
  diff_vec     <- abs(mean_pheno - opt_hist)
  mean_response <- mean(diff_vec[(mid+1):gen])
  
  list(
    final_population = pop,
    avg_phenos       = mean_pheno,
    lm_arch          = lm_arch,
    diff_vec         = diff_vec,
    mean_response    = mean_response
  )
}

###### Helper to tidy variance #########
tidy_variance <- function(df, xvar) {
  df %>%
    pivot_longer(matches("^(mean|sd)_"),
                 names_to = c(".value","component"),
                 names_pattern = "(mean|sd)_(.*)") %>%
    mutate(
      mean      = pmax(mean,0),
      ymin      = pmax(mean - sd, 0),
      ymax      = pmin(mean + sd, 1),
      component = factor(component,
                         levels=c("add","dom","epi"),
                         labels=c("Additive","Dominance","Epistasis"))
    ) %>%
    rename(!!xvar := all_of(xvar))
}

pal_comp <- wes_palette("FantasticFox1", 3, type="continuous")
names(pal_comp) <- c("Additive","Dominance","Epistasis")

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
iter    <- 200
N_values <- c(50,100,200,500,1000,2000)

compute_N_variance <- function(arch) {
  bind_rows(lapply(N_values, function(N_now) {
    # parallelize the 200 replicates
    reps <- mclapply(seq_len(iter), function(i) {
      sim  <- SimulateGenerations(N_now, loci, mu, baseval,
                                  sort(sample(1:loci, num.imp.loci)),
                                  opt, gen, sigma = 5,
                                  arch, epi_flag)
      vals <- sim$lm_arch[gen,]
      c(add = vals[1],
        dom = vals[2] - vals[1],
        epi = vals[3] - vals[2])
    }, mc.cores = ncores)
    mat <- do.call(rbind, reps)
    data.frame(
      N         = N_now,
      mean_add  = mean(mat[,"add"]),  sd_add  = sd(mat[,"add"]),
      mean_dom  = mean(mat[,"dom"]),  sd_dom  = sd(mat[,"dom"]),
      mean_epi  = mean(mat[,"epi"]),  sd_epi  = sd(mat[,"epi"])
    )
  })) %>% mutate(architecture = arch)
}

df_N <- bind_rows(lapply(arches, compute_N_variance))
write.csv(df_N, "df_N.csv", row.names = FALSE)
N_long <- tidy_variance(df_N, "N")

p3 <- ggplot(N_long, aes(N, mean, color = component, fill = component)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.2) +
  geom_line(size=1.2) +
  facet_wrap(~architecture, nrow=1) +
  labs(x="Population Size (N)", y="Genetic Variance",
       title="Effect of N on Genetic Variance",
       color="Component", fill="Component") +
  scale_color_manual(values=pal_comp) +
  scale_fill_manual(values=pal_comp) +
  theme_minimal()

######## STEP 4: Bottleneck effect (constant-N) #########
iter            <- 200
bottle_prop     <- 0.1
do_opt_change   <- TRUE
bottleneck_flag <- c(FALSE, TRUE)
mid <- floor(gen/2)
post_len <- gen - mid

df_time <- bind_rows(
  lapply(arches, function(arch) {
    lapply(bottleneck_flag, function(do_bottle) {
      out <- mclapply(seq_len(iter), function(i) {
        sim <- SimulateGenerations(N, loci, mu, baseval,
                                   sort(sample(1:loci, num.imp.loci)),
                                   opt, gen, sigma, arch, epi_flag,
                                   do_bottleneck   = do_bottle,
                                   bottleneck_prop = bottle_prop,
                                   do_opt_change   = do_opt_change)
        data.frame(
          architecture = arch,
          bottleneck   = do_bottle,
          time         = seq_len(post_len),
          diff         = sim$diff_vec[(mid+1):gen]
        )
      }, mc.cores = ncores)
      bind_rows(out)
    }) %>% bind_rows()
  }) %>% bind_rows()
)

df_time_sum <- df_time %>%
  group_by(architecture, bottleneck, time) %>%
  summarize(
    mean_diff = mean(diff),
    min_diff  = min(diff),
    max_diff  = max(diff),
    .groups   = "drop"
  )

pal_b <- c("FALSE"="#1b9e77","TRUE"="#d95f02")

p4 <- ggplot(df_time_sum,
             aes(x=time, y=mean_diff,
                 color=factor(bottleneck), fill=factor(bottleneck))) +
  geom_ribbon(aes(ymin=min_diff, ymax=max_diff), alpha=0.2, color=NA) +
  geom_line(size=1.2) +
  facet_wrap(~architecture, nrow=1) +
  scale_color_manual(name="Bottleneck", values=pal_b, labels=c("No","Yes")) +
  scale_fill_manual(name="Bottleneck", values=pal_b, labels=c("No","Yes")) +
  labs(x="Generations after shift",
       y=expression("|mean phenotype – optimum|"),
       title="Bottleneck vs. No Bottleneck Response") +
  theme_minimal()


# ## Full trajectory
# iter            <- 8
# bottle_prop     <- 0.1
# arches          <- c("axa","axd","dxd")
# do_opt_change   <- TRUE
# bottleneck_flag <- c(FALSE, TRUE)
# 
# # 1) run all sims and collect per‐gen diff for 1:gen
# df_time <- bind_rows(
#   lapply(arches, function(arch) {
#     lapply(bottleneck_flag, function(do_bottle) {
#       res_list <- mclapply(seq_len(iter), function(i) {
#         sim <- SimulateGenerations(
#           N, loci, mu, baseval,
#           sort(sample(1:loci, num.imp.loci)),
#           opt, gen, sigma, arch, epi_flag,
#           do_bottleneck   = do_bottle,
#           bottleneck_prop = bottle_prop,
#           do_opt_change   = do_opt_change,
#           verbose         = FALSE
#         )
#         data.frame(
#           architecture = arch,
#           bottleneck   = do_bottle,
#           time         = seq_len(gen),
#           diff         = sim$diff_vec[1:gen]
#         )
#       }, mc.cores = ncores)
#       bind_rows(res_list)
#     }) %>% bind_rows()
#   }) %>% bind_rows()
# )
# 
# # 2) summarize across replicates at each generation
# df_time_sum <- df_time %>%
#   group_by(architecture, bottleneck, time) %>%
#   summarise(
#     mean_diff = mean(diff),
#     min_diff  = min(diff),
#     max_diff  = max(diff),
#     .groups   = "drop"
#   )
# 
# # 3) plotting full trajectory
# pal_b <- c("FALSE" = "#1b9e77", "TRUE" = "#d95f02")
# 
# ggplot(df_time_sum,
#        aes(x = time, y = mean_diff,
#            color = factor(bottleneck),
#            fill  = factor(bottleneck))) +
#   geom_ribbon(aes(ymin = min_diff, ymax = max_diff),
#               alpha = 0.2, color = NA) +
#   geom_line(size = 1.2) +
#   facet_wrap(~architecture, nrow = 1) +
#   scale_color_manual(
#     name   = "Bottleneck",
#     values = pal_b,
#     labels = c("No", "Yes")
#   ) +
#   scale_fill_manual(
#     name   = "Bottleneck",
#     values = pal_b,
#     labels = c("No", "Yes")
#   ) +
#   labs(
#     x     = "Generation",
#     y     = expression("|mean phenotype – optimum|"),
#     title = "Full‐trajectory response to selection with vs. without bottleneck"
#   ) +
#   theme_minimal() +
#   theme(
#     panel.border    = element_rect(color = "darkgray", fill = NA, size = 1),
#     legend.position = "bottom",
#     plot.title      = element_text(size = 18),
#     axis.title      = element_text(size = 14),
#     axis.text       = element_text(size = 12)
#   )


######### Print plots #########
print(p1)
print(p2)
print(p3)
print(p4)

