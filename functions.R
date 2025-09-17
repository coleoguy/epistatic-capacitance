# Andres Barboza
# 
# Functions

GetPopulation <- function(N, loci) {
  # pop <- matrix(sample(1:4, N*loci, replace=T), N, loci)
  pop <- matrix(rep(4, N * loci), nrow = N, ncol = loci)
  return(pop)
}

MutatePop <- function(pop, mu) {
  mut.prob <- 1 - (1 - mu)^1547
  coords   <- which(matrix(runif(nrow(pop) * ncol(pop)), nrow = nrow(pop)) < mut.prob, arr.ind = TRUE)
  for(i in seq_len(nrow(coords))) {
    r <- coords[i,1]; c <- coords[i,2]
    pop[r,c] <- switch(pop[r,c], sample(c(2,3),1), sample(c(1,4),1), sample(c(1,4),1), sample(c(2,3),1))
  }
  return(pop)
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
    c5 <- (dir==0 & inv==1)
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
  w <- exp(-((obs-opt)^2)/(2*sigma)^2)
  return(w)
}

Reproduction <- function(pop, N, w, loci) {
  parents <- sample(nrow(pop), 2*N, TRUE, prob = w)
  unique  <- length(unique(parents))
  reco    <- pmax(pmin(rpois(2*N,10), loci),1)
  haps    <- t(sapply(reco, function(r) rep(sample(1:2,r,TRUE), each=ceiling(loci/r), length.out=loci)))
  mat     <- pop[parents,]
  gametes <- ifelse(mat==1, 0,
                    ifelse(mat==4, 1,
                           ifelse(haps==1, 0, 1)))
  m_gametes <- gametes[seq(1,2*N,2),]
  p_gametes <- gametes[seq(2,2*N,2),]
  npop <- 1 + m_gametes + 2*p_gametes
  out <- list(npop, unique)
  return(out)
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
  archs <- c(R2a, R2ad, R2ade)
  return(archs)
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
  
  pop         <- GetPopulation(N = N, loci = loci)
  opt_current <- opt
  mean_pheno  <- numeric(total_gen)
  opt_history <- numeric(total_gen)
  lm_arch     <- matrix(NA, nrow = total_gen, ncol = 3)
  
  for (g in seq_len(total_gen)) {
    # 1) mutation
    pop <- MutatePop(pop = pop, mu = mu)
    
    # 2) bottleneck (only once, at generation = pre_gen)
    if (do_bottleneck && g == event_gen) {
      pop <- ApplyBottleneck(pop = pop, bottleneck_prop = bottleneck_prop)
      N   <- nrow(pop)  # recovery draws from survivors
    }
    
    # 3) optional optimal shift at the same event
    if (do_opt_change && g == event_gen) {
      opt_current <- sample(c(baseval, length(loci.imp) * 2), 1)
    }
    
    # 4) phenotype, fitness & record
    phenos        <- GetPheno(pop = pop, loci.imp = loci.imp, baseval = baseval,
                              arch = arch, epi_flag = epi_flag)
    mean_pheno[g] <- mean(phenos)
    opt_history[g]<- opt_current
    lm_arch[g, ]  <- GetArch(pop = pop, loci.imp = loci.imp, phenos = phenos)
    
    # 5) reproduction back up to N
    w      <- GetFit(obs =phenos, opt = opt_current, sigma = sigma)
    reprod <- Reproduction(pop = pop, N = N, w = w, loci = loci)
    pop    <- reprod[[1]]
    
    # store vector of unique parents
    uniq.parents <- c()
    uniq.parents[g] <- reprod[[2]]
    
    if (verbose) message("Gen ", g)
  }
  
  # 6) build the post-event deviation vector
  diff_vec     <- abs(mean_pheno - opt_history)
  mean_response <- mean(
    diff_vec[(event_gen + 1):(event_gen + gen_post)]
  )
  
  sim <- list(
    final_population = pop,
    avg_phenos       = mean_pheno,
    lm_arch          = lm_arch,
    diff_vec         = diff_vec,
    mean_response    = mean_response,
    uniq.parents    = reprod[[2]]
  )
  return(sim)
}
