# Heath Blackmon
# May 30 2024
# coleoguy@gmail.com
# script to compare response to selection between
# true generating architecture and inferred architecture


# get sim results
fulldat <- read.csv("../results/sim.results.csv")
starts <- seq(from=1, by=8, length.out=6000)
# cycle through all experiments
for(i in 1:6000){
  curdat <- fulldat[starts[i]:(starts[i]+7),]
  # create genome for for simulation:
  # 20 biallelic unlinked loci
  # 0 = [0,0]
  # 1 = [0,1] phase doesn't matter since unlinked?
  # 2 = [1,1]
  # loci 1-10 - trait one
  # loci 11-20 - trait two
  pop <- matrix(sample(0:2, 21000, replace = T), 1000, 20)
  # half female
  pop[1:500, 1] <- 0
  # half male
  pop[501:1000, 1] <- 1
  T.meanpheno <- c()
  # begin simulation for one trait pop of the true architecutre
  for(j in 1:100){
    ## phenotype
    trait <- GetPhenoTrueArch(pop, curdat, "t3")
    if(j==1){
      opt <- trait[[2]]
    }
    w <- GetFit(obs = trait[[1]], opt=opt, sigma=2)
    parents <- sample(1:1000, size=2000, prob = w, replace=T)
    newpop <- pop
    pairs <- (1:2000)[(1:2000) %% 2 == 1]
    for(k in 1:nrow(newpop)){
      newpop[k,] <- GetOffspring(pop[parents[pairs[k]],], pop[parents[pairs[k]+1],])
    }
    pop <- newpop
    T.meanpheno[j] <- mean(GetPhenoTrueArch(pop, curdat, "t3")[[1]])
  }
  
}

# create genome for inferred

# cycle through for 100 generations
## measure fitness
## select parents
## make next generation


# store results as a list of 2x100 matrices with phenotypes


###### Functions #########
GetPhenoTrueArch <- function(pop, curdat, trait){
  ### elemental 1
  # set base level of trait to twice beta to avoid negative phenotypes
  base <- 2 * curdat[1, 2]
  # calculating to be able to know max possible value for t3
  maxt1 <- base + curdat[1, 2] 
  # calculate the impact of genotype on trait
  t1 <- base + curdat[1, 2] * rowSums(pop[,1:10])/20
  ### elemental 2
  # set base level of trait to twice beta to avoid negative phenotypes
  base <- 2 * curdat[2, 2]
  # calculating to be able to know max possible value for t3
  mint2 <- base  
  maxt3 <- maxt1/mint2
  # calculate the impact of genotype on trait
  t2 <- base + curdat[2, 2] * rowSums(pop[,11:20])/20
  ### compound traits
  t3 <- t1/t2
  return(list(t3, maxt3))
}

GetFit <- function(obs, opt, sigma){
  numer <- (obs - opt)^2
  denom <- (2 * sigma)^2
  w <- exp(-(numer / denom))
  return(w)
}

GetOffspring <- function(parent1, parent2){
  gam <- parent1
  gam[gam==1] <- sample(x = 0:1, size=sum(gam==1), replace=T) 
  gam[gam==2] <- 1
  gam2 <- parent2
  gam2[gam2==1] <- sample(x = 0:1, size=sum(gam2==1), replace=T) 
  gam2[gam2==2] <- 1
  offspring <- gam + gam2
  return(offspring)
}
###### End Functions #########










cmat <- matrix(NA, 5, 5)
cmat[,1] <- c(1.0, 0.5, 0.0, -0.5,-1.0)
cmat[,2] <- c(0.0, 0.5, 1.0,  0.5, 0.0)
cmat[,3] <- cmat[,1] * cmat[,1]
cmat[,4] <- cmat[,1] * cmat[,2]
cmat[,5] <- cmat[,2] * cmat[,2]
row.names(cmat) <- c("P1", "BC1", "F1", "BC2", "P2")

