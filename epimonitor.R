N <- 100 #population
loci <- 100 #positions on the genome 
mu <- 10^-5 #human mutation rate 10^-9 for an individual nucleotide
baseval <- 0 # this is a base minimum value for our phenotype
loci.imp <- sort(sample(2:loci, loci/10))
opt <- 20
sigma <- 5
gen <- 50
arch <- "sign" # add, sign, inc, dec
sag <- 2
sign_flag <- "half" # half, alter


pop <- matrix(data = sample(1:4, N*loci, replace=T), N, loci)
phenos <- GetPheno(pop, loci.imp, baseval, arch, sign_flag)
pop <- pop[,loci.imp]
pop[,1:5] <- 1
a_mat <- matrix(data=rep(0, length(pop)), N, length(loci.imp))
a_mat[pop == 1] <- +1
a_mat[pop == 4] <- -1
fit <- lm(phenos ~ a_mat)
R2a <-summary(fit)$adj.r.squared

d_mat <- matrix(data=rep(1, length(pop)), N, length(loci.imp))
d_mat[pop == 1] <- 0
d_mat[pop == 4] <- 0
fit <- lm(phenos ~ a_mat + d_mat)
R2ad <-summary(fit)$adj.r.squared


e_mat_aa <- matrix(data=NA, N, length(loci.imp)/2)
curcol <- 1
if(sign_flag == "alter"){
  for(i in seq(from=1, to=ncol(pop),by=2)){
    e_mat_aa[,curcol] <- a_mat[ ,i] * a_mat[,i+1]
    curcol <- curcol + 1
  }
}else{
  for(i in 1:(ncol(pop)/2)){
    e_mat_aa[,curcol] <- a_mat[ ,i] * a_mat[,i+(ncol(pop)/2)]
    curcol <- curcol + 1
  }
}
#TODO figure out if this AD vs DA is an issue
e_mat_ad <- matrix(data=NA, N, length(loci.imp)/2)
curcol <- 1
if(sign_flag == "alter"){
  for(i in seq(from=1, to=ncol(pop),by=2)){
    e_mat_ad[,curcol] <- a_mat[ ,i] * d_mat[,i+1]
    curcol <- curcol + 1
  }
}else{
  for(i in 1:(ncol(pop)/2)){
    e_mat_ad[,curcol] <- a_mat[ ,i] * d_mat[,i+(ncol(pop)/2)]
    curcol <- curcol + 1
  }
}
e_mat_dd <- matrix(data=NA, N, length(loci.imp)/2)
curcol <- 1
if(sign_flag == "alter"){
  for(i in seq(from=1, to=ncol(pop),by=2)){
    e_mat_dd[,curcol] <- d_mat[ ,i] * d_mat[,i+1]
    curcol <- curcol + 1
  }
}else{
  for(i in 1:(ncol(pop)/2)){
    e_mat_dd[,curcol] <- d_mat[ ,i] * d_mat[,i+(ncol(pop)/2)]
    curcol <- curcol + 1
  }
}
fit <- lm(phenos ~ a_mat + d_mat + e_mat_aa + e_mat_ad + e_mat_dd)
R2ade <-summary(fit)$adj.r.squared





