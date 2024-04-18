MutatePop <- function(pop, mu){
  mut.prob <- 1-(1-mu)^1547 # 1547 is the number of base pairs in the coding region of our genome (treating it as the length of a gene)
  mutations.sites <- which(matrix(runif(nrow(pop) * ncol(pop)), nrow = nrow(pop), ncol = ncol(pop)) <   mut.prob, arr.ind = TRUE)
  # not finished but will work on it
  Mutate <- function(index) {
    switch(pop[index,2],
           sample(c(2,3), 1),  #for 1
           sample(c(1,4), 1),  #for 2
           sample(c(1,4), 1),  #for 3
           sample(c(2,3), 1))  #for 4
  }
  
  mutants <- sapply(mutations[,1], Mutate)
}



