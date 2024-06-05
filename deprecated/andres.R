MutatePop <- function(pop, mu){
  mut.prob <- 1-(1-mu)^1547 # 1547 is the number of base pairs in the coding region of our genome (treating it as the length of a gene)
  mut.coord <- which(matrix(runif(nrow(pop) * ncol(pop)), nrow = nrow(pop), ncol = ncol(pop)) < mut.prob, arr.ind = TRUE)
  apply_switch <- function(pop, coordinates) {
    # Create a function to apply the switch statement to a single element
    apply_switch_single <- function(row_index, col_index) {
      value <- pop[row_index, col_index]
      pop[row_index, col_index] <- switch(value,
                                          sample(c(2,3), 1),  #for 1
                                          sample(c(1,4), 1),  #for 2
                                          sample(c(1,4), 1),  #for 3
                                          sample(c(2,3), 1))  #for 4
      return(pop[row_index, col_index])  # return the updated value
    }
    # Apply the function to each row of coordinates and collect results
    updated_values <- mapply(apply_switch_single, coordinates[,1], coordinates[,2])
    # Update the pop matrix with the updated values
    pop[coordinates] <- updated_values
    
    return(pop)
  }
  mutants <- apply_switch(pop, mut.coord)
  
  return(mutants)
}


