# Gaussian fitness function
obs <- rnorm(1000, mean=10,sd=5)
opt <- 10
# sigma determines the strength of selection
# IMHO it should be adjusted such that the most extreme
# possible phenotypes produce values close to zero
GetFit <- function(obs, opt, sigma){
  numer <- (obs - opt)^2
  denom <- (2 * sigma)^2
  w <- exp(-(numer / denom))
  return(w)
}
w <- GetFit(obs=obs, opt=opt, sigma=1)
plot(w~obs)
