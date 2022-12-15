
get_prior <- function(trophy, n = 1){
  if (trophy == 'oligo'){
    
    Flux <- rnorm(n, mean = -0.32, sd = 0.096)  # (g / m2 / d)
    Khalf <- rnorm(n, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(n, mean = 1.07, sd = 0.03)
    
  } else if (trophy == 'eutro'){
    
    Flux <- rnorm(n, mean = -3.2, sd = 0.096)  # (g / m2 / d)
    Khalf <- rnorm(n, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(n, mean = 1.08, sd = 0.03)
    
  }
  return(tibble(Flux = Flux, Khalf = Khalf, Theta = Theta, trophic_state = trophy))
}