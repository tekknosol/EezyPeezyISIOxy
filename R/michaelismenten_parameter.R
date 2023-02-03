
get_prior <- function(trophy, n = 1){
  if (trophy == 'oligo'){
    
    Flux <- rnorm(n, mean = -0.2731324, sd = 0.1669265)  # (g / m2 / d)
    Khalf <- rnorm(n, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(n, mean = 1.07, sd = 0.03)
    
  } else if (trophy == 'eutro'){
    
    Flux <- rnorm(n, mean = -170.1525, sd = 873.4926)  # (g / m2 / d)
    Khalf <- rnorm(n, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(n, mean = 1.08, sd = 0.03)
    
  }
  return(tibble(Flux = Flux, Khalf = Khalf, Theta = Theta, trophic_state = trophy))
}