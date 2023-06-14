
get_prior <- function(trophy, n = 1){
  sdf <- .7
  
  if (trophy == 'oligo'){
    
    # Flux <- rnorm(n, mean = -0.2731324, sd = 0.1669265)  # (g / m2 / d)
    # Flux <- -rlnorm(n, meanlog = -1.247376, sdlog = 0.01631389)  # (g / m2 / d)

    #truncated
    Flux <- -rlnorm(n, meanlog = -1.688881, sdlog = 1.474197 * sdf)  # (g / m2 / d)

    #truncated < 0.3
    # Flux <- -rlnorm(n, meanlog = -1.919972, sdlog = 1.680658)  # (g / m2 / d)

    # original
    # Flux <- -rlnorm(n, meanlog = -1.547376, sdlog = 1.417046)  # (g / m2 / d)

    # Flux <- -rlnorm(n, meanlog = -5.247376, sdlog = 0.0041631389)  # (g / m2 / d)
    # Flux <- rnorm(n, mean = -0.32, sd = 0.096)  # (g / m2 / d)
    Khalf <- rnorm(n, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(n, mean = 1.07, sd = 0.03)
    
  } else if (trophy == 'eutro'){
    
    # Flux <- abs(rnorm(n, mean = -170.1525, sd = 873.4926*1e-1))  # (g / m2 / d)
    # Flux <- -rlnorm(n, meanlog = 1.082377, sdlog = 0.0346636)  # (g / m2 / d)
    
    #truncated
    Flux <- -rlnorm(n, meanlog = 0.6602776, sdlog = 1.502307 * sdf)  # (g / m2 / d)

    #truncated > 0.45
    # Flux <- -rlnorm(n, meanlog = 0.8511423, sdlog = 1.486916)  # (g / m2 / d)

    #original
    # Flux <- -rlnorm(n, meanlog = 0.9098668, sdlog = 2.2648594)  # (g / m2 / d)
    
    # Flux <- -rlnorm(n, meanlog = -0.4558641, sdlog = 0.8448238)  # (g / m2 / d)
    # Flux <- -rlnorm(n, meanlog = -0.4558641, sdlog = 0.08448238)  # (g / m2 / d)
    # Flux <- rnorm(n, mean = -3.2, sd = 0.096)
    Khalf <- rnorm(n, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(n, mean = 1.08, sd = 0.03)
    
  } else if (trophy == 'eutro2'){
    Flux <- -rlnorm(n, meanlog = 0.9098668, sdlog = 0.4346636)  # (g / m2 / d)
    # Flux <- -rlnorm(n, meanlog = -0.4558641, sdlog = 0.08448238)  # (g / m2 / d)
    # Flux <- rnorm(n, mean = -3.2, sd = 0.096)
    Khalf <- rnorm(n, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(n, mean = 1.08, sd = 0.03)
  }
  # Flux <- sapply(Flux, FUN=min, 0)
  return(tibble(Flux = Flux, Khalf = Khalf, Theta = Theta, trophic_state = trophy))
}
