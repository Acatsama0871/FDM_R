# Auxiliary functions
# Define payoff functions
# The call payoff
# The simple payoff function for call options
# Args:
#  strike: the strike price of the option
#  spot: the spot price
# Return:
#  The payoff of the call option
Payoff_call <- function(strike, spot) {
  return(max(spot - strike, 0.0))
}

# The put payoff
# The simple payoff function for the put option
# Args:
#  strike: the strike price of the option
#  spot: the spot price
# Return:
#  The payoff of the put option
Payoff_put <- function(strike, spot) {
  return(max(strike - spot, 0.0))
}


# Define the pre-final valuation function
# The pre-final valuation function for the European function
# Args:
#  strike: the strike price of the opton
#  futureValue: the discounted expectation for the future nodes
#  spot: the current spot price
#  thePayoff: the payoff function
# Return:
#  the option price for current node
European_preFinal <- function(futureValue, strike, spot, thePayoff) {
  return(futureValue)
}

# The pre-final valuation function for the American function
# Args:

#  strike: the strike price of the option
#  futureValue: the discounted expectation for the future nodes
#  spot: the current spot price
#  thePayoff: the payoff function
# Return:
#  the option price for current node
American_preFinal <- function(futureValue, strike, spot, thePayoff) {
  return(max(thePayoff(strike, spot), futureValue))
}



# Finite difference method
# The explicit finite difference method
# Args:
#  k: strike price
#  t: time to maturity
#  s: spot price
#  r: risk-free rate
#  sig: volatility
#  q: dividen yield
#  N: divide time to maturity to N intervals
#  Nj: number of dx
# Return:
#  the option price
Explicit_FDM.Eur.call <- function(k, t, s, r, sig, q, N, Nj) {
  # Precompute constants
  dt <- t / N
  dx <- sig * sqrt(3 * dt)
  log_s <- log(s)
  a = r - q - 0.5 * sig^2
  pu <- dt * (sig^2 / (2 * dx^2) + a / (2 * dx))
  pm <- 1- dt * sig^2 / dx^2 - r * dt
  pd <- dt * (sig^2 / (2 * dx^2) - a / (2 * dx))
  v <- matrix(data = 0, nrow = 2 * Nj + 1, ncol = N)
  
  # Initialise the prices at the maturity
  st <- vector(mode = "double", length = 2 * Nj + 1) # Create an empty vector
  st[length(st)] <- log_s - Nj * dx
  for (i in (length(st) - 1):1) {
    st[i] <- st[i + 1] + dx
  }
  st <- exp(st)
  
  # Initialise option values at maturity
  v[, N] <- mapply(Payoff_call, strike = k, spot = st)
  
  # Iterate back through lattice
  for (j in (N - 1):1) {
    for (i in 2:(2 * Nj)) {
      v[i, j] <- pu * v[i - 1, j + 1] + pm * v[i, j + 1] + pd * v[i + 1, j + 1]
    }
    
    # Boundary condition
    v[1, j] <- v[2, j] + exp(dx)
    v[(2 * Nj + 1), j] <- v[(2 * Nj), j]
  }
  
  return(v[Nj + 1, 1])
}



Explicit_FDM.Eur.call(k = 65, t = 0.25, s = 60, r = 0.08, sig = 0.3, q = 0, N = 10000, Nj = 30000)
