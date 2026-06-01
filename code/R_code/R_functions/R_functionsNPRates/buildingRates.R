getRates <- function(Q) {
  n <- nrow(Q)
  n_rates <- n * (n - 1)
  Q_vec <- rep(0, n_rates)
  
  counter <- 1
  for (i in 1:n) {
    for (j in i:n) {
      if (i != j) {
        Q_vec[counter] <- Q[i, j]
        counter <- counter + 1
      }
    }
  }
  for (j in 1:n) {
    for (i in j:n) {
      if (i != j) {
        Q_vec[counter] <- Q[i, j]
        counter <- counter + 1
      }
    }
  }
  return(Q_vec)
}

getQ <- function(rates) {
  n <- 0.5 * (1 + sqrt(1 + 4 * length(rates)))
  Q <- matrix(0, n, n)
  counter <- 1
  for (i in 1:n) {
    for (j in i:n) {
      if (i != j) {
        Q[i, j] <- rates[counter]
        counter <- counter + 1
      }
    }
  }
  for (j in 1:n) {
    for (i in j:n) {
      if (i != j) {
        Q[i, j] <- rates[counter]
        counter <- counter + 1
      }
    }
  }
  for (i in 1:n) {
    Q[i, i] <- -sum(Q[i, -i])
  }
  return(Q)
}
# # Example usage
# frequencies <- c(A = 1, C = 1, G = 1, T = 1)
# kappa <- 5.0
# 
# Q <- build_HKY_Q_matrix(frequencies, kappa)
# print(Q)
# log(getRates(Q))

# for (i in 3:8) {
#   x <- 1:i
#   outer(x, x, function(x, y) abs(x - y)) |> getRates() |> print()
# }
if ("5"=="4") nstates <- 6
computeRatesAlph <- function(nstates, predictorName) {
  if ( typeof(nstates) == "character" ) nstates <- nstates %>% as.numeric()
  x <- 1:nstates
  if ( predictorName == "L1") {
    Q <- outer(x, x, function(x, y) 1/abs(x - y)) 
  } else if ( predictorName == "L2") {
    Q <- outer(x, x, function(x, y) 1/(x - y)^2) 
  } else if ( predictorName == "L1Circle") {
    Q <- outer(x, x, function(x, y) 1/pmin(abs(x - y), nstates - abs(x - y))) 
  } else if ( predictorName == "L2Circle") {
    Q <- 2*outer(x, x, function(x, y) 1/pmin((x - y)^2, (nstates - abs(x - y))^2)) 
  } else {
    stop("Unknown predictor name")
  }
  return(getRates(Q))
}

# Q has not diagonal elements!!
# rates <- Q |> getRates()
# cat(rates)



