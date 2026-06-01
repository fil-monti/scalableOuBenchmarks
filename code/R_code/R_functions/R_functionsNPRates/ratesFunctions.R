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
  for (i in 1:n) {
    for (j in 1:i) {
      if (i != j) {
        Q_vec[counter] <- Q[i, j]
        counter <- counter + 1
      }
    }
  }
  return(Q_vec)
}