build_HKY_Q_matrix <- function(frequencies, kappa) {
  # frequencies is a named vector with elements A, C, G, T
  # kappa is the transition/transversion rate ratio
  
  # Ensure frequencies sum to 1
  if (sum(frequencies) != 1) {
    print("Frequencies are being normalized 1")
    frequencies <- frequencies / sum(frequencies)
  }
  
  # Ensure all frequencies are named and in the order A, C, G, T
  if (!all(names(frequencies) %in% c("A", "C", "G", "T"))) {
    stop("Frequencies vector must have names A, C, G, T")
  }
  
  # Extract nucleotide frequencies
  pi_A <- frequencies["A"]
  pi_C <- frequencies["C"]
  pi_G <- frequencies["G"]
  pi_T <- frequencies["T"]
  
  # Initialize the Q matrix
  Q <- matrix(0, 4, 4)
  rownames(Q) <- colnames(Q) <- c("A", "C", "G", "T")
  
  # Fill the Q matrix with the appropriate rates
  Q["A", "C"] <- pi_C
  Q["A", "G"] <- kappa * pi_G
  Q["A", "T"] <- pi_T
  
  Q["C", "A"] <- pi_A
  Q["C", "G"] <- pi_G
  Q["C", "T"] <- kappa * pi_T
  
  Q["G", "A"] <- kappa * pi_A
  Q["G", "C"] <- pi_C
  Q["G", "T"] <- pi_T
  
  Q["T", "A"] <- pi_A
  Q["T", "C"] <- kappa * pi_C
  Q["T", "G"] <- pi_G
  
  # Set the diagonal elements
  for (i in 1:4) {
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


# rates <- ratesSimul
# eigen(getQ(rates))
# eigen(t(getQ(rates)))
# branchLength=0.02865092600135935
# exp_sub=0.001
# expm::expm(exp_sub*branchLength*getQ(rates))


ratesToEigen <- function(rates, branchLength=1, exp_sub=1, normalize=FALSE) {
  cat("rates=", rates, "\n")
  cat("logRates=", log(rates))
  Q <- getQ(rates)
  print(Q)
  eigenQ <- eigen(Q)
  if (normalize == TRUE) norm = - sum(diag(Q) * c(0.33333, 0.66666)) else norm = 1
  print(eigenQ)
  print(norm)
  return(expm::expm(exp_sub * branchLength * Q / norm))
}


hhhhhhh <- function(x) {
  

ratesAna <- c(1.2370433983690592, 0.7228641596600762, 0.8518582935889486, 1.4070408050329466, 1.18306706956355, 0.8962899369542581) 
ratesToEigen(ratesAna)
ratesSimul <- exp(c(0, .6931471805, 0, 0, .6931471805, 0)) 
branchLength <- 0.02865092600135935
ratesToEigen(ratesSimul, branchLength)

a <- "0 .6931471805 1.0986122886 1.3862943611 1.6094379124 1.7917594692 1.9459101490 0 .6931471805 1.0986122886 1.3862943611 1.6094379124 1.7917594692 0 .6931471805 1.0986122886 1.3862943611 1.6094379124 0 .6931471805 1.0986122886 1.3862943611 0 .6931471805 1.0986122886 0 .6931471805 0 0 .6931471805 0 1.0986122886 .6931471805 0 1.3862943611 1.0986122886 .6931471805 0 1.6094379124 1.3862943611 1.0986122886 .6931471805 0 1.7917594692 1.6094379124 1.3862943611 1.0986122886 .6931471805 0 1.9459101490 1.7917594692 1.6094379124 1.3862943611 1.0986122886 .6931471805 0"
as.numeric(strsplit(a, " ")[[1]]) |> exp() |>  getQ()



ratesGlm <- c(0.24999999027997286, 1.2500000608572392, 0.24999999027997286, 0.24999999027997286, 1.2500000608572392, 0.24999999027997286, 0.24999999027997286, 1.2500000608572392, 0.24999999027997286, 0.24999999027997286, 1.2500000608572392, 0.24999999027997286)

B <- eigen(getQ(ratesGlm))
A <- matrix(c(-2.3852447794681098E-17, 0.7071067811865475, -0.5000000000000002, 0.5000000000000002, 0.7071067811865476, 1.1102230246251565E-16, 0.5000000000000002, 0.4999999999999997, -2.3852447794681098E-17, -0.7071067811865476, -0.5000000000000002, 0.5000000000000002, -0.7071067811865476, 0.0, 0.5000000000000003, 0.5000000000000001), 4, 4, byrow=TRUE) 

getQ(ratesGlm) %*% A[,1] * eigen(getQ(ratesGlm))$values[3]

A %*% diag(c(-3.000000102274425, -3.0000001022744245, -0.9999999611198924, 0.0)) %*% A

A %*% diag(c(-3.000000102274425, -3.0000001022744245,  -0.9999999611198924, 0.0)) %*% solve(A)
# the third and forth are in the right location

B$vec %*% diag(B$val) %*% solve(B$vec)

cbind(getQ(ratesGlm)%*%A[,1], A[,1]*(-3.000000102274425))

ratesToEigen(as.vector(rep(1,12)))

# $values
# -3.000000102274425, -3.0000001022744245,  -0.9999999611198924, 0.0
# $vectors
# [,1]          [,2] [,3] [,4]
# [1,] -2.385245e-17  7.071068e-01 -0.5  0.5
# [2,]  7.071068e-01  1.110223e-16  0.5  0.5
# [3,] -2.385245e-17 -7.071068e-01 -0.5  0.5
# [4,] -7.071068e-01  0.000000e+00  0.5  0.5



}