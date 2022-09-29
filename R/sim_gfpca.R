# simulate bfpca data
# simulated data is a function wrapped from refund package


# N: number of subjects
# J: dimension of the data
# if case = 1 true eigenfunctions are based on alternating sine and cos
# if case = 2 true eigenfunctions are based on sqrt functions
bfpca_sim <- function(N = 500, J = 100, case = 1){
  sind <- (1:J)/J
  K    <- 4 #number of eigenfunctions
  lambdaTrue <- c(1,0.5,0.5^2,0.5^3) # True eigenvalues

  if(case == 1){
    psi <- sqrt(2)*cbind(sin(2*pi*sind),cos(2*pi*sind),
                         sin(4*pi*sind),cos(4*pi*sind))
  }else if(case == 2){
    psi <- cbind(rep(1,J),sqrt(3)*(2*sind-1),
                 sqrt(5)*(6*sind^2-6*sind+1),
                 sqrt(7)*(20*sind^3-30*sind^2+12*sind-1))
  }

  # simulate \xi_ik
  xi <- matrix(rnorm(N*K),N,K)
  xi <- xi %*% diag(sqrt(lambdaTrue))
  # calculate linear predictor \eta_i(s)
  X <- xi %*% t(psi)

  # store in a matrix
  Y <- matrix(rbinom(N*J, size=1, prob=plogis(X)), N, J, byrow=F)

  df_bfpca <- data.frame("id"=rep(1:N, each=J),
                         "index" = rep(sind, N),
                         "value" = as.vector(t(Y)))

  list(
    df_bfpca = df_bfpca,
    xB = X,
    psi = psi,
    lambda = lambdaTrue
  )
}



