#' Simulate GFPCA
#'
#' This function generates simulated generalized functional principal components data.
#' Simulation setup is adapted from code in the help file for the function fpca.face
#' from the refund package. 4 true principal components are generated.
#'
#' @param N Number of subjects.
#' @param J Number of timepoints per subject.
#' @param case Takes on values of 1 or 2. If case = 1 true eigenfunctions
#' are based on alternating sine and cosines. If case = 2 true eigenfunctions are
#' based on sqrt functions.
#' @param family Defines exponential family for generating data. Options are binomial (default),
#' poisson, or gaussian.
#' @param sigma Defaults to 2. Error variance for Y when family = "gaussian".
#'
#' @importFrom stats plogis binomial rnorm rbinom
#'
#' @export
#'
#' @examples
#' # generate binary data
#' df <- sim_gfpca(N = 200, J = 100, case = 1, family = "binomial")$df_gfpca
#'
#' # generate Poisson data
#' df <- sim_gfpca(N = 200, J = 100, case = 1, family = "poisson")$df_gfpca
#'
#' # generate Gaussian data
#' df <- sim_gfpca(N = 200, J = 100, case = 1, family = "gaussian")$df_gfpca
#'
sim_gfpca <- function(N = 500,
                      J = 100,
                      case = 1,
                      family = "binomial",
                      sigma = 2){
  sind <- (1:J)/J
  K    <- 4 #number of eigenfunctions
  lambdaTrue <- c(1,0.5,0.5^2,0.5^3) # True eigenvalues

  if(case == 1){
    phi <- sqrt(2)*cbind(sin(2*pi*sind),cos(2*pi*sind),
                         sin(4*pi*sind),cos(4*pi*sind))
  }else if(case == 2){
    phi <- cbind(rep(1,J),sqrt(3)*(2*sind-1),
                 sqrt(5)*(6*sind^2-6*sind+1),
                 sqrt(7)*(20*sind^3-30*sind^2+12*sind-1))
  }

  # simulate \xi_ik
  xi <- matrix(rnorm(N*K),N,K)
  xi <- xi %*% diag(sqrt(lambdaTrue))
  # calculate linear predictor \eta_i(s)
  X <- xi %*% t(phi)

  # store in a matrix
  # change this part for simulating data from other EF families
  if(family == "binomial"){
    Y <- matrix(rbinom(N*J, size=1, prob=plogis(X)), N, J, byrow=F)
  }else if(family == "poisson"){
    Y <- matrix(rpois(N*J, lambda=exp(X)), N, J, byrow=F)
  }else if(family == "gaussian"){
    Y <- X + sigma*matrix(rnorm(N*J),N,J)
  }

  df_gfpca <- data.frame(id = rep(1:N, each=J),
                         index = rep(seq(0, 1, length.out = J), N),
                         value = as.vector(t(Y)),
                         eta = as.vector(t(X)))

  list(
    df_gfpca = df_gfpca,
    phi = phi,
    lambda = lambdaTrue,
    scores = xi
  )
}



