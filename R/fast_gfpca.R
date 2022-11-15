#' Fast generalized functional principal components analysis
#'
#' This is the main function for the fastGFPCA package. The function requires input data \code{Y} to be a dataframe in long format with variables
#' \code{id}, \code{index}, and \code{value} to indicate subject IDs,
#' observation times on the domain, and observations, respectively.
#' The \code{index} must contain the same, equally spaced grid points for each subject.

#' The number of functional principal components (FPCs) can either be specified
#' directly (argument \code{npc}) or chosen based on the explained share of
#' variance (\code{pve}). Families and link functions can be specified using the same
#' syntax as the \code{family} argument from the \code{lme4::glmer()} function.
#'
#'
#' @author Andrew Leroux \email{andrew.leroux@@cuanschutz.edu},
#' Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' @import dplyr
#' @importFrom stats coef predict binomial lm median as.formula
#' @importFrom refund fpca.face
#' @importFrom lme4 glmer
#' @importFrom utils txtProgressBar setTxtProgressBar data
#' @importFrom mgcv bam predict.bam
#'
#' @return An object of class \code{fpca} containing:
#' \item{efunctions}{\eqn{D \times npc} matrix of estimated FPC basis functions.}
#' \item{evalues}{Estimated variance of the FPC scores.}
#' \item{npc}{number of FPCs.}
#' \item{scores}{\eqn{I \times npc} matrix of estimated FPC scores.}
#' \item{mu}{Estimated population-level mean.}
#' \item{Yhat}{FPC approximation of subject-specific means, before applying the
#' response function.}
#' \item{Y}{The observed data.}
#' \item{family}{for compatibility with \code{refund.shiny} package.}
#' @export
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.

#' @examples
#' # simulate data
#' set.seed(1001)
#'
#' # binomial data, with bins that do not overlap
#' df_gfpca <- sim_gfpca(N = 200, J = 200, case = 2)$df_gfpca
#' gfpca_mod <- fast_gfpca(df_gfpca, overlap = FALSE, binwidth = 10, family = binomial)
#'
#' # Poisson data, overlapping bins
#' df_gfpca <- sim_gfpca(N = 100, J = 100, case = 1, family = "poisson")$df_gfpca
#' gfpca_mod <- fast_gfpca(df_gfpca, overlap = TRUE, binwidth =6, family = "poisson")
#'
#' @param Y Dataframe. Should have variables id, value, index.
#' @param argvals numeric; grid over which functions are observed.  If null defaults to unique values of index.
#' @param overlap Logical; indicates whether or not to construct overlapping bins. Defaults to FALSE
#' @param binwidth controls the width of the bins for step 1. Must have an even integer value (defaults to 10). Can be no more than J/10.
#' @param pve proportion of variance explained: used to choose the number of
#' principal components unless `npc` is specified.
#' @param npc how many smooth PCs to try to extract, if \code{NULL} (the
#' default) then npc is chosen based on `pve`
#' @param family exponential family to be passed to \code{glmer} and \code{bam}.
#' @param periodicity Option for a periodic spline basis. Defaults to FALSE.
#' @param ... Additional arguments passed to or from other functions
#'@export

fast_gfpca <- function(Y,
                       argvals = NULL, # grid for functional observations
                       overlap = FALSE,
                       binwidth = 10, # must be even number
                       pve = NULL,
                       npc = NULL,
                       family,
                       periodicity = FALSE,
                       ...){


  if(!all(c("id", "value", "index") %in% ls(Y))){
    stop('Y must be a dataframe containing variables "id", "index", and "value".')
  }

  # sort data by id and index
  Y = Y[order(Y$id, Y$index),]

  if(is.null(npc) & is.null(pve)){
    pve = 0.95
    message("Setting pve = 0.95")
  }

  N <- length(unique(Y$id))
  J <- length(unique(Y$index)) # assumes all subjects are on same even grid

  # add check that binwidth is even. If not, it will bet converted to an even number
  if(J/binwidth < 10){
    binwidth = J/10
    message(paste0("binwidth should be no more than J/10. Converting to a new binwidth of ", binwidth, "."))
  }
  if((binwidth %% 2) != 0) {
    binwidth <- 2 * round(binwidth/2)
    message(paste0("binwidth should have an even integer value. Converting to a new binwidth of ", binwidth, "."))
  }


  if(is.null(argvals)){
    argvals <- sort(unique(Y$index))
  }

  if(overlap){

    fit_fastgfpca <- vector(mode="list",length=J)
    #pb <- txtProgressBar(0, J, style=3)

    # define midpoints for overlapping bins
    s_m  <- 1:J

    for(j in s_m){
      sind_j <- (j-binwidth/2):(j+binwidth/2) %% J
      sind_j[sind_j == 0] <- J
      df_j <- filter(Y, index %in% argvals[sind_j])
      fit_j <- glmer(value ~ 1 + (1|id), data=df_j, family=family, nAGQ = 0)
      fit_fastgfpca[[j]] <- data.frame("id" = 1:N,
                                         "eta_i" = coef(fit_j)$id[[1]],
                                         "sind_bin" = j)
      #setTxtProgressBar(pb, j)
    }
    fit_fastgfpca <- bind_rows(fit_fastgfpca)


    if(J/binwidth < 20){
      knots <- ceiling(J/binwidth/2)
    }else{
      knots <- 20
    }
    knots <- 20

    fastgfpca <- fpca.face(matrix(fit_fastgfpca$eta_i, N, J, byrow=FALSE),
                           npc=npc,
                           pve=pve,
                           argvals = argvals,
                           knots=knots,lower=0,
                           periodicity = periodicity)

    if(is.null(npc)){
      npc = fastgfpca$npc
    }

  }else{


    # define midpoints for non overlapping bins
    s_m <- seq(1, J, by=binwidth+1)
    s_m[1] <- median(1:(binwidth/2 + 1))
    last_bin <- (s_m[length(s_m)]-binwidth/2):J
    s_m[length(s_m)] <- median(last_bin)

    # create bins, which will have the following widths
    # 1st bin: always binwidth/2 + 1
    # last bin: between binwidth/2 + 1 and binwidth + 1 + (binwidth/2)
    # other bins: binwidth + 1
    bins <- rep(s_m[-length(s_m)], each = binwidth + 1)
    bins <- bins[-c(1:(binwidth/2))]
    bins <- c(bins, rep(s_m[length(s_m)], length.out = length(last_bin)))

    df_bin <- mutate(Y, sind_bin = rep(bins, N))

    # fit local model
    fit_fastgfpca <- df_bin %>%
      nest_by(sind_bin) %>%
      mutate(fit = list(glmer(value ~ 1 + (1|id), data=data, family=family, nAGQ = 0))) %>%
      summarize(eta_i = coef(fit)$id[[1]],
                b_i = eta_i - fit@beta,
                id = 1:N) %>%
      ungroup()

    # do FPCA on the local estimates \tilde{\eta_i(s)}
    # knots will be given by bindwidth for smaller values of D
    # need to edit number of knots here
    if(J/binwidth <= 40 + 10){
      knots <- ceiling(J/binwidth/2)
    }else{
      knots <- 40
    }

    argvals_bin <- sort(argvals[unique(bins)])
    fastgfpca <- fpca.face(matrix(fit_fastgfpca$eta_i, N, length(argvals_bin), byrow=FALSE),
                           npc = npc,
                           pve=pve,
                           argvals=argvals_bin,
                           knots=knots,
                           periodicity = periodicity)

    if(is.null(npc)){
      npc = fastgfpca$npc
    }

    fastgfpca$efunctions <- reeval_efunctions(knots, argvals_bin, argvals,
                                              fastgfpca$efunctions, npc)
  }





  Y$id_fac <- factor(Y$id)
  gam_formula = "value ~ s(index, k=10)"
  for(i in 1:npc){
    Y[[paste0("Phi", i)]] <- rep(fastgfpca$efunctions[,i], N)
    gam_formula = paste0(gam_formula, " + s(id_fac, by=Phi",i,", bs= 're')")
  }

  # fit model using eigenfunctions as covariates to update scores
  fit_fastgfpca <- bam(formula = as.formula(gam_formula),
                       method="fREML", data=Y, family=family, discrete=TRUE,
                       ...
  )

  # next:return proper mu and scores
  eta_hat <- predict(fit_fastgfpca, newdata=Y, type='link')
  score_hat <- coef(fit_fastgfpca)
  fastgfpca$mu <- (predict(fit_fastgfpca, type = "terms")[,1] + score_hat[1])[1:J]

  fastgfpca$Yhat <- matrix(eta_hat, N, J, byrow = TRUE)
  fastgfpca$family <- fit_fastgfpca$family
  fastgfpca$Y <- NULL # do not store Y

  fastgfpca$scores <- matrix(score_hat[grep("Phi",names(score_hat))], N, npc)


  fastgfpca

}
