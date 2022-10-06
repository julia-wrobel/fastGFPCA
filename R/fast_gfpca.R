#' Fast generalized functional principal components analysis
#'
#' I fast implementation of GFPCA that uses the ...
#' The number of functional principal components (FPCs) can either be specified
#' directly (argument \code{npc}) or chosen based on the explained share of
#' variance (\code{npc_varExplained}). In the latter case, the explained share of
#' variance and accordingly the number of FPCs is estimated before the main
#' estimation step by once running the FPCA with \code{npc = 20} (and
#' correspondingly \code{Kt = 20}). Doing so, we approximate the overall
#' variance in the data \code{Y} with the variance represented by the FPC basis
#' with 20 FPCs.
#'

#'
#' @author Andrew Leroux \email{andrew.leroux@@cuanschutz.edu},
#' Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' @import dplyr
#' @importFrom stats approx binomial coef predict binomial
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
#' df_gfpca <- sim_gfpca(N = 50, J = 100, case = 1)$df_gfpca
#' gfpca_mod <- fast_gfpca(df_gfpca, overlap = TRUE, binwidth = 10, family = "binomial")
#'
#' @param Y dataframe with very specific column
#' @param argvals numeric; grid over which functions are observed.  If null defaults to unique values of index.
#' @param pve proportion of variance explained: used to choose the number of
#' principal components unless `npc` is specified.
#' @param npc how many smooth PCs to try to extract, if \code{NULL} (the
#' default) then npc is chosen based on `pve`
#' @param family exponential family to be passed to \code{glmer} and \code{bam}.
#' @param periodicity Option for a periodic spline basis. Defaults to FALSE.
#' .
#'@export

fast_gfpca <- function(Y,
                       argvals = NULL, # grid for functional observations
                       overlap = TRUE,
                       binwidth = 10,
                       pve = 0.99,
                       npc = NULL,
                       family = "binomial",
                       periodicity = FALSE,
                       ...){

  # add check that binwidth is even. If not, it will bet converted to an even number


  N <- length(unique(Y$id))
  J <- length(unique(Y$index)) # assumes all subjects are on same even grid

  if(is.null(argvals)){
    argvals = sort(unique(Y$index))
  }

  if(overlap){
    fit_fastgfpca <- vector(mode="list",length=J)
    pb <- txtProgressBar(0, J, style=3)
    for(j in 1:J){
      sind_j <- (j-binwidth/2):(j+binwidth/2) %% J + 1
      df_j <-Y %>%
        filter(index %in% argvals[sind_j])
        fit_j <- glmer(value ~ 1 + (1|id), data=df_j, family=binomial)
        fit_fastgfpca[[j]] <- data.frame("id" = 1:N,
                                         "eta_i" = coef(fit_j)$id[[1]],
                                         "sind_bin" = j)
      setTxtProgressBar(pb, j)
    }
    fit_fastgfpca <- bind_rows(fit_fastgfpca)

    if(J/binwidth < 20){
      knots <- ceiling(J/binwidth/2)
    }else{
      knots <- 20
    }

    fastgfpca <- fpca.face(matrix(fit_fastgfpca$eta_i, N, J, byrow=FALSE),
                           npc=npc, pve=0.99,
                           argvals = argvals, knots=knots,lower=0,
                           periodicity = periodicity)

    if(is.null(npc)){
      npc = fastgfpca$npc
    }

  }else{
    # create indicator for asymmetric bins
    bins = c(rep(1, ceiling(binwidth/2)), rep(seq(binwidth, (J-binwidth), by = binwidth),
                                              each = binwidth),
             rep(J, ceiling(binwidth/2)))

    df_bin <- Y %>% mutate(sind_bin = rep(bins, N))

    # fit local model
    fit_fastgfpca <- df_bin %>%
      nest_by(sind_bin) %>%
      mutate(fit = list(glmer(value ~ 1 + (1|id), data=data, family=family))) %>%
      summarize(eta_i = coef(fit)$id[[1]],
                b_i = eta_i - fit@beta,
                id = 1:N) %>%
      ungroup()

    # do FPCA on the local estimates \tilde{\eta_i(s)}
    # knots will be given by bindwidth for smaller values of D
    # need to edit number of knots here
    if(J/binwidth < 40){
      knots <- ceiling(J/binwidth/2)
    }else{
      knots <- 40
    }

    sind_bin  <- sort(unique(fit_fastgfpca$sind_bin))
    fastgfpca <- fpca.face(matrix(fit_fastgfpca$eta_i, N, length(sind_bin), byrow=FALSE),
                           npc = npc, pve=0.99,
                           argvals=sind_bin, knots=knots,
                           periodicity = periodicity)

    if(is.null(npc)){
      npc = fastgfpca$npc
    }

    fastgfpca$efunctions <- reeval_efunctions(knots, sind_bin, argvals, fastgfpca$efunctions, npc)
  }

  ## re-scale eigenfunctions to have the correct magnitude
  fastgfpca$efunctions <- fastgfpca$efunctions*sqrt(J)


  Y$id_fac <- factor(Y$id)
  gam_formula = "value ~ s(index, k=10)"
  for(i in 1:npc){
    Y[[paste0("Phi", i)]] <- fastgfpca$efunctions[,i]
    gam_formula = paste0(gam_formula, " + s(id_fac, by=Phi",i,", bs= 're')")
  }


  # fit model using eigenfunctions as covariates to update scores
  fit_fastgfpca <- bam(formula = as.formula(gam_formula),
                       method="fREML", data=Y, family=family, discrete=TRUE,
                       ...
  )

  # next:return proper mu and scores
  # for now return scores from face and bam
  eta_hat <- predict(fit_fastgfpca, newdata=Y, type='link')
  score_hat <- coef(fit_fastgfpca)
  fastgfpca$mu <- (predict(fit_fastgfpca, type = "terms")[,1] + score_hat[1])[1:J]

  fastgfpca$Yhat <- matrix(eta_hat, N, J, byrow = TRUE)
  fastgfpca$family <- family

  fastgfpca$scores_face <- fastgfpca$scores
  fastgfpca$scores <- matrix(score_hat[grep("Phi",names(score_hat))], N, npc)

  #fastgfpca$fit <- fit_fastgfpca

  fastgfpca

}
