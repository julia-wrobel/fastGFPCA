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
#' @importFrom stats approx
#' @importFrom refund fpca.face
#' @importFrom lme4 glmer
#' @importFrom utils txtProgressBar
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
#' df_gfpca <- sim_gfpca(N = 200, J = 200, case = 1)$df_gfpca
#' gfpca_mod <- fast_gfpca(df_gfpca, overlap = TRUE, binwidth = 10, family = "binomial")
#'
#'@export

fast_gfpca <- function(df, overlap = TRUE, binwidth = 10,
                       npc = 4, # need to build in npc argument
                       family = "binomial",
                       extrapolate = FALSE,
                       ...){
  # create indicator function for bin
  N <- length(unique(df$id))
  J <- length(unique(df$index)) # assumes all subjects are on same even grid
  sind <- (1:J)/J

  if(overlap){
    fit_fastgfpca <- vector(mode="list",length=J)
    pb <- txtProgressBar(0, J, style=3)
    for(j in 1:J){
      sind_j <- (j-binwidth/2):(j+binwidth/2) %% J + 1
      df_j <-df %>%
        filter(index %in% sind[sind_j])
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
                             npc=npc, argvals=sind, knots=knots,lower=0)

  }else{
    # create indicator for bin
    df_bin <-df %>%
      mutate(sind_bin = rep(floor((1:J)/binwidth)*binwidth + binwidth/2, N)) %>%
      filter(sind_bin < J)

    # fit local model
    fit_fastgfpca <- df_bin %>%
      nest_by(sind_bin) %>%
      mutate(fit = list(glmer(value ~ 1 + (1|id), data=data, family=family))) %>%
      summarize(eta_i = coef(fit)$id[[1]],
                b_i = eta_i - fit@beta,
                id = 1:N) %>%
      ungroup()

    if(extrapolate){
      fit2 = df_bin %>%
        filter(sind_bin %in% c(max(sind_bin), min(sind_bin))) %>%
        mutate(index_int = index * J) %>%
        nest_by(sind_bin) %>%
        mutate(fit = list(glmer(value ~ index + (1|id), data=data, family=family,
                                ...))) %>%
        summarize(eta_i = predict(fit, type = "link"),
                  id = data$id,
                  sind_bin = data$index * J) %>%
        ungroup()

      fit_fastgfpca <- full_join(fit_fastgfpca, (fit2 %>%
                                                   filter(sind_bin < min(fit_fastgfpca$sind_bin) | sind_bin > max(fit_fastgfpca$sind_bin))))
    }

    # do FPCA on the local estimates \tilde{\eta_i(s)}
    # knots will be given by bindwidth for smaller values of D
      # need to edit number of knots here
    if(J/binwidth < 40){
      knots <- ceiling(J/binwidth/2)
    }else{
      knots <- 40
    }
    sind_bin  <- sort(unique(df_bin$sind_bin))
    fastgfpca <- fpca.face(matrix(fit_fastgfpca$eta_i, N, J/binwidth, byrow=FALSE),
                           pve=0.99, argvals=sind_bin, knots=knots)

    ## approximate eigenfunctions via linear interpolation
    ## on the original grid
    ## move this to utils
    efuncs_approx <- matrix(NA, J, 4)
    for(k in 1:4){
      efuncs_approx[,k] <- approx(sind_bin, fastgfpca$efunctions[,k], xout=1:J)$y
    }
    ## re-scale eigenfunctions to have the correct magnitude
    efuncs_approx <- efuncs_approx*sqrt(J/binwidth)
    fastgfpca$efunctions <- efuncs_approx
  }

  # generalize this code to more than 4 eigenfunctions
  df$id_fac <- factor(df$id)
  df$Phi1 <- fastgfpca$efunctions[,1]
  df$Phi2 <- fastgfpca$efunctions[,2]
  df$Phi3 <- fastgfpca$efunctions[,3]
  df$Phi4 <- fastgfpca$efunctions[,4]

  # fit model using eigenfunctions as covariates to update scores
  fit_fastgfpca <- bam(value ~ s(index, k=10) +
                              s(id_fac, by=Phi1, bs="re") + s(id_fac, by=Phi2, bs="re") +
                              s(id_fac, by=Phi3, bs="re") + s(id_fac, by=Phi4, bs="re"),
                            method="fREML", data=df, family=family, discrete=TRUE,
                       ...)



  eta_hat <- predict(fit_fastgfpca, newdata=df, type='link')
  # return something like a standard refund FPCA list object
  # return Yhat matrix (really eta_hat) where each row is a subject and each column a point in time

  # next:return proper mu and scores
  fastgfpca$Yhat <- matrix(eta_hat, N, J, byrow = TRUE)

  list(pca = fastgfpca, gam_mod = fit_fastgfpca)

}
