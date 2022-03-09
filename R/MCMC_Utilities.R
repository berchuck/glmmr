###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples <- function(DatObj, OmegaSamples) {

  ###Set data objects
  P <- DatObj$P
  Q <- DatObj$Q
  NL <- DatObj$NL
  
  ###Format raw samples
  OmegaSamples <- t(OmegaSamples)
  NSims <- nrow(OmegaSamples)
  Beta <- OmegaSamples[, 1:P, drop = FALSE]
  l <- OmegaSamples[, (P + NL):(P + NL), drop = FALSE]
  d <- OmegaSamples[, (P + NL + 1):(P + NL + Q), drop = FALSE]
  colnames(Beta) <- paste0("Beta", 0:(P - 1))
  colnames(l) <- paste0("l", 1:NL)
  colnames(d) <- paste0("d", 1:Q)
  Upsilon <- matrix(nrow = NSims, ncol = NL)
  Sigma <- matrix(nrow = NSims, ncol = NL + Q)
  for (s in 1:NSims) {
    L <- GetL(GetZ(l[s, ], Q), Q)
    UpsilonMat <- L %*% t(L)
    D <- diag(exp(d[s, ]))
    SigmaMat <- D %*% UpsilonMat %*% D
    Upsilon[s, ] <- UpsilonMat[lower.tri(UpsilonMat, diag = FALSE)]
    Sigma[s, ] <- SigmaMat[lower.tri(SigmaMat, diag = TRUE)]
  }
  SigmaInd <- which(lower.tri(apply(matrix(1:Q, ncol = 1), 1, function(x) paste0(paste0("Sigma", 1:Q), x)), diag = TRUE), arr.ind = TRUE)
  if (Q == 1) colnames(Sigma) <- apply(matrix(1:Q, ncol = 1), 1, function(x) paste0(paste0("Sigma", 1:Q), x))[SigmaInd[order(SigmaInd[, 1]), ]][1]
  if (Q > 1) colnames(Sigma) <- apply(matrix(1:Q, ncol = 1), 1, function(x) paste0(paste0("Sigma", 1:Q), x))[SigmaInd[order(SigmaInd[, 1]), ]]
  
  UpsilonInd <- which(lower.tri(apply(matrix(1:Q, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:Q, "_"), x)), diag = FALSE), arr.ind = TRUE)
  colnames(Upsilon) <- paste0("Upsilon", paste0(UpsilonInd, collapse = ""))
  Out <- list(Beta = Beta, l = l, d = d, Upsilon = Upsilon, Sigma = Sigma)
  return(Out)
}



###Function for creating a data object that contains objects needed for ModelFit-----------------------------------
OutputDatObj <- function(DatObj) {

  ###Collect needed objects
  # DatObjOut <- list(M = DatObj$M,
  #                   Nu = DatObj$Nu,
  #                   AdjacentEdgesBoolean = DatObj$AdjacentEdgesBoolean,
  #                   W = DatObj$W,
  #                   EyeM = DatObj$EyeM,
  #                   EyeNu = DatObj$EyeNu,
  #                   OneM = DatObj$OneM,
  #                   OneN = DatObj$OneN,
  #                   OneNu = DatObj$OneNu,
  #                   YStarWide = DatObj$YStarWide,
  #                   Rho = DatObj$Rho,
  #                   FamilyInd = DatObj$FamilyInd,
  #                   ScaleY = DatObj$ScaleY,
  #                   YObserved = DatObj$YObserved,
  #                   ScaleDM = DatObj$ScaleDM,
  #                   Time = DatObj$Time,
  #                   TimeVec = DatObj$TimeVec,
  #                   YObserved = DatObj$YObserved,
  #                   tNu = DatObj$tNu,
  #                   t1 = DatObj$t1,
  #                   XThetaInd = DatObj$XThetaInd,
  #                   N = DatObj$N,
  #                   EyeN = DatObj$EyeN)
  DatObjOut <- DatObj
  return(DatObjOut)

}



###Verify the class of our regression object------------------------------------------------------------------------
#' is.glmmr
#'
#' \code{#' is.glmmr} is a general test of an object being interpretable as a
#' \code{\link{glmmr}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{glmmr}} class is defined as the regression object that
#'  results from the \code{\link{glmmr}} regression function.
#'
#' @return \code{is.glmmr} returns a logical, depending on whether the object is of class \code{\link{glmmr}}.
#'   
#' @examples 
#' ###Load pre-computed results
#' data(reg_glmmr)
#' 
#' ###Test function
#' is.glmmr(reg_glmmr)
#'
#' @export
is.glmmr <- function(x) {
  identical(attributes(x)$class, "glmmr")
}



###GetL
#' GetL
#'
#' \code{GetL} is a function to computes the Cholesky factor for a correlation matrix from the intermediate matrix Z.
#'
#' @param Z Z matrix
#' @param Q dimension of the correlation matrix
#'
#' @return \code{GetL} returns a matrix of size Q
#'   
#' @examples 
#' ###Test function
#' l <- -0.34
#' Z <- GetZ(l, 2)
#' L <- GetL(Z, 2)
#'
#' @export
GetL <- function(Z, Q) {
  GetL(Z, Q)
}



###GetZ
#' GetZ
#'
#' \code{GetZ} is a function to compute the intermediate step for the Cholesky factor of a correlation matrix.
#'
#' @param l l vector of unconstrained parameters
#' @param Q dimension of the correlation matrix
#'
#' @return \code{GetZ} returns a matrix of size Q
#'   
#' @examples 
#' ###Test function
#' l <- -0.34
#' Z <- GetZ(l, 2)
#' L <- GetL(Z, 2)
#'
#' @export
GetZ <- function(l, Q) {
  GetZ(l, Q)
}