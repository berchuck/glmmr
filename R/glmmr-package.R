#' glmmr
#'
#' \code{glmmr} is a package for Bayesian generalized linear mixed models using stochastic gradient Markov chain Monte Carlo.
#' 
#' @author Samuel I. Berchuck \email{sib2@duke.edu}
#'
#' @name glmmr
#' @docType package
#' @import Rcpp
#' @importFrom stats quantile rnorm runif median model.frame model.matrix
#' @importFrom graphics abline axis layout par plot points polygon title segments symbols rect text lines
#' @importFrom grDevices col2rgb colorRampPalette
#' @importFrom utils tail
#' @importFrom stats dnorm lm sd var pnorm
#' @importFrom msm rtnorm
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm pmvnorm
#' @importFrom pgdraw pgdraw
#' @useDynLib glmmr
NULL
