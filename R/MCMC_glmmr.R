#' Scalable Bayesian generalized linear mixed models using stochastic gradient Markov chain Monte Carlo
#'
#' \code{glmmr} is a stochastic gradient Markov chain Monte Carlo (MCMC) sampler for a Bayesian generalized linear mixed model.
#'
#' @param pformula A \code{formula} object, corresponding to the population level model. The response must be on the left of a \code{~} operator, and the terms on the right 
#'                 must indicate the covariates to be included in the fixed effects.

#' @param gformula A \code{formula} object, corresponding to the group level model. There is no response variable to the left of the \code{~} operator. The terms on the right 
#'                 must indicate the covariates to be included in the mixed effects.
#'
#' @param group A character string, the name of the variable that defines the group levels.
#'                 
#' @param data A required \code{data.frame} containing the variables in the model. The data frame should be sorted by the \code{group} variable.
#'
#' @param family Character string indicating the distribution of the observed data. Options
#'  include: \code{"normal"}, \code{"binomial"}, and \code{"poisson"}.
#'  
#' @param algorithm Character string indicating the algorithm to be used. Options
#'  include: \code{"sgd"}, \code{"sgld"}, \code{"sgld_corrected"}, and \code{"gibbs"}.
#'  
#' @param starting Either \code{NULL} or a \code{list} containing starting values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the starting values may be specified.
#'
#'  When \code{NULL} is chosen then default starting values are automatically generated.
#'  Otherwise a \code{list} must be provided with names \code{Beta}, \code{l}, or \code{d}, containing appropriate objects. 
#'  \code{Beta} must either be a \code{P} dimensional vector or a scalar (the scalar populates the entire vector), where \code{P} is the number of population level regression parameters.
#'  \code{l} must either be a \code{NL} dimensional vector or a scalar (the scalar populates the entire vector), where \code{NL} is the number of unconstrained parameters in the Cholesky factor for a correlation matrix.
#'  \code{d} must either be a \code{Q} dimensional vector or a scalar (the scalar populates the entire vector), where \code{Q} is the number of group level regression parameters.
#'
#' @param hypers Either \code{NULL} or a \code{list} containing hyperparameter values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the hyperparameter values may be specified.
#'
#'  When \code{NULL} is chosen then default hyperparameter values are automatically
#'  generated. These default hyperparameters are described in detail in \code{brms}.
#'  Otherwise a \code{list} must be provided with names \code{l} or \code{d} containing further hyperparameter information. These objects are themselves
#'  \code{lists} and may be constructed as follows.
#'
#'  \code{l} is a \code{list} with one object, \code{Eta}. This value represents the prior degree of freedom 
#'  parameter for the LKJ prior.
#'  
#'  \code{d} is a \code{list} with one object, \code{Nu}. This value represents the prior degree of freedom parameter for a half-t distribution.
#'  
#' @param tuning Either \code{NULL} or a \code{list} containing tuning values
#'  to be specified for the optimization and sampling routines. If \code{NULL} is not chosen then none, some or all
#'  of the tuning values may be specified.
#'
#'  \code{EpsilonNADAM}: Epsilon value in the NADAM algorithm. (default = 1e-8)
#'  
#'  \code{MuNADAM}: Mu value in the NADAM algorithm. (default = 0.975)
#'  
#'  \code{AlphaNADAM}: Alpha value in the NADAM algorithm. (default = 0.002)
#'  
#'  \code{NuNADAM}: Nu value in the NADAM algorithm. (default = 0.999)
#'  
#'  \code{S}: The size of the mini-batches. (default = 10 percent of the groups or 10, whichever is smaller)
#'  
#'  \code{NEpochs}: The number of epochs in the NADAM phase. If using the \code{gibbs} algorithm, this number is the size of the burn-in samples. (default =
#'  \code{1,000})
#'  
#'  \code{R}: The number of MCMC samples to draw to sample the random effects. (default = 500)
#'
#'  \code{EpsilonSGLD}: Epsilon value in the SGLD algorithm. (default = 0.01)
#'  
#'  \code{EpsilonSGLDCorrected}: Positive number that is used to divide the epsilon value in the SGLD corrected algorithm. (default = 100)
#'  
#'  \code{S_SGLD}: The number of groups that are randomly sampled to compute the SGLD correction. (default = total number of groups)
#'  
#'  \code{NSims}: The number of SGLD scans for which to perform the
#'   sampler. (default = \code{5,000})
#'
#'  \code{NThin}: Value such that during the SGLD phase, only every
#'  \code{NThin}-th scan is recorded for use in posterior inference (For return values
#'  we define, NKeep = NSims / NThin (default = \code{1}).
#'  
#'  \code{NPilot}: If using the \code{gibbs} algorithm, the number of times during the burn-in phase that pilot adaptation
#'  is performed. (default = \code{20})
#'  
#'  \code{TuneL}: If using the \code{gibbs} algorithm, the tuning variances for \code{l} in the Metropolis step. Either a scalar or an \code{NL} dimensional vector. (default = 1)
#'  
#'  \code{TuneD}: If using the \code{gibbs} algorithm, the tuning variances for \code{d} in the Metropolis step. Either a scalar or an \code{Q} dimensional vector. (default = 1)
#'  
#'  \code{NTune}: The number of iterations to tune the \code{EpsilonSGLD} value to avoid singularity of \code{L} (default = \code{100}). 
#'  This is only performed for \code{algorithm = "sgld"} or \code{"sgld_corrected"}. For \code{"sgld"}, the initial value of epsilon is 
#'  determined using \code{EpsilonSGLD}. Then, at each iteration, if a value of \code{l} is not able to produce a non-singular value of \code{L} 
#'  within \code{NTune_seconds}, \code{EpsilonSGLD} is divided by 10. This repeats for \code{NTune} iterations. For \code{"sgld_corrected"},
#'  the intial value of epsilon is determined as the minimum of \code{EpsilongSGLD} and the maximum eigen value of the correction matrix.
#'  
#'  \code{NTune_seconds}: The number of seconds that a non-singular value of \text{l} are sampled for for each interation of \code{NTune} (default = 1).
#'  
#' @param seed An integer value used to set the seed for the random number generator
#'  (default = 54).
#'  
#' @param timer An integer value indicating the number of seconds the algorithm should be run (default = NULL). If a value is specificed the algorithm will disregard the \code{NSims} and run until the value is reached.
#'  
#' @details Details of the underlying statistical model proposed by
#'  Berchuck et al. 2019. are forthcoming.
#'
#' @return \code{glmmr} returns a list containing the following objects. Note that 
#' if the sgd algorithm was used the output are not samples, but represent the path to the MAP.
#'
#'   \describe{
#'
#'   \item{\code{beta}}{\code{NKeep x P} \code{matrix} of posterior samples for \code{beta}.}
#'
#'   \item{\code{l}}{\code{NKeep x NL} \code{matrix} of posterior samples for the unconstrained parameters of the Cholesky factor for a correlation matrix.}
#'
#'   \item{\code{d}}{\code{NKeep x Q} \code{matrix} of posterior samples for the standard deviations for the random effects.}
#'
#'   \item{\code{upsilon}}{\code{NKeep x NL} \code{matrix} of posterior samples for the correlation matrix \code{upsilon}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{Upsilon32} refers to the entry in row \code{3}, column \code{2}.}
#'   
#'   \item{\code{sigma}}{\code{NKeep x NL + Q} \code{matrix} of posterior samples for the covariance matrix \code{sigma}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{sigma32} refers to the entry in row \code{3}, column \code{2}.}
#'
#'   \item{\code{map}}{A \code{vector} containing the maximum a posteriori estimate.}
#'
#'   \item{\code{runtime}}{A \code{character} string giving the runtime of the MCMC sampler.}
#'
#'   \item{\code{datobj}}{A \code{list} of data objects that are used in future \code{bfa_sp} functions
#'   and should be ignored by the user.}
#'
#'   }
#'
#' @examples
#' \donttest{
#' ###Load simulated data
#' data(dat_sim)
#' 
#' ###Format data for glmmr
#' dat <- data.frame(y = dat_sim$y.1, time = dat_sim$time, id = dat_sim$id)
#' reg <- glmmr(y ~ 1 + time, ~ 1 + time, group = "id", data = dat, family = "bernoulli", 
#'              tuning = list(NSims = 1000, NEpochs = 1000, S = 1, NThin = 10))
#' 
#' }
#'
# @author Samuel I. Berchuck
#' @references Reference for Berchuck et al. 2022 is forthcoming.
#' @export
glmmr <- function(pformula, gformula, group, data, family = "binomial", algorithm = "sgld_corrected",
                  starting = NULL, hypers = NULL, tuning = NULL, seed = 54, timer = NULL) {
  
  ###Function Inputs
  # pformula = y ~ 1 + time
  # gformula = ~ 1 + time
  # group = "id"
  # data = dat
  # family = "bernoulli"
  # algorithm = "sgld_corrected"
  # starting = NULL
  # hypers = NULL
  # tuning = list(NSims = 1000, NEpochs = 1000, S = 1, NTune = 100, NTune_seconds = 1)
  # seed = 54
  # timer = NULL

  ###Check for missing objects
  if (missing(pformula)) stop("pformula: missing")
  if (missing(gformula)) stop("gformula: missing")
  if (missing(group)) stop("group: missing")
  if (missing(data)) stop("data: missing")

  ###Check model inputs
  CheckInputs(pformula, gformula, group, data, family, algorithm, starting, hypers, tuning, seed, timer)

  ####Set seed for reproducibility
  set.seed(seed)

  ###Check to see if the job is interactive
  Interactive <- interactive()

  ###Create objects for use in sampler
  DatObj <- CreateDatObj(pformula, gformula, group, data, family, algorithm, timer)
  HyPara <- CreateHyPara(hypers) 
  TuningObj <- CreateTuningObj(tuning, DatObj)
  Para <- CreatePara(starting, DatObj)

  ###Time MCMC sampler
  BeginTime <- Sys.time()

  ###Run SGD in Rcpp to obtain MAP
  RegObj <- glmmr_Rcpp(DatObj, HyPara, TuningObj, Para, Interactive)

  ###End time
  FinishTime <- Sys.time()
  RunTime <- FinishTime - BeginTime

  ###Set regression objects
  OmegaMAP <- RegObj$map
  Omega <- RegObj$samples
  TuningObj$EpsilonSGLDFinal <- RegObj$epsilon
  TuningObj$Counter <- RegObj$counter
  TuningObj$EpsilonTuneCounter <- RegObj$epsilontunecounter
  TuningObj$LambdaMax <- TuningObj$EpsilonSGLDFinal * 10^TuningObj$EpsilonTuneCounter
  
  ###Set metropolis objec
  MetropRcpp <- RegObj$metropolis
  
  ###Timer object
  Timer <- RegObj$timer
  LastIndex <- which.max(Timer == 0) - 1
  TimerOut <- NULL
  if (!is.null(timer)) {
    if (algorithm == "gibbs") {
      TimerOut <- data.frame(stage = c(rep("burnin", TuningObj$NEpochs), rep("sampling", TuningObj$NSims))[1:LastIndex],
                          iteration = 1:LastIndex,
                          seconds = Timer[1:LastIndex])
    }
    if (algorithm %in% c("sgld", "sgld_corrected")) {
      TimerOut <- data.frame(stage = c(rep("sgd", TuningObj$NEpochs), rep("tuning", TuningObj$NTune), rep("sampling", TuningObj$NSims))[1:LastIndex],
                             iteration = 1:LastIndex,
                             seconds = Timer[1:LastIndex])
    }
    if (algorithm == "sgd") {
      TimerOut <- data.frame(stage = c(rep("sgd", TuningObj$NEpochs))[1:LastIndex],
                          iteration = 1:LastIndex,
                          seconds = Timer[1:LastIndex])
    }
  }
  
  ###Collect output to be returned
  DatObjOut <- OutputDatObj(DatObj)
  if (!is.null(timer)) Omega <- Omega[, 1:LastIndex]
  Samples <- FormatSamples(DatObj, Omega)
  Metropolis <- NULL
  if (algorithm == "gibbs") {
    Metropolis <- SummarizeMetropolis(DatObj, TuningObj, MetropRcpp)
    if (!is.null(timer)) Metropolis[, 1] <- (Metropolis[, 1] * TuningObj$NSims) / (LastIndex - TuningObj$NEpochs)
  }
  
  ###Return spBFA object
  glmmr <- list(beta = Samples$Beta,
                l = Samples$l,
                d = Samples$d,
                upsilon = Samples$Upsilon,
                sigma = Samples$Sigma,
                map = OmegaMAP,
                datobj = DatObjOut,
                tuningobj = TuningObj,
                hypara = HyPara,
                starting = Para$Omega,
                metropolis = Metropolis,
                runtime = RunTime,
                timer = TimerOut)
  glmmr <- structure(glmmr, class = "glmmr")
  return(glmmr)

###End sampler
}
