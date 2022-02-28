#' @export
glmmr <- function(pformula, gformula, group, data, family = "binomial", 
                  starting = NULL, hypers = NULL, sgd = NULL, tuning = NULL, mcmc = NULL, seed = 54) {
  
  # ###Function Inputs
  # pformula = y ~ 1 + time
  # gformula = ~ 1 + time
  # group = "id"
  # data = dat
  # family = "binomial"
  # starting = NULL
  # hypers = NULL
  # sgd = NULL
  # tuning = NULL
  # mcmc = NULL
  # seed = 54
  
  ###Check for missing objects
  if (missing(pformula)) stop("pformula: missing")
  if (missing(gformula)) stop("gformula: missing")
  if (missing(group)) stop("group: missing")
  if (missing(data)) stop("data: missing")

  ###Check model inputs
  # CheckInputs(formula, data, family, starting, hypers, tuning, mcmc, seed)

  ####Set seed for reproducibility
  set.seed(seed)

  ###Check to see if the job is interactive
  Interactive <- interactive()

  ###Create objects for use in sampler
  DatObj <- CreateDatObj(pformula, gformula, group, data, family)
  HyPara <- CreateHyPara(hypers, DatObj) 
  SgdObj <- CreateSgdObj(sgd, DatObj)
  Para <- CreatePara(starting, DatObj, HyPara)

  # ###Time MCMC sampler
  # BeginTime <- Sys.time()

  ###Run SGD in Rcpp
  Omega <- glmmr_sgd_Rcpp(DatObj, HyPara, SgdObj, Para, Interactive)

  # ###Run SGMCMC sampler in Rcpp
  # RegObj <- glmmr_sgmcmc_Rcpp(DatObj, HyPara, MetrObj, Para, DatAug, McmcObj, RawSamples, Interactive)
  # 
  # ###Set regression objects
  # RawSamples <- RegObj$rawsamples
  # MetropRcpp <- RegObj$metropolis
  # 
  # ###End time
  # FinishTime <- Sys.time()
  # RunTime <- FinishTime - BeginTime
  # 
  # ###Collect output to be returned
  # DatObjOut <- OutputDatObj(DatObj)
  # DatAugOut <- OutputDatAug(DatAug)
  # Metropolis <- SummarizeMetropolis(DatObj, MetrObj, MetropRcpp, McmcObj)
  # Samples <- FormatSamples(DatObj, RawSamples)
  # 
  # ###Return spBFA object
  # spBFA <- list(lambda = Samples$Lambda,
  #               eta = Samples$Eta,
  #               beta = Samples$Beta,
  #               sigma2 = Samples$Sigma2,
  #               kappa = Samples$Kappa,
  #               delta = Samples$Delta,
  #               tau = Samples$Tau,
  #               upsilon = Samples$Upsilon,
  #               psi = Samples$Psi,
  #               xi = Samples$Xi,
  #               rho = Samples$Rho,
  #               metropolis = Metropolis,
  #               datobj = DatObjOut,
  #               dataug = DatAugOut,
  #               runtime = paste0("Model runtime: ",round(RunTime, 2), " ", attr(RunTime, "units")))
  # spBFA <- structure(spBFA, class = "spBFA")
  return(Omega)

###End sampler
}
