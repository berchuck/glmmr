###Function for reading in sampler inputs and creating a list object that contains all relevant data objects--------------------
CreateDatObj <- function(pformula, gformula, group, data, family, algorithm) {

  ###Extract group and convert to ordered integers
  dat_group <- dat[, group]
  dat_group <- as.numeric(as.factor(dat_group))
  
  ###Data objects
  N <- nrow(data) # total observations
  NUnits <- length(unique(dat_group)) #number of spatial locations
  data <- data[order(dat_group), ]
  Y <- matrix(data[, all.vars(pformula)[1]], ncol = 1)
  Group <- dat_group

  ###Covariates
  X <- model.matrix(pformula, data = data)
  ZMat <- model.matrix(gformula, data = data)
  P <- ncol(X)
  Q <- ncol(ZMat)
  X <- matrix(X, nrow = N, ncol = P)
  Z <- bdiag(lapply(split(ZMat, Group), function(x) matrix(x, ncol = Q)))
  Z <- matrix(Z, nrow = N, ncol = Q * NUnits)
  Group2 <- rep(1:NUnits, each = Q)
  NL <- choose(Q, 2)
  NOmega <- P + NL + Q
  
  ###Matrix Objects
  EyeQ <- diag(Q)
  EyeNOmega <- diag(NOmega)
  
  ###Family indicator
  FamilyInd <- -1
  if (family == "normal") FamilyInd <- 0
  if (family == "bernoulli") FamilyInd <- 1
  if (family == "poisson") FamilyInd <- 2
  
  ###Algorithm indicator
  AlgorithmInd <- -1
  if (algorithm == "sgd") AlgorithmInd <- 0
  if (algorithm == "sgld") AlgorithmInd <- 1
  if (algorithm == "sgld_corrected") AlgorithmInd <- 2
  
  ###Make parameters global
  DatObj <- list()
  DatObj$Y <- Y
  DatObj$X <- X
  DatObj$Z <- Z
  DatObj$Group <- Group
  DatObj$Group2 <- Group2
  DatObj$N <- N
  DatObj$NOmega <- NOmega
  DatObj$NUnits <- NUnits
  DatObj$P <- P
  DatObj$Q <- Q
  DatObj$NL <- NL
  DatObj$FamilyInd <- FamilyInd
  DatObj$AlgorithmInd <- AlgorithmInd
  DatObj$EyeQ <- EyeQ
  DatObj$EyeNOmega <- EyeNOmega
  return(DatObj)

}



###Function to create Hyper-parameter Object------------------------------------------------------------------------------------
CreateHyPara <- function(hypers) {

  ###Which parameters are user defined?
  Userhypers <- names(hypers)

  ###Set hyper-parameters for L
  if ("L" %in% Userhypers) {
    Eta <- hypers$l$Eta
  }
  if (!("L" %in% Userhypers)) {
    Eta <- 1
  }
  
  ###Set hyper-parameters for D
  if ("D" %in% Userhypers) {
    Nu <- hypers$d$Nu
  }
  if (!("D" %in% Userhypers)) {
    Nu <- 3
  }

  ###Create object for hyper-parameters
  HyPara <- list()
  HyPara$Eta <- Eta
  HyPara$Nu <- Nu
  return(HyPara)

}



###Function for creating an object containing relevant tuning information---------------------------------------------------
CreateTuningObj <- function(tuning, DatObj) {
  
  ###Set data objects
  P <- DatObj$P
  Q <- DatObj$Q
  NL <- DatObj$NL
  NUnits <- DatObj$NUnits
  AlgorithmInd <- DatObj$AlgorithmInd
  
  ###Which parameters are user defined?
  UserTuning <- names(tuning)
  
  ###Optimization
  if ("EpsilonNADAM" %in% UserTuning) EpsilonNADAM <- tuning$EpsilonNADAM
  if (!("EpsilonNADAM" %in% UserTuning)) EpsilonNADAM <- 1e-8
  if ("MuNADAM" %in% UserTuning) MuNADAM <- tuning$MuNADAM
  if (!("MuNADAM" %in% UserTuning)) MuNADAM <- 0.975
  if ("AlphaNADAM" %in% UserTuning) AlphaNADAM <- tuning$AlphaNADAM
  if (!("AlphaNADAM" %in% UserTuning)) AlphaNADAM <- 0.002
  if ("NuNADAM" %in% UserTuning) NuNADAM <- tuning$NuNADAM
  if (!("NuNADAM" %in% UserTuning)) NuNADAM <- 0.999
  if ("S" %in% UserTuning) S <- tuning$S
  if (!("S" %in% UserTuning)) S <- min(10, NUnits)
  if ("NEpochs" %in% UserTuning) NEpochs <- tuning$NEpochs
  if (!("NEpochs" %in% UserTuning)) NEpochs <- 1000
  if ("R" %in% UserTuning) R <- tuning$R
  if (!("R" %in% UserTuning)) R <- 500
  
  ###Sampling
  if ("EpsilonSGLD" %in% UserTuning) EpsilonSGLD <- tuning$EpsilonSGLD
  if (!("EpsilonSGLD" %in% UserTuning)) EpsilonSGLD <- 0.01
  if ("NSims" %in% UserTuning) NSims <- tuning$NSims
  if (!("NSims" %in% UserTuning)) NSims <- 5000
  if ("NThin" %in% UserTuning) NThin <- tuning$NThin
  if (!("NThin" %in% UserTuning)) NThin <- 1
  
  ###One last check of MCMC user inputs
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!(is.wholenumber(NSims / NThin))) stop('tuning: "NThin" must be a factor of "NSims"')

  ###Burn-in progress bar
  WhichKeep <- NEpochs + (1:(NSims / NThin)) * NThin
  NKeep <- length(WhichKeep)
  BarLength <- 50 #Burn-in bar length (arbitrary)
  ProgressBar <- seq(1 / BarLength, 1, 1 / BarLength)
  SamplerProgress <- seq(0.1, 1.0, 0.1) #Intervals of progress update (arbitrary)
  WhichMAPProgress <- sapply(ProgressBar, function(x) tail(which(1 : NEpochs <= x * NEpochs), 1))
  WhichSGLDProgress <- sapply(ProgressBar, function(x) tail(which(1 : NUnits <= x * NUnits), 1))
  WhichSamplerProgress <- sapply(ProgressBar, function(x) tail(which(1 : NSims <= x * NSims), 1)) + NEpochs
  WhichMAPProgressInt <- sapply(SamplerProgress, function(x) tail(which(1:NEpochs <= x * NEpochs), 1))
  WhichSGLDProgressInt <- sapply(SamplerProgress, function(x) tail(which(1:NUnits <= x * NUnits), 1))
  WhichSamplerProgressInt <- sapply(SamplerProgress, function(x) tail(which(1:NSims <= x * NSims), 1)) + NEpochs
  if (AlgorithmInd == 0) NSims <- 0
  if (AlgorithmInd == 0) NKeep <- NEpochs
  
  ###Return SGD object
  TuningObj <- list()
  TuningObj$EpsilonNADAM <- EpsilonNADAM
  TuningObj$MuNADAM <- MuNADAM
  TuningObj$AlphaNADAM <- AlphaNADAM
  TuningObj$NuNADAM <- NuNADAM
  TuningObj$MNADAM <- matrix(0, nrow = P + NL + Q)
  TuningObj$NNADAM <- matrix(0, nrow = P + NL + Q)
  TuningObj$S <- S
  TuningObj$NEpochs <- NEpochs
  TuningObj$R <- R
  TuningObj$EpsilonSGLD <- EpsilonSGLD
  TuningObj$NSims <- NSims
  TuningObj$NThin <- NThin
  TuningObj$NTotal <- NSims + NEpochs
  TuningObj$WhichKeep <- WhichKeep
  TuningObj$NKeep <- NKeep
  TuningObj$BarLength <- BarLength
  TuningObj$WhichMAPProgress <- WhichMAPProgress
  TuningObj$WhichMAPProgressInt <- WhichMAPProgressInt
  TuningObj$WhichSamplerProgress <- WhichSamplerProgress
  TuningObj$WhichSamplerProgressInt <- WhichSamplerProgressInt
  TuningObj$WhichSGLDProgress <- WhichSGLDProgress
  TuningObj$WhichSGLDProgressInt <- WhichSGLDProgressInt
  return(TuningObj)
  
}



###Function for creating initial parameter object-------------------------------------------------------------------------------
CreatePara <- function(starting, DatObj) {

  ###Set data objects
  P <- DatObj$P
  Q <- DatObj$Q
  NL <- DatObj$NL
  NOmega <- DatObj$NOmega
  EyeQ <- DatObj$EyeQ

  ###Which parameters are user defined?
  UserStarters <- names(starting)

  ###Set initial values of Beta
  if ("Beta" %in% UserStarters) Beta <- matrix(starting$Beta, nrow = P, ncol = 1)
  if ((!"Beta" %in% UserStarters)) Beta <- matrix(rnorm(P), nrow = P, ncol = 1)

  ###Set initial values of l
  if ("l" %in% UserStarters) l <- matrix(starting$l, nrow = NL, ncol = 1)
  if ((!"l" %in% UserStarters)) l <- matrix(rnorm(NL), nrow = NL, ncol = 1)
  
  ###Set initial values of d
  if ("d" %in% UserStarters) d <- matrix(starting$d, nrow = Q, ncol = 1)
  if ((!"d" %in% UserStarters)) d <- matrix(rnorm(Q), nrow = Q, ncol = 1)
  
  ###Functions of parameters
  Z <- GetZ(l, Q)
  L <- GetL(Z, Q)
  LoverZ <- vecLT(L / Z)
  D <- diag(exp(as.numeric(d)))
  Upsilon <- L %*% t(L)
  Sigma <- D %*% Upsilon %*% D
  LInv <- forwardsolve(l = L, x = EyeQ)
  tLInv <- t(LInv)
  UpsilonInv <- tLInv %*% LInv
  DInv <- diag(1 / diag(D))
  SigmaInv <- DInv %*% UpsilonInv %*% DInv
  GradLZ <- diag(as.numeric(LoverZ), nrow = NL, ncol = NL)
  GradZl <- diag(as.numeric(1 / cosh(l)^2), nrow = NL, ncol = NL)
  GradLl <- GradLZ %*% GradZl
  
  ###Save parameter objects
  Para <- list()
  Para$Beta <- Beta
  Para$l <- l
  Para$d <- d
  Para$Omega <- rbind(Beta, l, d)
  Para$Z <- Z
  Para$L <- L
  Para$LoverZ <- LoverZ
  Para$D <- D
  Para$Upsilon <- Upsilon
  Para$Sigma <- Sigma
  Para$LInv <- LInv
  Para$tLInv <- tLInv
  Para$UpsilonInv <- UpsilonInv
  Para$DInv <- DInv
  Para$SigmaInv <- SigmaInv
  Para$GradLZ <- GradLZ
  Para$GradZl <- GradZl
  Para$GradLl <- GradLl
  Para$SigmaPrime <- matrix(0, nrow = NOmega, ncol = NOmega)
  return(Para)

}