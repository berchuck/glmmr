###Function for reading in sampler inputs and creating a list object that contains all relavent data objects--------------------
CreateDatObj <- function(pformula, gformula, group, data, family) {

  ###Data objects
  N <- nrow(data) # total observations
  n <- length(unique(dat[, group])) #number of spatial locations
  Y <- matrix(data[, all.vars(pformula)[1]], ncol = 1)
  group <- data[, group] # group must be sorted or we have to sort it here.

  ###Covariates
  X <- model.matrix(pformula, data = data)
  Z_temp <- model.matrix(gformula, data = data)
  p <- ncol(X)
  q <- ncol (Z_temp)
  X <- matrix(X, nrow = N, ncol = p)
  Z <- bdiag(lapply(split(Z_temp, group), function(x) matrix(x, ncol = q)))
  Z <- matrix(Z, nrow = N, ncol = q * n)
  group2 <- rep(1:n, each = q)
    
  ###Matrix Objects
  EyeQ <- diag(q)
  
  ###Family indicator
  FamilyInd <- -1
  if (family == "normal") FamilyInd <- 0
  if (family == "binomial") FamilyInd <- 1
  if (family == "poisson") FamilyInd <- 2
  
  ###Make parameters global
  DatObj <- list()
  DatObj$Y <- Y
  DatObj$X <- X
  DatObj$Z <- Z
  DatObj$group <- group
  DatObj$group2 <- group2
  DatObj$N <- N
  DatObj$n <- n
  DatObj$p <- p
  DatObj$q <- q
  DatObj$n_L <- choose(q, 2)
  DatObj$FamilyInd <- FamilyInd
  DatObj$EyeQ <- EyeQ
  return(DatObj)

}



###Function to create Hyperparameter Object------------------------------------------------------------------------------------
CreateHyPara <- function(hypers, DatObj) {

  ###Set data objects
  # p <- DatObj$p
  
  ###Which parameters are user defined?
  Userhypers <- names(hypers)

  ###Set hyperparameters for Beta
  # if ("Beta" %in% Userhypers) {
  #   SigmaBeta <- hypers$Beta$SigmaBeta
  # }
  # if (!("Beta" %in% Userhypers)) {
  #   SigmaBeta <- diag(p)
  # }

  ###Set hyperparameters for L
  if ("L" %in% Userhypers) {
    Eta <- hypers$L$Eta
  }
  if (!("L" %in% Userhypers)) {
    Eta <- 1
  }
  
  ###Set hyperparameters for D
  if ("D" %in% Userhypers) {
    Nu <- hypers$D$Nu
  }
  if (!("D" %in% Userhypers)) {
    Nu <- 3
  }

  ###Create object for hyperparameters
  HyPara <- list()
  HyPara$Eta <- Eta
  HyPara$Nu <- Nu
  return(HyPara)

}



###Function for creating an object containing relevant SGD information---------------------------------------------------
CreateSgdObj <- function(sgd, DatObj) {
  
  ###Set data objects
  p <- DatObj$p
  q <- DatObj$q
  n_L <- DatObj$n_L
  
  ###Which parameters are user defined?
  UserSgds <- names(sgd)
  
  ###Set Epsilon
  if ("Epsilon" %in% UserSgds) Epsilon <- sgd$Epsilon
  if (!("Epsilon" %in% UserSgds)) Epsilon <- 1e-8
  
  ###Set Mu_nadam
  if ("Mu_nadam" %in% UserSgds) Mu_nadam <- sgd$Mu_nadam
  if (!("Mu_nadam" %in% UserSgds)) Mu_nadam <- 0.975
  
  ###Set Alpha_nadam
  if ("Alpha_nadam" %in% UserSgds) Alpha_nadam <- sgd$Alpha_nadam
  if (!("Alpha_nadam" %in% UserSgds)) Alpha_nadam <- 0.002
  
  ###Set Nu_nadam
  if ("Nu_nadam" %in% UserSgds) Nu_nadam <- sgd$Nu_nadam
  if (!("Nu_nadam" %in% UserSgds)) Nu_nadam <- 0.999
  
  ###Set S # size of the minibatches 
  if ("S" %in% UserSgds) S <- sgd$S
  if (!("S" %in% UserSgds)) S <- 10
  
  ###Set n_epochs
  if ("n_epochs" %in% UserSgds) n_epochs <- sgd$n_epochs
  if (!("n_epochs" %in% UserSgds)) n_epochs <- 50000

  ###Set R
  if ("R" %in% UserSgds) R <- sgd$R
  if (!("R" %in% UserSgds)) R <- 500
  
  ###Burn-in progress bar
  BarLength <- 50 #Burn-in bar length (arbitrary)
  BurnInProgress <- seq(1 / BarLength, 1, 1 / BarLength)
  WhichBurnInProgress <- sapply(BurnInProgress, function(x) tail(which(1 : n_epochs <= x * n_epochs), 1))
  SamplerProgress <- seq(0.1, 1.0, 0.1) #Intervals of progress update (arbitrary)
  WhichBurnInProgressInt <- sapply(SamplerProgress, function(x) tail(which(1:n_epochs <= x * n_epochs), 1))
  
  ###Return SGD object
  SgdObj <- list()
  SgdObj$Epsilon <- Epsilon
  SgdObj$M_nadam <- matrix(0, nrow = p + n_L + q)
  SgdObj$N_nadam <- matrix(0, nrow = p + n_L + q)
  SgdObj$Mu_nadam <- Mu_nadam
  SgdObj$Alpha_nadam <- Alpha_nadam
  SgdObj$Nu_nadam <- Nu_nadam
  SgdObj$S <- 10
  SgdObj$n_epochs <- 50000
  SgdObj$R <- 500
  SgdObj$BarLength <- BarLength
  SgdObj$WhichBurnInProgress <- WhichBurnInProgress
  SgdObj$WhichBurnInProgressInt <- WhichBurnInProgressInt
  return(SgdObj)
  
}



###Function for creating initial parameter object-------------------------------------------------------------------------------
CreatePara <- function(starting, DatObj, HyPara) {

  ###Set data objects
  p <- DatObj$p
  q <- DatObj$q
  n_L <- DatObj$n_L

  ###Set hyperparameter objects
  Eta <- HyPara$Eta
  Nu <- HyPara$Nu

  ###Which parameters are user defined?
  UserStarters <- names(starting)

  ###Set initial values of Beta
  if ("Beta" %in% UserStarters) Beta <- matrix(starting$Beta, nrow = p, ncol = 1)
  if ((!"Beta" %in% UserStarters)) Beta <- matrix(rnorm(p), nrow = p, ncol = 1)
  
  ###Set initial values of l
  if ("l" %in% UserStarters) l <- matrix(starting$l, nrow = n_L, ncol = 1)
  if ((!"l" %in% UserStarters)) l <- matrix(rnorm(n_L), nrow = n_L, ncol = 1)
  
  ###Set initial values of d
  if ("d" %in% UserStarters) d <- matrix(starting$d, nrow = q, ncol = 1)
  if ((!"d" %in% UserStarters)) d <- matrix(rnorm(q), nrow = q, ncol = 1)
  
  ###Functions of parameters
  z <- GetZ(l, q)
  L <- GetL(z, q)
  Loverz <- vecLT(L / z)
  D <- diag(exp(as.numeric(d)))
  Upsilon <- L %*% t(L)
  Sigma <- D %*% Upsilon %*% D
  LInv <- forwardsolve(l = L, x = diag(ncol(L)))
  tLInv <- t(LInv)
  UpsilonInv <- tLInv %*% LInv
  DInv <- diag(1 / diag(D))
  SigmaInv <- DInv %*% UpsilonInv %*% DInv
  gradLz <- diag(as.numeric(Loverz), nrow = n_L, ncol = n_L)
  gradzl <- diag(as.numeric(1 / cosh(l)^2), nrow = n_L, ncol = n_L)
  gradLl <- gradLz %*% gradzl
  
  ###Save parameter objects
  Para <- list()
  Para$Beta <- Beta
  Para$l <- l
  Para$d <- d
  Para$Omega <- rbind(Beta, l, d)
  Para$z <- z
  Para$L <- L
  Para$Loverz <- Loverz
  Para$D <- D
  Para$Upsilon <- Upsilon
  Para$Sigma <- Sigma
  Para$LInv <- LInv
  Para$tLInv <- t(LInv)
  Para$UpsilonInv <- UpsilonInv
  Para$DInv <- DInv
  Para$SigmaInv <- SigmaInv
  Para$gradLz <- gradLz
  Para$gradzl <- gradzl
  Para$gradLl <- gradLl
  return(Para)

}