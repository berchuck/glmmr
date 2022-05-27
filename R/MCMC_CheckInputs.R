CheckInputs <- function(pformula, gformula, group, data, family, algorithm, starting, hypers, tuning, seed, timer) {
  
  ###Data dimensions
  N <- nrow(data)
  
  ###Algorithm
  if (!algorithm %in% c("sgd", "sgld", "sgld_corrected", "gibbs")) stop('family: Must be one of "sgd", "sgld", "sgld_corrected", or "gibbs"')
  
  ###Family
  if (!family %in% c("normal", "bernoulli", "poisson")) stop('family: Must be one of "normal", "bernoulli", or "poisson"')

  ###Data and Formula
  if (!(class(pformula) == "formula")) stop('"pformula" must be of class formula')
  if (!(class(gformula) == "formula")) stop('"gformula" must be of class formula')
  if (!is.data.frame(data)) stop('"data" must be of class data.frame')
  pformula.test <- model.frame(pformula, data) # check that the variables in formula are in data (will return error)
  gformula.test <- model.frame(gformula, data) # check that the variables in formula are in data (will return error)
  if (!is.character(group)) stop('"group" must be a character string')
  if (length(group) > 1) stop('"group" must be length 1')
  group_dat <- data[, group] # will return error if group is not included in data
  NUnits <- length(unique(group_dat))
  
  ###Data checks for Y
  outcome <- all.vars(pformula)[1] # outcome variable
  Y <- pformula.test[, outcome]
  if (length(Y) != N) stop(paste0('Y must have exactly ', N, 'values'))
  if (any(is.na(Y))) stop("Y may have no missing values")
  if (any(!is.finite(Y))) stop("Y must have strictly finite entries")
  if (family == "bernoulli") if (any(Y != 0 & Y != 1)) stop('Y: for "bernoulli" observed data must contain only 0 and 1')
  if (family == "poisson") if (any(Y < 0)) stop('Y: for "poisson" observed data must be non-negative')
  if (family == "poisson") if (any(!is.integer(Y))) stop('Y: for "poisson" observed data must be an integer')
  
  ###Family indicator
  if (family == "normal") FamilyInd <- 0
  if (family == "bernoulli") FamilyInd <- 1
  if (family == "poisson") FamilyInd <- 2
  
  ###Data checks for covariates
  covariates <- all.vars(pformula)[-1]
  X <- as.matrix(pformula.test[, covariates])
  P <- dim(X)[2]
  if (any(is.na(X))) stop('Covariates in pformula cannot contain missing values')
  if (any(is.infinite(X))) stop('Covariates in pformula cannot contain infinite values')
  covariates <- all.vars(gformula)
  Z <- as.matrix(gformula.test[, covariates])
  Q <- dim(Z)[2]
  if (any(is.na(Z))) stop('Covariates in gformula cannot contain missing values')
  if (any(is.infinite(Z))) stop('Covariates in gformula cannot contain infinite values')

  ###Data checks for group variable
  if (any(is.na(group_dat))) stop('"group" variable cannot contain missing values')
  if (any(is.infinite(group_dat))) stop('"group" variable cannot contain infinite values')
  
  ###Timer
  if (!is.null(timer)) {
    if (!is.scalar(timer)) stop('timer must be a scalar')
    if (is.na(timer)) stop('timer cannot be NA')
    if (!is.finite(timer)) stop('timer cannot be infinite')
    if (!is.wholenumber(timer) | timer <= 0) stop('timer must be a positive integer')
  }
    
  ###hypers
  if (!is.null(hypers)) {
    if (!is.list(hypers)) stop('hypers must be a list')
    if (!all(names(hypers) %in% c("l", "d"))) stop('hypers: Can only contain lists with names "l" or "d"')

    ###If l hyperparameters are provided
    if ("l" %in% names(hypers)) {
      if (!is.list(hypers$l)) stop('hypers: "l" must be a list')
      if (!"Eta" %in% names(hypers$l)) stop('hypers: "Eta" value missing')
      if (!is.scalar(hypers$l$Eta)) stop('hypers: "Eta" must be a scalar')
      if (is.na(hypers$l$Eta)) stop('hypers: "Eta" cannot be NA')
      if (!is.finite(hypers$l$Eta)) stop('hypers: "Eta" cannot be infinite')
      if (hypers$l$Eta <= 0) stop('hypers: "Eta" must be strictly positive')
    }

    ###If d hyperparameters are provided
    if ("d" %in% names(hypers)) {
      if (!is.list(hypers$d)) stop('hypers: "d" must be a list')
      if (!"Nu" %in% names(hypers$d)) stop('hypers: "Nu" value missing')
      if (!is.scalar(hypers$d$Nu)) stop('hypers: "Nu" must be a scalar')
      if (is.na(hypers$d$Nu)) stop('hypers: "Nu" cannot be NA')
      if (!is.finite(hypers$d$Nu)) stop('hypers: "Nu" cannot be infinite')
      if (hypers$d$Nu <= 0) stop('hypers: "Nu" must be strictly positive')
    }

  ###End Hyperparameters
  }

  ###starting Values
  if (!is.null(starting)) {
    if (!is.list(starting)) stop('starting must be a list')
    if (!all(names(starting) %in% c("Beta", "l", "d"))) stop('starting: Can only contain objects with names "Beta", "l", or "d"')

    ###If Beta starting values is provided
    if ("Beta" %in% names(starting)) {
      if (is.vector(starting$Beta)) {
        if (!is.numeric(starting$Beta)) stop('starting: "Beta" must be a vector')
        if (length(starting$Beta) != P) stop('starting: "Beta" must be length P')
        if (!all(!is.na(starting$Beta))) stop('starting: "Beta" cannot have missing values')
        if (!all(is.finite(starting$Beta))) stop('starting: "Beta" cannot have infinite values')
      }
      if (is.scalar(starting$Beta)) {
        if (is.na(starting$Beta)) stop('starting: "Beta" cannot be NA')
        if (!is.finite(starting$Beta)) stop('starting: "Beta" cannot be infinite')
      }
      if ((!is.vector(starting$Beta)) & (!is.scalar(starting$Beta))) stop('starting: "Beta" must be a scalar or a vector')
    }
    
    ###If NL starting values is provided
    NL <- choose(Q, 2)
    if ("l" %in% names(starting)) {
      if (is.vector(starting$l)) {
        if (!is.numeric(starting$l)) stop('starting: "l" must be a vector')
        if (length(starting$l) != NL) stop('starting: "l" must be length NL')
        if (!all(!is.na(starting$l))) stop('starting: "l" cannot have missing values')
        if (!all(is.finite(starting$l))) stop('starting: "l" cannot have infinite values')
      }
      if (is.scalar(starting$l)) {
        if (is.na(starting$l)) stop('starting: "l" cannot be NA')
        if (!is.finite(starting$l)) stop('starting: "l" cannot be infinite')
      }
      if ((!is.vector(starting$l)) & (!is.scalar(starting$l))) stop('starting: "l" must be a scalar or a vector')
    }
    
    ###If d starting values is provided
    if ("d" %in% names(starting)) {
      if (is.vector(starting$d)) {
        if (!is.numeric(starting$d)) stop('starting: "d" must be a vector')
        if (length(starting$d) != Q) stop('starting: "d" must be length Q')
        if (!all(!is.na(starting$d))) stop('starting: "d" cannot have missing values')
        if (!all(is.finite(starting$d))) stop('starting: "d" cannot have infinite values')
      }
      if (is.scalar(starting$d)) {
        if (is.na(starting$d)) stop('starting: "d" cannot be NA')
        if (!is.finite(starting$d)) stop('starting: "d" cannot be infinite')
      }
      if ((!is.vector(starting$d)) & (!is.scalar(starting$d))) stop('starting: "d" must be a scalar or a vector')
    }
    
  ###End starting Values
  }

  ###Tuning Values
  if (!is.null(tuning)) {
    if (!is.list(tuning)) stop('tuning must be a list')
    if (!all(names(tuning) %in% c("NPilot", "TuneD", "TuneL", "EpsilonNADAM", "MuNADAM", "AlphaNADAM", "NuNADAM", "S", "S_SGLD", "NEpochs", "R", "EpsilonSGLD", "NSims", "NThin", "NTune", "NTune_seconds"))) stop('tuning: Can only contain objects with names "EpsilonNADAM", "MuNADAM", "AlphaNADAM", "NuNADAM", "S", "S_SGLD, "NEpochs", "R", "EpsilonSGLD", "NSims", "NThin", "NPilot", "TuneL", "TuneD", "NTune", or "NTune_seconds"')
    
    ###If L tuning values are provided
    if ("TuneL" %in% names(tuning)) {
      if (!is.numeric(tuning$TuneL) | !is.scalar(tuning$TuneL)) stop('tuning: "TuneL" must a vector or scalar')
      if (length(tuning$TuneL) > 1) {
        if (length(tuning$TuneL) != M) stop(paste0('tuning: "TuneL" must have length ', NL))
        if (!all(!is.na(tuning$TuneL))) stop('tuning: "TuneL" cannot have missing values')
        if (!all(is.finite(tuning$TuneL))) stop('tuning: "TuneL" cannot have infinite values')
        if (any(tuning$TuneL < 0)) stop('tuning: "TuneL" must have non-negative components')
      }
      if (is.scalar(tuning$TuneL)) {
        if (!is.scalar(tuning$TuneL)) stop('tuning: "TuneL" must be a scalar')
        if (is.na(tuning$TuneL)) stop('tuning: "TuneL" cannot be NA')
        if (!is.finite(tuning$TuneL)) stop('tuning: "TuneL" cannot be infinite')
        if (tuning$TuneL < 0) stop('tuning: "TuneL" must be non-negative')
      }
    }
    
    ###If D tuning values are provided
    if ("TuneD" %in% names(tuning)) {
      if (!is.numeric(tuning$TuneD) | !is.scalar(tuning$TuneD)) stop('tuning: "TuneD" must a vector or scalar')
      if (length(tuning$TuneD) > 1) {
        if (length(tuning$TuneD) != M) stop(paste0('tuning: "TuneD" must have length ', Q))
        if (!all(!is.na(tuning$TuneD))) stop('tuning: "TuneD" cannot have missing values')
        if (!all(is.finite(tuning$TuneD))) stop('tuning: "TuneD" cannot have infinite values')
        if (any(tuning$TuneD < 0)) stop('tuning: "TuneD" must have non-negative components')
      }
      if (is.scalar(tuning$TuneD)) {
        if (!is.scalar(tuning$TuneD)) stop('tuning: "TuneD" must be a scalar')
        if (is.na(tuning$TuneD)) stop('tuning: "TuneD" cannot be NA')
        if (!is.finite(tuning$TuneD)) stop('tuning: "TuneD" cannot be infinite')
        if (tuning$TuneD < 0) stop('tuning: "TuneD" must be non-negative')
      }
    }
    
    ###If EpsilonNADAM tuning value is provided
    if ("EpsilonNADAM" %in% names(tuning)) {
      if (!is.scalar(tuning$EpsilonNADAM)) stop('tuning: "EpsilonNADAM" must be a scalar')
      if (is.na(tuning$EpsilonNADAM)) stop('tuning: "EpsilonNADAM" cannot be NA')
      if (!is.finite(tuning$EpsilonNADAM)) stop('tuning: "EpsilonNADAM" cannot be infinite')
      if (tuning$EpsilonNADAM < 0) stop('tuning: "EpsilonNADAM" must be non-negative')
    }
    
    ###If MuNADAM tuning value is provided
    if ("MuNADAM" %in% names(tuning)) {
      if (!is.scalar(tuning$MuNADAM)) stop('tuning: "MuNADAM" must be a scalar')
      if (is.na(tuning$MuNADAM)) stop('tuning: "MuNADAM" cannot be NA')
      if (!is.finite(tuning$MuNADAM)) stop('tuning: "MuNADAM" cannot be infinite')
      if (tuning$MuNADAM < 0 | tuning$MuNADAM > 1) stop('tuning: "MuNADAM" must be in [0, 1]')
    }
    
    ###If NuNADAM tuning value is provided
    if ("NuNADAM" %in% names(tuning)) {
      if (!is.scalar(tuning$NuNADAM)) stop('tuning: "NuNADAM" must be a scalar')
      if (is.na(tuning$NuNADAM)) stop('tuning: "NuNADAM" cannot be NA')
      if (!is.finite(tuning$NuNADAM)) stop('tuning: "NuNADAM" cannot be infinite')
      if (tuning$NuNADAM < 0 | tuning$NuNADAM > 1) stop('tuning: "NuNADAM" must be in [0, 1]')
    }
    
    ###If AlphaNADAM tuning value is provided
    if ("AlphaNADAM" %in% names(tuning)) {
      if (!is.scalar(tuning$AlphaNADAM)) stop('tuning: "AlphaNADAM" must be a scalar')
      if (is.na(tuning$AlphaNADAM)) stop('tuning: "AlphaNADAM" cannot be NA')
      if (!is.finite(tuning$AlphaNADAM)) stop('tuning: "AlphaNADAM" cannot be infinite')
      if (tuning$AlphaNADAM < 0) stop('tuning: "AlphaNADAM" must be non-negative')
    }
    
    ###If EpsilonSGLD tuning value is provided
    if ("EpsilonSGLD" %in% names(tuning)) {
      if (!is.scalar(tuning$EpsilonSGLD)) stop('tuning: "EpsilonSGLD" must be a scalar')
      if (is.na(tuning$EpsilonSGLD)) stop('tuning: "EpsilonSGLD" cannot be NA')
      if (!is.finite(tuning$EpsilonSGLD)) stop('tuning: "EpsilonSGLD" cannot be infinite')
      if (tuning$EpsilonSGLD < 0) stop('tuning: "EpsilonSGLD" must be non-negative')
    }
    
    ###If NBurn is provided
    if ("NEpochs" %in% names(tuning)) {
      if (!is.scalar(tuning$NEpochs)) stop('tuning: "NEpochs" must be a scalar')
      if (is.na(tuning$NEpochs)) stop('tuning: "NEpochs" cannot be NA')
      if (!is.finite(tuning$NEpochs)) stop('tuning: "NEpochs" cannot be infinite')
      if (!is.wholenumber(tuning$NEpochs) | tuning$NEpochs < 0) stop('tuning: "NEpochs" must be a non-negative integer')
      if (tuning$NEpochs < 1) stop('tuning: "NEpochs" must be at least 1')
    }
    
    ###If NSims is provided
    if ("NSims" %in% names(tuning)) {
      if (!is.scalar(tuning$NSims)) stop('tuning: "NSims" must be a scalar')
      if (is.na(tuning$NSims)) stop('tuning: "NSims" cannot be NA')
      if (!is.finite(tuning$NSims)) stop('tuning: "NSims" cannot be infinite')
      if (!is.wholenumber(tuning$NSims) | tuning$NSims <= 0) stop('tuning: "NSims" must be a positive integer')
      if (tuning$NSims < 100) stop('tuning: "NSims" must be at least 100')
    }
    
    ###If NTune is provided
    if ("NTune" %in% names(tuning)) {
      if (!is.scalar(tuning$NTune)) stop('tuning: "NTune" must be a scalar')
      if (is.na(tuning$NTune)) stop('tuning: "NTune" cannot be NA')
      if (!is.finite(tuning$NTune)) stop('tuning: "NTune" cannot be infinite')
      if (!is.wholenumber(tuning$NTune) | tuning$NTune <= 0) stop('tuning: "NTune" must be a positive integer')
      if (tuning$NTune < 100) stop('tuning: "NTune" must be at least 100')
    }
    
    ###If NTune_seconds tuning value is provided
    if ("NTune_seconds" %in% names(tuning)) {
      if (!is.scalar(tuning$NTune_seconds)) stop('tuning: "NTune_seconds" must be a scalar')
      if (is.na(tuning$NTune_seconds)) stop('tuning: "NTune_seconds" cannot be NA')
      if (!is.finite(tuning$NTune_seconds)) stop('tuning: "NTune_seconds" cannot be infinite')
      if (tuning$NTune_seconds <= 0) stop('tuning: "NTune_seconds" must be positive')
    }
    
    ###If NThin is provided
    if ("NThin" %in% names(tuning)) {
      if (!is.scalar(tuning$NThin)) stop('tuning: "NThin" must be a scalar')
      if (is.na(tuning$NThin)) stop('tuning: "NThin" cannot be NA')
      if (!is.finite(tuning$NThin)) stop('tuning: "NThin" cannot be infinite')
      if (!is.wholenumber(tuning$NThin) | tuning$NThin <= 0) stop('tuning: "NThin" must be a positive integer')
      # if (!is.wholenumber(tuning$NSims / tuning$NThin)) stop('tuning: "NThin" must be a factor of "NSims"') enforced in Createtuning();
    }
    
    ###If NPilot is provided
    if ("NPilot" %in% names(tuning)) {
      if (!is.scalar(tuning$NPilot)) stop('tuning: "NPilot" must be a scalar')
      if (is.na(tuning$NPilot)) stop('tuning: "NPilot" cannot be NA')
      if (!is.finite(tuning$NPilot)) stop('tuning: "NPilot" cannot be infinite')
      if (!is.wholenumber(tuning$NPilot) | tuning$NPilot < 0) stop('tuning: "NPilot" must be a positive integer')
      # if (!is.wholenumber(tuning$NBurn / tuning$NPilot)) stop('tuning: "NPilot" must be a factor of "NBurn"') enforced in Createmcmc();
    }
    
    ###If R is provided
    if ("R" %in% names(tuning)) {
      if (!is.scalar(tuning$R)) stop('tuning: "R" must be a scalar')
      if (is.na(tuning$R)) stop('tuning: "R" cannot be NA')
      if (!is.finite(tuning$R)) stop('tuning: "R" cannot be infinite')
      if (!is.wholenumber(tuning$R) | tuning$R <= 0) stop('tuning: "R" must be a positive integer')
      if (tuning$R < 10) stop('tuning: "R" must be at least 10')
    }
    
    ###If S is provided
    if ("S" %in% names(tuning)) {
      if (!is.scalar(tuning$S)) stop('tuning: "S" must be a scalar')
      if (is.na(tuning$S)) stop('tuning: "S" cannot be NA')
      if (!is.finite(tuning$S)) stop('tuning: "S" cannot be infinite')
      if (!is.wholenumber(tuning$S) | tuning$S <= 0) stop('tuning: "S" must be a positive integer')
      if (!is.wholenumber(tuning$S) | tuning$S > NUnits) stop('tuning: "S" must be less than the number of unique values in group')
    }
    
    ###If S_SGLD is provided
    if ("S_SGLD" %in% names(tuning)) {
      if (!is.scalar(tuning$S_SGLD)) stop('tuning: "S_SGLD" must be a scalar')
      if (is.na(tuning$S_SGLD)) stop('tuning: "S_SGLD" cannot be NA')
      if (!is.finite(tuning$S_SGLD)) stop('tuning: "S_SGLD" cannot be infinite')
      if (!is.wholenumber(tuning$S_SGLD) | tuning$S_SGLD <= 0) stop('tuning: "S_SGLD" must be a positive integer')
      if (!is.wholenumber(tuning$S_SGLD) | tuning$S_SGLD > NUnits) stop('tuning: "S_SGLD" must be less than the number of unique values in group')
    }
    
  ###End tuning Values
  }

}

###Helper Functions
is.scalar <- function(x) ((is.numeric(x)) & (length(x) == 1))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
