##
## Time Series Stationary Test Functions
##
## Authors: Luis Gustavo Nardin
##          Vasile Alexandru Suchar
##
## Last Modification: 10/08/2016
##
rpackages <- commandArgs()

if(is.null(rpackages)){
  library(splus2R)
  library(ifultools)
  library(fractal)
  library(timeDate)
  library(timeSeries)
  library(fBasics)
  library(urca)
  library(fUnitRoots)
  library(wavethresh)
  library(locits)
  library(quadprog)
  library(zoo)
  library(tseries)
  library(costat)
  library(smoothmest)
  library(changepoint)
  library(lmtest)
  library(FinTS)
	library(strucchange)
} else {
  library(splus2R, lib=rpackages)
  library(ifultools, lib=rpackages)
  library(fractal, lib=rpackages)
  library(timeDate, lib=rpackages)
  library(timeSeries, lib=rpackages)
  library(fBasics, lib=rpackages)
  library(urca, lib=rpackages)
  library(fUnitRoots, lib=rpackages)
  library(wavethresh, lib=rpackages)
  library(locits, lib=rpackages)
  library(quadprog, lib=rpackages)
  library(zoo, lib=rpackages)
  library(tseries, lib=rpackages)
  library(costat, lib=rpackages)
  library(smoothmest, lib=rpackages)
	library(changepoint, lib=rpackages)
  library(lmtest, lib=rpackages)
  library(FinTS, lib=rpackages)
	library(strucchange, lib=rpackages)
}


###############
## CONSTANTS
###############
PARAM_ARMA      <- "arma"
PARAM_TREND     <- "trend"
PARAM_BREAK     <- "break"
PARAM_OPTION    <- c(PARAM_ARMA, PARAM_TREND, PARAM_BREAK)

ERROR_NORMAL    <- "normal"
ERROR_EXP       <- "exp"
ERROR_TSTUDENT  <- "studentt"
ERROR_OPTION    <- c(ERROR_NORMAL, ERROR_EXP, ERROR_TSTUDENT)

NONSTATIONARY   <- 0
STATIONARY      <- 1

ST              <- 0
NS_UNIT         <- 1
NS_MEAN         <- 2
NS_VAR          <- 3
NS_AC           <- 4


##
## AutoRegression Move Average ARMA(p,q) Time Series Generator (Normal Error PDF)
##
## Parameters
##    n      Number of elements
##    phi    Order of the autoregressive part
##    theta  Order of the moving average part
##    eCoeff Error multiplicative coefficient
##    miu    Mean value
##    sigma  Standard deviation
##
ARMANgenerator <- function(n, phi, theta, eCoeff, miu, sigma){
  l <- max(length(phi), length(theta))
  
  if(is.null(eCoeff)){
    epsilon <- rnorm(n, miu, sigma)
  } else {
    epsilon <- eCoeff * rnorm(n, miu, sigma)
  }
  
  ARMA <- vector()
  for (i in 1:n){
    if (i <= l){
      ARMA[i] <- epsilon[i]
    } else {
      s <- 0.0
      for (j in 1:length(phi)){
        s <- s + phi[j] * ARMA[i - j]
      }
      for (j in 1:length(theta)){
        s <- s + theta[j] * epsilon[i - j]
      }
      ARMA[i] <- s + epsilon[i]
    }
  }
  
  return(ARMA)
}


##
## AutoRegression Move Average ARMA(p,q) Time Series Generator (Double Exponential Error PDF)
##
## Parameters
##    n      Number of elements
##    phi    Order of the autoregressive part
##    theta  Order of the moving average part
##    eCoeff Error multiplicative coefficient
##    miu    Mean value
##    lambda Scale value
##
ARMAEgenerator <- function(n, phi, theta, eCoeff, miu, lambda){
  l <- max(length(phi), length(theta))
  
  if(is.null(eCoeff)){
    epsilon <- rdoublex(n, miu, lambda)
  } else {
    epsilon <- eCoeff * rdoublex(n, miu, lambda)
  }
  
  ARMA <- vector()
  for (i in 1:n){
    if (i <= l){
      ARMA[i] <- epsilon[i]
    } else {
      s <- 0.0
      for (j in 1:length(phi)){
        s <- s + phi[j] * ARMA[i - j]
      }
      for (j in 1:length(theta)){
        s <- s + theta[j] * epsilon[i - j]
      }
      ARMA[i] <- s + epsilon[i]
    }
  }
  
  return(ARMA)
}


##
## AutoRegression Move Average ARMA(p,q) Time Series Generator (T Test Error PDF)
##
## Parameters
##    n      Number of elements
##    phi    Order of the autoregressive part
##    theta  Order of the moving average part
##    eCoeff Error multiplicative coefficient
##    df     Degrees of freedom
##
ARMATgenerator <- function(n, phi, theta, eCoeff, df){
  l <- max(length(phi), length(theta))
  
  if(is.null(eCoeff)){
    epsilon <- rt(n, df)
  } else {
    epsilon <- eCoeff * rt(n, df)
  }
  
  ARMA <- vector()
  for (i in 1:n){
    if (i <= l){
      ARMA[i] <- epsilon[i]
    } else {
      s <- 0.0
      for (j in 1:length(phi)){
        s <- s + phi[j] * ARMA[i - j]
      }
      for (j in 1:length(theta)){
        s <- s + theta[j] * epsilon[i - j]
      }
      ARMA[i] <- s + epsilon[i]
    }
  }
  
  return(ARMA)
}


##
## ARMA(p,q)
##
tsARMA <- function(N, TS, p, q, error, seeds){
  ts <- array(0, dim=c(TS, N))
  for(n in 1:N){
    set.seed(seeds[n])
    if(error == ERROR_NORMAL){
      ts[,n] <- ARMANgenerator(TS, p, q, NULL, 0, 1)
    } else if(error == ERROR_EXP){
      ts[,n] <- ARMAEgenerator(TS, p, q, NULL, 0, 1)
    } else if (error == ERROR_TSTUDENT){
      ts[,n] <- ARMATgenerator(TS, p, q, NULL, 4)
    }
  }
  
  return(ts)
}


##
## Trend Mean ARMA(phi,theta) Time Series
##
## eparam: NORMAL   (1: Mean, 2: Variance)
##         EXP      (1: Mean, 2: Lambda)
##         TSTUDENT (1: Degree of Freedom)
##
tsTrendMean <- function(N, TS, trend, phi, theta, error, eparam, seeds){
  l <- max(length(phi), length(theta))
  
  ts <- array(0, dim=c(TS, N))
  for(n in 1:N){
    set.seed(seeds[n])
    
    if(error == ERROR_NORMAL){
      ts[1, n] <- rnorm(1, eparam[1], eparam[2])
      for(t in 2:TS){
        if(t <= l){
          ts[t, n] <- rnorm(1, eparam[1], eparam[2])
        } else {
          s <- 0.0
          for (j in 1:length(phi)){
            s <- s + phi[j] * ts[t - j, n]
          }
          for (j in 1:length(theta)){
            s <- s + theta[j] * rnorm(1, eparam[1], eparam[2])
          }
          ts[t, n] <- ((trend * t) / TS) + s + rnorm(1, eparam[1], eparam[2])
        }
      }
    } else if(error == ERROR_EXP){
      ts[1, n] <- rdoublex(1, eparam[1], eparam[2])
      for(t in 2:TS){
        if(t <= l){
          ts[t, n] <- rdoublex(1, eparam[1], eparam[2])
        } else {
          s <- 0.0
          for (j in 1:length(phi)){
            s <- s + phi[j] * ts[t - j, n]
          }
          for (j in 1:length(theta)){
            s <- s + theta[j] * rdoublex(1, eparam[1], eparam[2])
          }
          ts[t, n] <- ((trend * t) / TS) + s + rdoublex(1, eparam[1], eparam[2])
        }
      }
    } else if (error == ERROR_TSTUDENT){
      ts[1, n] <- rt(1, eparam[1])
      for(t in 2:TS){
        if(t <= l){
          ts[t, n] <- rt(1, eparam[1])
        } else {
          s <- 0.0
          for (j in 1:length(phi)){
            s <- s + phi[j] * ts[t - j, n]
          }
          for (j in 1:length(theta)){
            s <- s + theta[j] * rt(1, eparam[1])
          }
          ts[t, n] <- ((trend * t) / TS) + s + rt(1, eparam[1])
        }
      }
    }
  }
  
  return(ts)
}


##
## Trend Variance Time Series
##
## eparam: NORMAL   (1: Mean, 2: Variance)
##         EXP      (1: Mean, 2: Lambda)
##         TSTUDENT (1: Degree of Freedom)
##
tsTrendVariance <- function(N, TS, trend, phi, theta, error, eparam, seeds){
  ts <- array(0, dim=c(TS, N))
  for(n in 1:N){
    set.seed(seeds[n])
    coeff <- 1 + ((trend * (1:TS - 1)) / (TS - 1))
    if(error == ERROR_NORMAL){
      ts[,n] <- ARMANgenerator(TS, phi, theta, coeff, eparam[1], eparam[2])
    } else if(error == ERROR_EXP){
      ts[,n] <- ARMAEgenerator(TS, phi, theta, coeff, eparam[1], eparam[2])
    } else if (error == ERROR_TSTUDENT){
      ts[,n] <- ARMATgenerator(TS, phi, theta, coeff, eparam[1])
    }
  }
  
  return(ts)
}


##
## Trend AutoCorrelation Time Series
##
## eparam: NORMAL   (1: Mean, 2: Variance)
##         EXP      (1: Mean, 2: Lambda)
##         TSTUDENT (1: Degree of Freedom)
##
tsTrendAutoCorrelation <- function(N, TS, p, error, eparam, seeds){
  ts <- array(0, dim=c(TS, N))
  for(n in 1:N){
    set.seed(seeds[n])
    if(error == ERROR_NORMAL){
      epsilon <- rnorm(TS, mean=eparam[1], sd=eparam[2])
    } else if(error == ERROR_EXP){
      epsilon <- rdoublex(TS, mu=eparam[1], lambda=eparam[2])
    } else if (error == ERROR_TSTUDENT){
      epsilon <- rt(TS, df=eparam[1])
    }
    
    phi <- -(p + ((2 * p) / (TS - 1))) + (((2 * p) * 1:TS) / (TS - 1))
    ts[1, n] <- epsilon[1]
    for(i in 2:TS){
      ts[i, n] <- phi[i] * ts[i - 1, n] + epsilon[i]
    }
  }
  
  return(ts)
}


##
## Structural Break Mean Time Series
##
## eparam: NORMAL   (1: Mean, 2: Variance)
##         EXP      (1: Mean, 2: Lambda)
##         TSTUDENT (1: Degree of Freedom)
##
tsSBMean <- function(N, TS, value, phi, theta, error, eparam, seeds){
  ts <- array(0, dim=c(TS, N))
  aux <- c(rep(value, as.integer(TS / 2)), rep(0, (TS - as.integer(TS / 2))))
  for(n in 1:N){
    set.seed(seeds[n])
    if(error == ERROR_NORMAL){
      ts[,n] <- aux + ARMANgenerator(TS, phi, theta, NULL, eparam[1], eparam[2])
    } else if(error == ERROR_EXP){
      ts[,n] <- aux + ARMAEgenerator(TS, phi, theta, NULL, eparam[1], eparam[2])
    } else if (error == ERROR_TSTUDENT){
      ts[,n] <- aux + ARMATgenerator(TS, phi, theta, NULL, eparam[1])
    }
  }
  
  return(ts)
}


##
## Structural Break Variance Time Series
##
## eparam: NORMAL   (1: Mean, 2: Variance)
##
tsSBVariance <- function(N, TS, value, phi, theta, error, eparam, seeds){
  ts <- array(0, dim=c(TS, N))
  for(n in 1:N){
    set.seed(seeds[n])
    if(error == ERROR_NORMAL){
      ts1 <- ARMANgenerator(as.integer(TS / 2), phi, theta, NULL, eparam[1], eparam[2] + value)
      ts2 <- ARMANgenerator((TS - as.integer(TS / 2)), phi, theta, NULL, eparam[1], eparam[2])
      ts[,n] <- c(ts1, ts2)
    }
  }
  
  return(ts)
}


##
## Structural Break AutoCorrelation Time Series
##
## eparam: NORMAL   (1: Mean, 2: Variance)
##         EXP      (1: Mean, 2: Lambda)
##         TSTUDENT (1: Degree of Freedom)
##
tsSBAutoCorrelation <- function(N, TS, phi, theta, error, eparam, seeds){
  ts <- array(0, dim=c(TS, N))
  for(n in 1:N){
    set.seed(seeds[n])
    if(error == ERROR_NORMAL){
      ts1 <- ARMANgenerator(as.integer(TS / 2), phi[1], theta, NULL, eparam[1], eparam[2])
      ts2 <- ARMANgenerator((TS - as.integer(TS / 2)), phi[2], theta, NULL, eparam[1], eparam[2])
    } else if(error == ERROR_EXP){
      ts1 <- ARMAEgenerator(as.integer(TS / 2), phi[1], theta, NULL, eparam[1], eparam[2])
      ts2 <- ARMAEgenerator((TS - as.integer(TS / 2)), phi[2], theta, NULL, eparam[1], eparam[2])
    } else if (error == ERROR_TSTUDENT){
      ts1 <- ARMATgenerator(as.integer(TS / 2), phi[1], theta, NULL, eparam[1])
      ts2 <- ARMATgenerator((TS - as.integer(TS / 2)), phi[2], theta, NULL, eparam[1])
    }
    ts[,n] <- rbind(ts1, ts2)
  }
  
  return(ts)
}


##
## Augmented Dickey-Fuller (ADF)
##
## Package..: tseries
## Function.: adf.test
##
df <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- adf.test(data, alternative="stationary", k=0)$p.value,
           error=function(e){return(NA)})
  if (!is.na(p) && !is.null(p) && !is.nan(p)){
    if (p <= alpha){
      return(STATIONARY)
    } else {
      return(NONSTATIONARY)
    }
  }
  return(NA)
}


##
## Augmented Dickey-Fuller (ADF)
##
## Package..: tseries
## Function.: adf.test
##
adf <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- adf.test(data, alternative="stationary")$p.value,
           error=function(e){return(NA)})
  if (!is.na(p) && !is.null(p) && !is.nan(p)){
    if (p <= alpha){
      return(STATIONARY)
    } else {
      return(NONSTATIONARY)
    }
  }
  return(NA)
}


##
## Phillips-Perron (PP)
##
## Package..: tseries
## Function.: pp.test
##
pp <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- pp.test(data, alternative="stationary")$p.value,
           error=function(e){return(NA)})
  if (!is.na(p) && !is.null(p) && !is.nan(p)){
    if (p <= alpha){
      return(STATIONARY)
    } else {
      return(NONSTATIONARY)
    }
  }
  return(NA)
}


##
## Schmidt-Phillips (SP)
##
## Package..: fUnitRoots
## Function.: urspTest
##
sp <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- urspTest(data, signif=alpha, doplot=FALSE),
           error=function(e){return(NA)})
  if (!is.null(p)){
    if (!is.na(p@test$test@teststat)){
      if (p@test$test@teststat <= p@test$test@cval){
        return(STATIONARY)
      } else {
        return(NONSTATIONARY)
      }
    }
  }
  return(NA)
}


##
## Zivot-Andrews (ZA)
##
## Package..: fUnitRoots
## Function.: urzaTest
##
za <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- urzaTest(data, doplot=FALSE),
           error=function(e){return(NA)})
  if (!is.null(p)){
    if (!is.na(p@test$test@teststat)){
      if (p@test$test@teststat <= p@test$test@cval[2]){
        return(STATIONARY)
      } else {
        return(NONSTATIONARY)
      }
    }
  }
  return(NA)
}


##
## Elliot-Rothenberg-Stock (ERS)
##
## Package..: fUnitRoots
## Function.: urersTest
##
ers <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- urersTest(data, type="DF-GLS", model="constant", doplot=FALSE),
           error=function(e){return(NA)})
  if (!is.null(p)){
    if (!is.na(p@test$test@teststat)){
      if (p@test$test@teststat <= p@test$test@cval[2]){
        return(STATIONARY)
      } else {
        return(NONSTATIONARY)
      }
    }
  }
  return(NA)
}


##
## KPSS
##
## Package..: tseries
## Function.: kpss.test
##
kpss <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- kpss.test(data)$p.value,
           error=function(e){return(NA)})
  if (!is.na(p) && !is.null(p) && !is.nan(p)){
    if (p <= alpha){
      return(NONSTATIONARY)
    } else {
      return(STATIONARY)
    }
  }
  return(NA)
}


##
## Priestley-Subba Rao (PSR)
##
## Package..: fractal
## Function.: stationarity
##
psr <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- stationarity(data, significance=alpha),
           error=function(e){return(NA)})
  if (!is.null(p)){
    if ((attr(p, "pvals")[1]) <= alpha){
      return(NONSTATIONARY)
    } else {
      return(STATIONARY)
    }
  }
  return(NA)
}


##
## Wavelet
##
## Package..: locits
## Function.: hwtos2
##
wavelet <- function(data, alpha, ...){
  num <- log2(length(data))
  if ((num %% as.integer(num)) == 0){
    p <- NULL
    tryCatch(p <- hwtos2(data, alpha=alpha),
             error=function(e){return(NA)})
    if (!is.null(p)){
      if (p$nreject < 1){
        return(STATIONARY)
      } else {
        return(NONSTATIONARY)
      }
    }
  }
  return(NA)
}


##
## Bootstrap
##
## Package..: costat
## Function.: BootTOS
##
bootstrap <- function(data, alpha, ...){
  num <- log2(length(data))
  if ((num %% as.integer(num)) == 0){
    p <- NULL
    tryCatch(p <- BootTOS(data)$p.value,
             error=function(e){return(NA)})
    if (!is.null(p)){
      if (p <= alpha){
        return(NONSTATIONARY)
      } else {
        return(STATIONARY)
      }
    }
  }
  return(NA)
}


##
## Changepoint Mean
##
## Package..: changepoint
## Function.: cpt.mean
##
cpt.mean.test <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- cpt.mean(data, penalty="MBIC", method="AMOC"))
  
  if (!is.null(p)){
    if (ncpts(p) != 0){
      return(NONSTATIONARY)
    } else {
      return(STATIONARY)
    }
  }
  return(NA)
}


##
## Changepoint Variance
##
## Package..: changepoint
## Function.: cpt.var
##
cpt.var.test <- function(data, alpha, ...){
  p <- NULL
  tryCatch(p <- cpt.var(data, penalty="MBIC", method="AMOC"))
  
  if (!is.null(p)){
    if (ncpts(p) != 0){
      return(NONSTATIONARY)
    } else {
      return(STATIONARY)
    }
  }
  return(NA)
}


##
## ARMA and ARMA test
##
## Package..: tseries
## Function.: arma
##
arma.arma.test <- function(data, alpha, ...){
  
  model <- try(arma(data, order=c(1, 0)))
  
  if(class(model) != "try-error"){
    p <- NULL
    tryCatch(p <- summary(arma(na.omit(model$residuals^2),
                order=c(1, 0), lag=NULL, coef=NULL)))
              
    if(!is.null(p)){
      if(p$coef[1,4] <= alpha){
        return(NONSTATIONARY)
      } else {
        return(STATIONARY)
      }
    }
  }
  return(NA)
}


##
## Ljung–Box test statistic
##
## Package..: tseries, stats
## Function.: arma, Box.test
##
ljung.box.test <- function(data, alpha, ...){
  
  model <- try(arma(data, order=c(1, 0)))
  
  if(class(model) != "try-error"){
    p <- NULL
    tryCatch(p <- Box.test(model$residuals, lag=1, type="Ljung-Box"))
    
    if(!is.null(p)){
      if (p$p.value <= alpha){
        return(NONSTATIONARY)
      } else {
        return(STATIONARY)
      }
    }
  }
  return(NA)
}


##
## Breusch-Godfrey
##
## Package..: lmtest
## Function.: bgtest
##
breusch.godfrey.test <- function(data, alpha, ...){
  
  lag <- filter(data, c(0,1), method="conv", sides=1)
  
  p <- NULL
  tryCatch(p <- bgtest(data ~ 1 + lag, order=1, type="F", fill=NA))
  
  if (!is.null(p)){
    if (p$p.value <= alpha){
      return(NONSTATIONARY)
    } else {
      return(STATIONARY)
    }
  }
  return(NA)
}


##
## ArchTest
##
## Package..: tseries, FinTS
## Function.: arma, ArchTest
##
arch.test <- function(data, alpha, ...){
  
  model <- try(arma(data, order=c(1, 0)))
  
  if(class(model) != "try-error"){
    p <- NULL
    tryCatch(p <- ArchTest(model$residuals, lags=1, demean=FALSE))
    
    if (!is.null(p)){
      if (p$p.value <= alpha){
        return(NONSTATIONARY)
      } else {
        return(STATIONARY)
      }
    }
  }
  return(NA)
}


##
## Breusch-Pagan
##
## Package..: tseries, lmtest
## Function.: arma, bptest
##
breusch.pagan.test <- function(data, alpha, ...){
  
  model <- try(arma(data, order=c(1,0)))
  
  if(class(model) != "try-error"){
    lag <- filter(model$residuals, c(0,1), method="conv", sides=1)
    
    p <- NULL
    tryCatch(p <- bptest(model$residuals ~ 1 + lag))
    
    if (!is.null(p)){
      if (p$p.value <= alpha){
        return(NONSTATIONARY)
      } else {
        return(STATIONARY)
      }
    }
  }
  return(NA)
}


##
## Structural Change
##
## Package..: strucchange
## Function.: efp, sctest
## Parameter: window - search window
##
structure.change.test <- function(data, alpha, ...){
  arguments <- list(...)
  
  if((is.null(arguments)) | (length(arguments) == 0)){
    window <- 0.01
  } else {
    window <- as.double(arguments[[1]])
    
    if (window < 0){
      window <- 0
    } else if(window > 1){
      window <- 1
    }
  }
  
  lag <- filter(data, c(0, 1), method= "conv", sides=1)
  
  p <- NULL
  p <- try(efp(data ~ lag, h=window, type="ME"))
  
  if (class(p) != "try-error"){
    if(sctest(p)$p.value <= alpha){
      return(NONSTATIONARY)
    } else {
      return(STATIONARY)
    }
  }
  return(NA)
}


##
## Perform the test
##
runTest <- function(test, data, alpha, ...){
	args <- list(...)
	
  result <- 0
  
  testName <- tolower(substr(test, 1, 2))
  
  if(testName == "df"){
    result <- df(data, alpha)
  } else if(testName == "adf"){
    result <- adf(data, alpha)
  } else if(testName == "pp"){
    result <- pp(data, alpha)
  } else if(testName == "sp"){
    result <- sp(data, alpha)
  } else if(testName == "za"){
    result <- za(data, alpha)
  } else if(testName == "ers"){
    result <- ers(data, alpha)
  } else if(testName == "kpss"){
    result <- kpss(data, alpha)
  } else if(testName == "psr"){
    result <- psr(data, alpha)
  } else if(testName == "wavelet"){
    result <- wavelet(data, alpha)
  } else if(testName == "bootstrap"){
    result <- bootstrap(data, alpha)
  } else if(testName == "cptm"){
    result <- cpt.mean.test(data, alpha)
  } else if(testName == "cptv"){
	  result <- cpt.var.test(data, alpha)
  } else if(testName == "aa"){
    result <- arma.arma.test(data, alpha)
  } else if(testName == "lb"){
    result <- ljung.box.test(data, alpha)
  } else if(testName == "bg"){
    result <- breusch.godfrey.test(data, alpha)
  } else if(testName == "at"){
    result <- arch.test(data, alpha)
  } else if(testName == "bp"){
    result <- breusch.pagan.test(data, alpha)
  } else if(testName == "sc"){
    window <- NULL
    tryCatch(window <- as.double(substr(test, 3, nchar(test))) / 100)
    
    if(((is.null(args)) | (length(args) == 0)) &
        (is.na(window))){
      result <- structure.change.test(data, alpha)
    } else {
      if(!is.na(window)){
        result <- structure.change.test(data, alpha, window)
      } else {
        result <- structure.change.test(data, alpha, unlist(args))
      }
    }
  }
  
  return(result)
}


##
## Stationarity 2-steps decision tree
##
## mode 1 - KPSS-PSR+ERS
##      2 - KPSS-BOOTSTRAP+ERS
##      3 - ERS-PSR
##      4 - ERS-BOOTSTRAP
##      5 - KPSS-PSR
##      6 - KPSS-BOOTSTRAP
##
decisionTree2Steps <- function(data, alpha, mode){
  steps <- NULL
  
  ## KPSS-PSR+ERS
  if(mode == 1){
    result <- kpss(data, alpha)
    steps <- rbind(steps, c("KPSS", result))
    if(is.na(result)){
      return(list(NA, steps))
    }
    
    ## Stationarity
    if(result == 1){
      result <- psr(data, alpha)
      steps <- rbind(steps, c("PSR", result))
      
    ## Non-Stationary
    } else if(result == 0){
      result <- ers(data, alpha)
      steps <- rbind(steps, c("ERS", result))
    }
    
  ## KPSS-BOOTSTRAP+ERS
  } else if(mode == 2){
    result <- kpss(data, alpha)
    steps <- rbind(steps, c("KPSS", result))
    if(is.na(result)){
      return(list(NA, steps))
    }
    
    ## Stationarity
    if(result == 1){
      result <- bootstrap(data, alpha)
      steps <- rbind(steps, c("BOOTSTRAP", result))
      
    ## Non-Stationary
    } else if(result == 0){
      result <- ers(data, alpha)
      steps <- rbind(steps, c("ERS", result))
    }
    
  ## ERS-PSR
  } else if(mode == 3){
    result <- ers(data, alpha)
    steps <- rbind(steps, c("ERS", result))
    if(is.na(result)){
      return(list(NA, steps))
    }
    
    ## Stationarity
    if(result == 1){
      result <- psr(data, alpha)
      steps <- rbind(steps, c("PSR", result))
    }
    
  ## ERS-BOOTSTRAP
  } else if(mode == 4){
    result <- ers(data, alpha)
    steps <- rbind(steps, c("ERS", result))
    if(is.na(result)){
      return(list(NA, steps))
    }
    
    ## Stationarity
    if(result == 1){
      result <- bootstrap(data, alpha)
      steps <- rbind(steps, c("BOOTSTRAP", result))
    }
    
  ## KPSS-PSR
  } else if(mode == 5){
    result <- kpss(data, alpha)
    steps <- rbind(steps, c("KPSS", result))
    if(is.na(result)){
      return(list(NA, steps))
    }
    
    ## Stationarity
    if(result == 1){
      result <- psr(data, alpha)
      steps <- rbind(steps, c("PSR", result))
    }
    
  ## KPSS-BOOTSTRAP
  } else if(mode == 6){
    result <- kpss(data, alpha)
    steps <- rbind(steps, c("KPSS", result))
    if(is.na(result)){
      return(list(NA, steps))
    }
    
    ## Stationarity
    if(result == 1){
      result <- bootstrap(data, alpha)
      steps <- rbind(steps, c("BOOTSTRAP", result))
    }
  }
  
  steps <- data.frame(as.character(steps[,1]),
                      as.integer(steps[,2]), stringsAsFactors=FALSE)
  names(steps) <- c("Test", "Result")
  
  return(list(result, steps))
}


##
## Stationarity 3-steps decision tree
##
## mode 1 - KPSS-PSR-BOOTSTRAP+ERS
##
decisionTree3Steps <- function(data, alpha, mode){
  steps <- NULL
  
  if(mode == 1){
    result <- kpss(data, alpha)
    steps <- rbind(steps, c("KPSS", result))
    if(is.na(result)){
      return(list(NA, steps))
    }
    
    ## Stationary
    if(result == 1){
      
      result <- psr(data, alpha)
      steps <- rbind(steps, c("PSR", result))
      if(is.na(result)){
        return(list(NA, steps))
      }
      
      ## Stationary
      if(result == 1){
        result <- bootstrap(data, alpha)
        steps <- rbind(steps, c("BOOTSTRAP", result))
      }
      
    ## Non-Stationary
    } else if(result == 0){
      result <- ers(data, alpha)
      steps <- rbind(steps, c("ERS", result))
    }
  }
  
  steps <- data.frame(as.character(steps[,1]),
                      as.integer(steps[,2]), stringsAsFactors=FALSE)
  names(steps) <- c("Test", "Result")
  
  return(list(result, steps))
}


##
## Stationarity detection algorithm
##
## mode 1 - ALGO-ERS
##
detectStationarity <- function(data, alpha, mode){
  
  result <- NA
  
  # ERS tests Positive Unit-Root
  #           Trend Mean
  #           Break of large time series
  test <- NA
  test <- ers(data, alpha)
  
  if(!is.na(test)){
    if(test == NONSTATIONARY){
      result <- NONSTATIONARY
    } else {
      # BG tests Trend/Break AutoCorrelation
      #          Missclassified Break Mean and Trend Variance
      bg <- breusch.godfrey.test(data, alpha)
      
      if(!is.na(bg)){
        if(bg == NONSTATIONARY){
          result <- NONSTATIONARY
        } else {
          # AT tests Trend/Break Variance
          at <- arch.test(data, alpha)
          
          if(!is.na(at)){
            if(at == NONSTATIONARY){
              result <- NONSTATIONARY
            } else {
              result <- STATIONARY
            }
          }
        }
      }
    }
  }
  return(result)
}
