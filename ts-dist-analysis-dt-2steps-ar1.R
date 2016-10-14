##
## Description: Execute 2 steps decision tree stationarity tests on different
##              time series ARMA (NORMAL, EXP, TSTUDENT), TREND and BREAK
##              (Distributed PBS_ARRAY)
##
## Authors: Luis Gustavo Nardin
##          Vasile Alexandru Suchar
##
## Last Modification: 10/07/2016
##
library(caTools)
library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(cores=8)


#############
## COMMAND-LINE
#############
args <- commandArgs(TRUE)

paramID <- as.numeric(as.character(args[1]))
errorID <- as.numeric(as.character(args[2]))
itemID <- as.numeric(as.character(args[3]))


#############
## PATHS
#############
baseDir <- "/scratch/nardluis"
rpackages <- paste0(baseDir, "/rpackages")
scriptDir <- paste0(baseDir, "/scripts/stationarity")
inputDir <- paste0(baseDir, "/data/stationarity")
outputDir <- paste0(baseDir, "/data/stationarity")


#############
## FUNCTIONS
#############
commandArgs <- function(){ rpackages }
source(paste0(scriptDir, "/functions.R"))


###################
## SET SEED
###################
fname <- paste0(inputDir, "/seeds.csv")
if(file.exists(fname)){
  s <- read.table(fname, sep=";", header=FALSE)
  seeds <- array(0, dim=nrow(s))
  for(i in 1:nrow(s)){
    seeds[i] <- s[i,]
  }
} else {
  seeds <- array(0, dim=N)
  s <- 1
  while(s <= N){
    seed <- as.integer(runif(1, min=0, max=(N*10)))
    seeds[s] <- seed
    s <- s + 1
  }
  
  s <- data.table(seeds)
  fname <- write.table(s, file=fname, quote=FALSE, sep=";",
                       row.names=FALSE, col.names=FALSE)
}


###################
## INPUT PARAMETERS
###################

## Number of replications
N <- 1000

## Largest size
TS <- 5000

## alpha
alpha <- 0.05

## Stationary Tests
tests <- c("DT-KPSS-PSR+ERS", "DT-KPSS-BOOTSTRAP+ERS", "DT-ERS-PSR",
           "DT-ERS-BOOTSTRAP", "DT-KPSS-PSR", "DT-KPSS-BOOTSTRAP")


##
## Stationary - ARMA(p,q)
##
if(PARAM_OPTION[paramID] == PARAM_ARMA){
  ## phi values
  values <- c(0, 0.99, 0.9, 0.8, 0.5, 0.2, -0.2, -0.5, -0.8, -0.9, -0.99, 1, -1)
  
  phi <- values[itemID]
  theta <- 0
  
  ## Window Size
  windows <- sort(c(20, seq(30, TS, 10),
                    c(32, 64, 128, 256, 512, 1024, 2048, 4096)))
  
  ## Header
  header <- paste0("W", as.character(windows))
  
  ## Error
  error <- ERROR_OPTION[errorID]
  
	## Create directory
  tsOutputDir <- paste0(outputDir, "/dt-2steps-arma-", ERROR_OPTION[errorID])
	dir.create(tsOutputDir, showWarnings=FALSE)
  
  armaResults <- array(NA, dim=c(length(tests), length(windows), N))
  armaTime <- array(NA, dim=c(length(tests), length(windows)))
  
  ts <- tsARMA(N, max(TS, windows), phi, theta, error, seeds)
  
  for(tindex in 1:length(tests)){
    for(windex in 1:length(windows)){
      beginT <- proc.time()
      aux <- foreach(n=1:N, .combine=c) %dopar%
        decisionTree2Steps(ts[1:windows[windex], n], alpha, tindex)[[1]]
      armaResults[tindex, windex,] <- aux
      endT <- proc.time()
      armaTime[tindex, windex] <- as.numeric((endT - beginT)[3])
    }
    
    x <- data.table(t(armaResults[tindex,,]))
    names(x) <- header
    fname <- paste0(tsOutputDir, "/arma-p", phi, "-q", theta, "-",
                    tests[tindex], "-", error, ".csv")
    write.table(x, file=fname, sep=";", quote=FALSE,
                row.names=FALSE, col.names=TRUE)
  }
  
  ## Write the execution time of the tests over Normal time series
  x <- data.table(cbind(tests, armaTime))
  names(x) <- c("Test", header)
  
  fname <- paste0(tsOutputDir, "/armaTime-p", phi, "-q", theta, "-", error,
                  ".csv")
  if (file.exists(fname)){
    write.table(x, file=fname, sep=";", quote=FALSE,
                row.names=FALSE, col.names=FALSE, append=TRUE)
  } else {
    write.table(x, file=fname, sep=";", quote=FALSE,
                row.names=FALSE, col.names=TRUE)
  }
}


##
## Non-Stationary - TREND
##
if(PARAM_OPTION[paramID] == PARAM_TREND){
  ## phi values
  phis <- c(0, 0.99, 0.9, 0.8, 0.5, 0.2, -0.2, -0.5, -0.8, -0.9, -0.99, 1, -1)
  phi <- phis[itemID]
  theta <- 0
  
  ## Window Size
  windows <- c(32, 64, 128, 256, 512, 1024, 2048, 4096)
  
  ## Header
  header <- paste0("W", as.character(windows))
  
  ## Error
  error <- ERROR_NORMAL
  
	## Create directory
  tsOutputDir <- paste0(outputDir, "/dt-2steps-", PARAM_OPTION[paramID])
	dir.create(tsOutputDir, showWarnings=FALSE)
  
	
  ##
  ## Trend MEAN
  ##
  values <- c(0, 0.5, 1, 2, 3, 4)
  
  for(trendM in values){
    nsResults <- array(NA, dim=c(length(tests), length(windows), N))
    nsTime <- array(NA, dim=c(length(tests), length(windows)))
    
    ts <- tsTrendMean(N, max(TS, windows), trendM, phi, theta, error,
                      c(0, 1), seeds)
    
    for(tindex in 1:length(tests)){
      for(windex in 1:length(windows)){
        beginT <- proc.time()
        aux <- foreach(n=1:N, .combine=c) %dopar%
          decisionTree2Steps(ts[1:windows[windex], n], alpha, tindex)[[1]]
        nsResults[tindex, windex,] <- aux
        endT <- proc.time()
        nsTime[tindex, windex] <- as.numeric((endT - beginT)[3])
      }
      
      x <- data.table(t(nsResults[tindex,,]))
      names(x) <- header
      fname <- paste0(tsOutputDir, "/trendM", trendM, "-phi", phi, "-",
                      tests[tindex], "-", error, ".csv")
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=TRUE)
    }
    
    ## Write the execution time of the tests over Normal time series
    x <- data.table(cbind(tests, nsTime))
    names(x) <- c("Test", header)
    
    fname <- paste0(tsOutputDir, "/trendMTime", trendM, "-phi", phi, "-",
                    error, ".csv")
    if (file.exists(fname)){
      write.table(x, file=fname,
                  sep=";", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    } else {
      write.table(x, file=fname,
                  sep=";", quote=FALSE, row.names=FALSE, col.names=TRUE)
    }
  }


  ##
  ## Trend VARIANCE
  ##
  values <- c(0, 0.5, 1, 2, 3, 4)
  
  for(trendV in values){
    nsResults <- array(NA, dim=c(length(tests), length(windows), N))
    nsTime <- array(NA, dim=c(length(tests), length(windows)))
    
    ts <- tsTrendVariance(N, max(TS, windows), trendV, phi, theta, error,
                          c(0, 1), seeds)
    
    for(tindex in 1:length(tests)){
      for(windex in 1:length(windows)){
        beginT <- proc.time()
        aux <- foreach(n=1:N, .combine=c) %dopar%
          decisionTree2Steps(ts[1:windows[windex], n], alpha, tindex)[[1]]
        nsResults[tindex, windex,] <- aux
        endT <- proc.time()
        nsTime[tindex, windex] <- as.numeric((endT - beginT)[3])
      }
      
      x <- data.table(t(nsResults[tindex,,]))
      names(x) <- header
      fname <- paste0(tsOutputDir, "/trendV", trendV, "-phi", phi, "-",
                      tests[tindex], "-", error, ".csv")
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=TRUE)
    }
    
    ## Write the execution time of the tests over Normal time series
    x <- data.table(cbind(tests, nsTime))
    names(x) <- c("Test", header)
    
    fname <- paste0(tsOutputDir, "/trendVTime", trendV, "-phi", phi, "-",
                    error, ".csv")
    if (file.exists(fname)){
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=FALSE, append=TRUE)
    } else {
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=TRUE)
    }
  }


  ##
  ## Trend AUTOCORRELATION
  ##
  if(itemID == 1){
    values <- c(0, 0.9, 0.45, -0.9, -0.45)
    
    for(trendAC in values){
      nsResults <- array(NA, dim=c(length(tests), length(windows), N))
      nsTime <- array(NA, dim=c(length(tests), length(windows)))
      
      for(tindex in 1:length(tests)){
        for(windex in 1:length(windows)){
          ts <- tsTrendAutoCorrelation(N, windows[windex], trendAC, error,
                                       c(0, 1), seeds)
          
          beginT <- proc.time()
          aux <- foreach(n=1:N, .combine=c) %dopar%
            decisionTree2Steps(ts[1:windows[windex],n], alpha, tindex)[[1]]
          nsResults[tindex, windex,] <- aux
          endT <- proc.time()
          nsTime[tindex, windex] <- as.numeric((endT - beginT)[3])
        }
        
        x <- data.table(t(nsResults[tindex,,]))
        names(x) <- header
        fname <- paste0(tsOutputDir, "/trendAC", trendAC, "-", tests[tindex],
                        "-", error, ".csv")
        write.table(x, file=fname, sep=";", quote=FALSE,
                    row.names=FALSE, col.names=TRUE)
      }
      
      ## Write the execution time of the tests over Normal time series
      x <- data.table(cbind(tests, nsTime))
      names(x) <- c("Test", header)
      
      fname <- paste0(tsOutputDir, "/trendACTime", trendAC, "-", error, ".csv")
      if (file.exists(fname)){
        write.table(x, file=fname, sep=";", quote=FALSE,
                    row.names=FALSE, col.names=FALSE, append=TRUE)
      } else {
        write.table(x, file=fname, sep=";", quote=FALSE,
                    row.names=FALSE, col.names=TRUE)
      }
    }
  }
}


##
## Non-Stationary - BREAK
##
if(PARAM_OPTION[paramID] == PARAM_BREAK){
  ## phi values
  phis <- c(0, 0.99, 0.9, 0.8, 0.5, 0.2, -0.2, -0.5, -0.8, -0.9, -0.99, 1, -1)
  phi <- phis[itemID]
  theta <- 0
  
  ## Window Size
  windows <- c(32, 64, 128, 256, 512, 1024, 2048, 4096)
  
  ## Header
  header <- paste0("W", as.character(windows))
  
  ## Error
  error <- ERROR_NORMAL
  
	## Create directory
  tsOutputDir <- paste0(outputDir, "/dt-2steps-", PARAM_OPTION[paramID])
	dir.create(tsOutputDir, showWarnings=FALSE)
	
	
  ##
  ## Break MEAN
  ##
  values <- c(0, 0.5, 1, 2, 3, 4)
  
  for(sbMean in values){
    nsResults <- array(NA, dim=c(length(tests), length(windows), N))
    nsTime <- array(NA, dim=c(length(tests), length(windows)))
    
    for(tindex in 1:length(tests)){
      for(windex in 1:length(windows)){
        ts <- tsSBMean(N, windows[windex], sbMean, phi, theta, error,
                       c(0, 1), seeds)
        
        beginT <- proc.time()
        aux <- foreach(n=1:N, .combine=c) %dopar%
          decisionTree2Steps(ts[1:windows[windex], n], alpha, tindex)[[1]]
        nsResults[tindex, windex,] <- aux
        endT <- proc.time()
        nsTime[tindex, windex] <- as.numeric((endT - beginT)[3])
      }
      
      x <- data.table(t(nsResults[tindex,,]))
      names(x) <- header
      fname <- paste0(tsOutputDir, "/breakM", sbMean, "-phi", phi, "-",
                      tests[tindex], "-", error, ".csv")
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=TRUE)
    }
    
    ## Write the execution time of the tests over Normal time series
    x <- data.table(cbind(tests, nsTime))
    names(x) <- c("Test", header)
    
    fname <- paste0(tsOutputDir, "/breakMTime", sbMean, "-phi", phi, "-",
                    error, ".csv")
    if (file.exists(fname)){
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=FALSE, append=TRUE)
    } else {
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=TRUE)
    }
  }
  
  
  ##
  ## Break VARIANCE
  ##
  values <- c(0, 0.5, 1, 2, 3, 4)
  
  for(sbVar in values){
    nsResults <- array(NA, dim=c(length(tests), length(windows), N))
    nsTime <- array(NA, dim=c(length(tests), length(windows)))
  
    for(tindex in 1:length(tests)){
      for(windex in 1:length(windows)){
        ts <- tsSBVariance(N, windows[windex], sbVar, phi, theta, error,
                           c(0, 1), seeds)
        
        beginT <- proc.time()
        aux <- foreach(n=1:N, .combine=c) %dopar%
          decisionTree2Steps(ts[1:windows[windex], n], alpha, tindex)[[1]]
        nsResults[tindex, windex,] <- aux
        endT <- proc.time()
        nsTime[tindex, windex] <- as.numeric((endT - beginT)[3])
      }
      
      x <- data.table(t(nsResults[tindex,,]))
      names(x) <- header
      fname <- paste0(tsOutputDir, "/breakV", sbVar, "-phi", phi, "-",
                      tests[tindex], "-", error, ".csv")
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=TRUE)
    }
    
    ## Write the execution time of the tests over Normal time series
    x <- data.table(cbind(tests, nsTime))
    names(x) <- c("Test", header)
    
    fname <- paste0(tsOutputDir, "/breakVTime", sbVar, "-phi", phi, "-",
                    error, ".csv")
    if (file.exists(fname)){
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=FALSE, append=TRUE)
    } else {
      write.table(x, file=fname, sep=";", quote=FALSE,
                  row.names=FALSE, col.names=TRUE)
    }
  }
  
  
  ##
  ## Break AUTOCORRELATION
  ##
  if(itemID == 1){
    values <- c(0.9, 0.45, 0, -0.45, -0.9)
    
    items <- combs(values, 2)
    
    for(item in 1:nrow(items)){
      trendAC1 <- items[item, 1]
      trendAC2 <- items[item, 2]
      
      nsResults <- array(NA, dim=c(length(tests), length(windows), N))
      nsTime <- array(NA, dim=c(length(tests), length(windows)))
      
      for(tindex in 1:length(tests)){
        for(windex in 1:length(windows)){
          ts <- tsSBAutoCorrelation(N, windows[windex], c(trendAC1,trendAC2),
                                    theta, error, c(0, 1), seeds)
          
          beginT <- proc.time()
          aux <- foreach(n=1:N, .combine=c) %dopar%
            decisionTree2Steps(ts[1:windows[windex], n], alpha, tindex)[[1]]
          nsResults[tindex, windex,] <- aux
          endT <- proc.time()
          nsTime[tindex, windex] <- as.numeric((endT - beginT)[3])
        }
        
        x <- data.table(t(nsResults[tindex,,]))
        names(x) <- header
        fname <- paste0(tsOutputDir, "/breakAC", trendAC1, "-", trendAC2, "-",
                        tests[tindex], "-", error, ".csv")
        write.table(x, file=fname, sep=";", quote=FALSE,
                    row.names=FALSE, col.names=TRUE)
      }
      
      ## Write the execution time of the tests over Normal time series
      x <- data.table(cbind(tests, nsTime))
      names(x) <- c("Test", header)
      
      fname <- paste0(tsOutputDir, "/breakACTime", trendAC1, "-", trendAC2, "-",
                      error, ".csv")
      if (file.exists(fname)){
        write.table(x, file=fname, sep=";", quote=FALSE,
                    row.names=FALSE, col.names=FALSE, append=TRUE)
      } else {
        write.table(x, file=fname, sep=";", quote=FALSE,
                    row.names=FALSE, col.names=TRUE)
      }
    }
  }
}
