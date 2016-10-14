##
## Description: Summarizes the results of 3 steps decision tree stationarity tests
##
## Authors: Luis Gustavo Nardin
##          Vasile Alexandru Suchar
##
## Last Modification: 10/03/2016
##
library(caTools)
library(data.table)


#############
## PATHS
#############
baseDir <- "/data/projects/current/cmci/stationarity/data/ar1"
scriptDir <- paste0(baseDir, "/code")
inputDir <- paste0(baseDir, "/raw")
outputDir <- paste0(baseDir, "/summary")


#############
## FUNCTIONS
#############
commandArgs <- function(){ NULL }
source(paste0(scriptDir, "/functions.R"))


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
tests <- c("DT-KPSS-PSR-BOOTSTRAP+ERS")


##
## Stationary - ARMA(p,q)
##

## Window Size
windows <- sort(c(20, seq(30, TS, 10),
                  c(32, 64, 128, 256, 512, 1024, 2048, 4096)))

phis <- c(0, 0.99, 0.9, 0.8, 0.5, 0.2, -0.2, -0.5, -0.8, -0.9, -0.99, 1, -1)
thetas <- c(0)

for(error in ERROR_OPTION){
  summary <- NULL
  for(p in phis){
    for(q in thetas){
      
      data <- NULL
      for(tindex in 1:length(tests)){
        armaResults <- read.csv(paste0(inputDir, "/dt-3steps-arma-", error,
                                       "/arma-p", p, "-q", q, "-", tests[tindex],
                                       "-",error, ".csv"),
                                header=TRUE, sep=";")
        
        for(windex in 1:length(windows)){
          num <- length(armaResults[,windex][!is.na(armaResults[,windex])])
          calc <- sum(armaResults[,windex], na.rm=TRUE)
          if(num == 0){
            data <- rbind(data, cbind(tests[tindex], windows[windex],
                                      NA, 0, 0, N))
          } else {
            data <- rbind(data, cbind(tests[tindex], windows[windex],
                                      calc / num, calc, num - calc, N - num))
          }
        }
      }
      
      summary <- rbind(summary, cbind(data[,1], p, q,
                                      as.integer(as.character(data[,2])),
                                      as.numeric(as.character(data[,3])),
                                      as.integer(as.character(data[,4])),
                                      as.integer(as.character(data[,5])),
                                      as.integer(as.character(data[,6]))))
      
      gc(reset=TRUE)
    }
  }
  
  summary <- data.table(summary)
  names(summary) <- c("test", "p", "q", "window", "value", "stationary",
                      "non-stationary", "NA")
  
  write.table(summary, file=paste0(outputDir, "/dt-3steps-arma-", error, "-summary.csv"),
              append=FALSE, quote=FALSE, sep=";", row.names=FALSE,
              col.names=TRUE)
}


## Time
for(error in ERROR_OPTION){
  summary <- NULL
  for(p in phis){
    for(q in thetas){
      
      data <- NULL
      armaTime <- read.csv(paste0(inputDir, "/dt-3steps-arma-", error,
                                  "/armaTime-p", p, "-q", q, "-", error, ".csv"),
                           header=TRUE, sep=";")
      for(tindex in 1:length(tests)){
        for(windex in 1:length(windows)){
          data <- rbind(data, cbind(tests[tindex], windows[windex],
                                    armaTime[tindex, windex+1]))
        }
      }
      
      summary <- rbind(summary, cbind(data[,1], p, q,
                                      as.integer(as.character(data[,2])),
                                      as.numeric(as.character(data[,3]))))
    }
  }
  
  summary <- data.table(summary)
  names(summary) <- c("test", "p", "q", "window", "value")
  
  write.table(summary, file=paste0(outputDir, "/dt-3steps-arma-", error, "-time.csv"),
              append=FALSE, quote=FALSE, sep=";", row.names=FALSE,
              col.names=TRUE)
}


##
## Non-Stationary - TREND
##

## Window Size
windows <- c(32, 64, 128, 256, 512, 1024, 2048, 4096)

## Error
error <- ERROR_NORMAL

##
## M Trend Mean
## V Trend Variance
## AC Auto Correlation
##
types <- c("M", "V", "AC")

summary <- NULL
for(type in types){
  if (type == "M"){
    values <- c(0, 0.5, 1)
  } else if (type == "V"){
    values <- c(0, 1, 3, 5)
  } else if (type == "AC"){
    values <- c(0, 0.9, 0.45, -0.9, -0.45)
  }
  
  if((type == "M") | (type == "V")){
    for(tindex in 1:length(tests)){
      for(phi in 1:length(phis)){
        for(value in 1:length(values)){
          armaResults <- read.csv(paste0(inputDir,"/dt-3steps-trend/trend", type,
                                         values[value], "-phi", phis[phi], "-",
                                         tests[tindex], "-",error, ".csv"),
                                  header=TRUE, sep=";")
          
          data <- NULL
          for(windex in 1:length(windows)){
            num <- length(armaResults[,windex][!is.na(armaResults[,windex])])
            calc <- sum(armaResults[,windex], na.rm=TRUE)
            if(num == 0){
              data <- rbind(data, cbind(tests[tindex], windows[windex],
                                        NA, 0, 0, N))
            } else {
              data <- rbind(data, cbind(tests[tindex], windows[windex],
                                        calc / num, calc, num - calc, N - num))
            }
          }
          
          summary <- rbind(summary, cbind(data[,1], type, values[value], phis[phi],
                                          as.integer(as.character(data[,2])),
                                          as.numeric(as.character(data[,3])),
                                          as.integer(as.character(data[,4])),
                                          as.integer(as.character(data[,5])),
                                          as.integer(as.character(data[,6]))))
        }
      }
    }
  } else if(type == "AC"){
    for(tindex in 1:length(tests)){
      for(value in 1:length(values)){
        armaResults <- read.csv(paste0(inputDir,"/dt-3steps-trend/trend", type,
                                       values[value], "-", tests[tindex], "-",
                                       error, ".csv"),
                                header=TRUE, sep=";")
        data <- NULL
        for(windex in 1:length(windows)){
          num <- length(armaResults[,windex][!is.na(armaResults[,windex])])
          calc <- sum(armaResults[,windex], na.rm=TRUE)
          if(num == 0){
            data <- rbind(data, cbind(tests[tindex], windows[windex],
                                      NA, 0, 0, N))
          } else {
            data <- rbind(data, cbind(tests[tindex], windows[windex],
                                      calc / num, calc, num - calc, N - num))
          }
        }
        
        summary <- rbind(summary, cbind(data[,1], type, values[value], 0,
                                        as.integer(as.character(data[,2])),
                                        as.numeric(as.character(data[,3])),
                                        as.integer(as.character(data[,4])),
                                        as.integer(as.character(data[,5])),
                                        as.integer(as.character(data[,6]))))
      }
    }
  }
}

summary <- data.table(summary)
names(summary) <- c("test", "type", "param", "phi", "window", "value",
                    "stationary", "non-stationary", "NA")

write.table(summary, file=paste0(outputDir, "/dt-3steps-trend-summary.csv"),
            append=FALSE, quote=FALSE, sep=";", row.names=FALSE,
            col.names=TRUE)


## Time
summary <- NULL
for(type in types){
  if (type == "M"){
    values <- c(0, 0.5, 1)
  } else if (type == "V"){
    values <- c(0, 1, 3, 5)
  } else if (type == "AC"){
    values <- c(0, 0.9, 0.45, -0.9, -0.45)
  }
  
  if((type == "M") | (type == "V")){
    for(phi in 1:length(phis)){
      for(value in 1:length(values)){
        armaTime <- read.csv(paste0(inputDir,"/dt-3steps-trend/trend", type,
                                    "Time", values[value], "-phi", phis[phi], "-",
                                    error, ".csv"),
                             header=TRUE, sep=";")
        data <- NULL
        for(tindex in 1:length(tests)){
          for(windex in 1:length(windows)){
            data <- rbind(data, cbind(tests[tindex], windows[windex],
                                      armaTime[tindex,windex+1]))
          }
        }
        
        summary <- rbind(summary, cbind(data[,1], type, values[value], phis[phi],
                                        as.integer(as.character(data[,2])),
                                        as.numeric(as.character(data[,3]))))
      }
    }
  } else if(type == "AC"){
    for(value in 1:length(values)){
      armaTime <- read.csv(paste0(inputDir,"/dt-3steps-trend/trend", type,
                                  "Time", values[value], "-", error, ".csv"),
                           header=TRUE, sep=";")
      data <- NULL
      for(tindex in 1:length(tests)){
        for(windex in 1:length(windows)){
          data <- rbind(data, cbind(tests[tindex], windows[windex],
                                    armaTime[tindex, windex+1]))
        }
      }
      
      summary <- rbind(summary, cbind(data[,1], type, values[value], 0,
                                      as.integer(as.character(data[,2])),
                                      as.numeric(as.character(data[,3]))))
    }
  }
}

summary <- data.table(summary)
names(summary) <- c("test", "type", "param", "phi", "window", "value")

write.table(summary, file=paste0(outputDir, "/dt-3steps-trend-time.csv"),
            append=FALSE, quote=FALSE, sep=";", row.names=FALSE,
            col.names=TRUE)


##
## Non-Stationary - BREAK
##

## Window Size
windows <- c(32, 64, 128, 256, 512, 1024, 2048, 4096)

## Error
error <- ERROR_NORMAL

##
## M Trend Mean
## V Trend Variance
## AC Auto Correlation
##
types <- c("M", "V", "AC")

summary <- NULL
for(type in types){
  if (type == "M"){
    values <- c(0, 1, 2, 3, 4)
  } else if (type == "V"){
    values <- c(0, 0.5, 1, 2, 3, 4)
  } else if (type == "AC"){
    values <- c(0.9, 0.45, 0, -0.45, -0.9)
  }
  
  if((type == "M") || (type == "V")){
    for(tindex in 1:length(tests)){
      for(phi in 1:length(phis)){
        for(value in 1:length(values)){
          armaResults <- read.csv(paste0(inputDir, "/dt-3steps-break/break",
                                         type, values[value], "-phi", phis[phi],
                                         "-", tests[tindex], "-", error, ".csv"),
                                  header=TRUE, sep=";")
          
          data <- NULL
          for(windex in 1:length(windows)){
            num <- length(armaResults[,windex][!is.na(armaResults[,windex])])
            calc <- sum(armaResults[,windex], na.rm=TRUE)
            if(num == 0){
              data <- rbind(data, cbind(tests[tindex], windows[windex],
                                        NA, 0, 0, N))
            } else {
              data <- rbind(data, cbind(tests[tindex], windows[windex],
                                        calc / num, calc, num - calc, N - num))
            }
          }
          
          summary <- rbind(summary, data.table(data[,1], type,
                                               as.character(values[value]),
                                               phis[phi],
                                               as.integer(as.character(data[,2])),
                                               as.numeric(as.character(data[,3])),
                                               as.integer(as.character(data[,4])),
                                               as.integer(as.character(data[,5])),
                                               as.integer(as.character(data[,6]))))
        }
      }
    }
  } else if(type == "AC"){
    for(tindex in 1:length(tests)){
      items <- combs(values, 2)
      
      for(item in 1:nrow(items)){
        trendAC1 <- items[item, 1]
        trendAC2 <- items[item, 2]
        
        armaResults <- read.csv(paste0(inputDir, "/dt-3steps-break/break",
                                       type, trendAC1, "-", trendAC2, "-",
                                       tests[tindex], "-", error, ".csv"),
                                header=TRUE, sep=";")
        
        data <- NULL
        for(windex in 1:length(windows)){
          num <- length(armaResults[,windex][!is.na(armaResults[,windex])])
          calc <- sum(armaResults[,windex], na.rm=TRUE)
          if(num == 0){
            data <- rbind(data, cbind(tests[tindex], windows[windex],
                                      NA, 0, 0, N))
          } else {
            data <- rbind(data, cbind(tests[tindex], windows[windex],
                                      calc / num, calc, num - calc, N - num))
          }
        }
        
        summary <- rbind(summary, data.table(data[,1], type,
                                        paste0(trendAC1,"-",trendAC2),
                                        0,
                                        as.integer(as.character(data[,2])),
                                        as.numeric(as.character(data[,3])),
                                        as.integer(as.character(data[,4])),
                                        as.integer(as.character(data[,5])),
                                        as.integer(as.character(data[,6]))))
      }
    }
  }
}

summary <- data.table(summary)
names(summary) <- c("test", "type", "param", "phi", "window", "value",
                    "stationary", "non-stationary", "NA")

write.table(summary, file=paste0(outputDir, "/dt-3steps-break-summary.csv"),
            append=FALSE, quote=FALSE, sep=";", row.names=FALSE,
            col.names=TRUE)


## Time
summary <- NULL
for(type in types){
  if (type == "M"){
    values <- c(0, 1, 2, 3, 4)
  } else if (type == "V"){
    values <- c(0, 0.5, 1, 2, 3, 4)
  } else if (type == "AC"){
    values <- c(0.9, 0.45, 0, -0.45, -0.9)
  }
  
  
  if ((type == "M") || (type == "V")){
    for(phi in 1:length(phis)){
      for(value in 1:length(values)){
        armaTime <- read.csv(paste0(inputDir, "/dt-3steps-break/break", type,
                                    "Time", values[value], "-phi", phis[phi],
                                    "-", error, ".csv"),
                             header=TRUE, sep=";")
        
        data <- NULL
        for(tindex in 1:length(tests)){
          for(windex in 1:length(windows)){
            data <- rbind(data, cbind(tests[tindex], windows[windex],
                                      armaTime[tindex, windex+1]))
          }
        }
        
        summary <- rbind(summary, cbind(data[,1], type, values[value], phis[phi],
                                        as.integer(as.character(data[,2])),
                                        as.numeric(as.character(data[,3]))))
      }
    }
  } else if (type == "AC"){
    items <- combs(values, 2)
    
    for(item in 1:nrow(items)){
      trendAC1 <- items[item, 1]
      trendAC2 <- items[item, 2]
      
      armaTime <- read.csv(paste0(inputDir, "/dt-3steps-break/break", type,
                                  "Time", trendAC1, "-", trendAC2, "-", error,
                                  ".csv"),
                           header=TRUE, sep=";")
      
      data <- NULL
      for(tindex in 1:length(tests)){
        for(windex in 1:length(windows)){
          data <- rbind(data, cbind(tests[tindex], windows[windex],
                                    armaTime[tindex, windex+1]))
        }
      }
      
      summary <- rbind(summary, cbind(data[,1], type,
                                      paste0(trendAC1,"-",trendAC2), 0,
                                      as.integer(as.character(data[,2])),
                                      as.numeric(as.character(data[,3]))))
    }
  }
}

summary <- data.table(summary)
names(summary) <- c("test", "type", "param", "phi", "window", "value")

write.table(summary, file=paste0(outputDir, "/dt-3steps-break-time.csv"),
            append=FALSE, quote=FALSE, sep=";", row.names=FALSE,
            col.names=TRUE)
