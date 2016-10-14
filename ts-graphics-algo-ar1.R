##
## Description: Plot the results of the stationarity algorithm
##
## Authors: Luis Gustavo Nardin
##          Vasile Alexandru Suchar
##
## Last Modification: 10/10/2016
##
library(caTools)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)


#############
## PATHS
#############
baseDir <- "/data/Box/Stationarity"
scriptDir <- paste0(baseDir, "/code")
inputDir <- paste0(baseDir, "/data/ar1/raw")
outputDir <- paste0(baseDir, "/data/ar1/figures")


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
tests <- c("ALGO-ERS")


###################
## GRAPHICS
###################

##
## Stationary - ARMA(p,q)
##

## Window Size
windows <- sort(c(20, seq(30, TS, 10),
                  c(32, 64, 128, 256, 512, 1024, 2048, 4096)))

phis <- c(0, 0.99, 0.9, 0.8, 0.5, 0.2, -0.2, -0.5, -0.8, -0.9, -0.99, 1, -1)
thetas <- c(0)

for(error in ERROR_OPTION){
  for(p in phis){
    for(q in thetas){
      
      data <- NULL
      for(test in 1:length(tests)){
        armaResults <- read.csv(paste0(inputDir, "/algo-arma-", error,
                                       "/arma-p", p, "-q", q, "-", tests[test],
                                       "-", error, ".csv"),
                                header=TRUE, sep=";")
        
        for(windex in 1:length(windows)){
          data <- rbind(data, cbind(tests[test], windows[windex], sum(armaResults[,windex]) / N))
        }
      }
      
      data <- data.table(data[,1],
                         as.integer(as.character(data[,2])),
                         as.numeric(as.character(data[,3])))
      names(data) <- c("test", "window", "value")
      
      pl <- ggplot(data[which(value > 0)], aes(x=window, y=value*100, group=test)) +
        xlab("Window Size") +
        ggtitle(substitute(paste("AR(1)   ", phi[1]," = ", p), list(p=p))) +
        geom_line(aes(color=test), size=1.5) +
        scale_y_continuous(name="% Detected Stationarity",
                           breaks=c(0, 20, 40, 60, 80, 100),
                           labels=c("0%", "20%", "40%", "60%", "80%", "100%"),
                           limits=c(0, 100)) +
        guides(color=guide_legend(override.aes=list(size=2))) +
        theme(axis.title.x=element_text(color='black', size=28, face='bold'),
              axis.title.y=element_text(color='black', size=28, face='bold'),
              axis.text.x=element_text(color='black', size=18, face='bold'),
              axis.text.y=element_text(color='black', size=18, face='bold'),
              axis.line.x=element_line(color='black', size=1.5, linetype='solid'),
              axis.line.y=element_line(color='black', size=1.5, linetype='solid'),
              panel.background=element_rect(fill="transparent", color=NA),
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title=element_text(color="black", size=24, face="bold"),
              plot.margin=unit(c(1, 8, 0.5, 0.5), "lines"),
              legend.position=c(1,1), 
              legend.justification=c(0, 1),
              legend.title=element_blank(),
              legend.text=element_text(color="black", size=18, face="bold"),
              legend.key=element_rect(fill="white"),
              legend.key.height=unit(2,"line"))
      
      ggsave(paste0(outputDir,"/algo-arma-", error, "/stationary-ar1-" ,
                    p, "-", error, ".png"), pl, width=20, height=9.5)
    }
  }
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

for(type in types){
  if (type == "M"){
    values <- c(0, 0.5, 1, 2, 3, 4)
    
    distortion <- c(rep(0, length(windows)),
                    rep(0.01, length(windows)),
                    rep(0.02, length(windows)))
    
  } else if (type == "V"){
    values <- c(0, 0.5, 1, 2, 3, 4)
    
    distortion <- c(rep(0, length(windows)),
                    rep(0.01, length(windows)),
                    rep(0.02, length(windows)),
                    rep(0.03, length(windows)))
    
  } else if (type == "AC"){
    values <- c(0, 0.9, 0.45, -0.9, -0.45)
    
    distortion <- c(rep(0, length(windows)),
                    rep(0.01, length(windows)),
                    rep(0.02, length(windows)),
                    rep(0.03, length(windows)),
                    rep(0.04, length(windows)))
  }
  
  if((type == "M") | (type == "V")){
    for(test in 1:length(tests)){
      for(phi in 1:length(phis)){
        data <- NULL
        for(value in 1:length(values)){
          armaResults <- read.csv(paste0(inputDir, "/algo-trend/trend", type,
                                         values[value], "-phi", phis[phi], "-",
                                         tests[test], "-", error, ".csv"),
                                  header=TRUE, sep=";")
          
          for(windex in 1:length(windows)){
            data <- rbind(data, cbind(values[value], windows[windex],
                                      sum(armaResults[,windex]) / N))
          }
        }
        
        data <- data.table(data[,1],
                           as.integer(as.character(data[,2])),
                           as.numeric(as.character(data[,3]))-distortion)
        names(data) <- c("m", "window", "value")
        
        pl <- ggplot(data[which(window <= 4096)], aes(x=as.numeric(as.character(window)),
                         y=as.numeric(as.character(value))*100,
                         group=as.factor(as.numeric(as.character(m))))) +
          geom_line(aes(color=as.factor(as.numeric(as.character(m)))),
                        size=2) +
          xlab("Window Size") +
          ggtitle(substitute(paste(t, "     ", phi[1], " = ", p), list(t=tests[test], p=phis[phi]))) +
          scale_y_continuous(name="% Detected Stationarity",
                             breaks=c(0, 20, 40, 60, 80, 100),
                             labels=c("0%", "20%", "40%", "60%", "80%", "100%"),
                             limits=c(-10, 100)) +
          scale_color_discrete(name="Trend\nMean") +
          guides(color=guide_legend(override.aes=list(size=2))) +
          theme(axis.title.x=element_text(color='black', size=28, face='bold'),
                axis.title.y=element_text(color='black', size=28, face='bold'),
                axis.text.x=element_text(color='black', size=18, face='bold'),
                axis.text.y=element_text(color='black', size=18, face='bold'),
                axis.line.x=element_line(color='black', size=1.5, linetype='solid'),
                axis.line.y=element_line(color='black', size=1.5, linetype='solid'),
                panel.background=element_rect(fill="transparent", color=NA),
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                plot.title=element_text(color="black", size=24, face="bold"),
                plot.margin=unit(c(1, 15, 0.5, 0.5), "lines"),
                legend.position="bottom",
                legend.justification=c(0, 1),
                legend.title=element_text(color="black", size=18, face="bold"),
                legend.text=element_text(color="black", size=18, face="bold"),
                legend.key=element_rect(fill="white"),
                legend.key.height=unit(2,"line")
                )
        
        ggsave(paste0(outputDir, "/algo-trend/trend", type, "-phi",
                      phis[phi], "-", tests[test], ".png"),
               pl, width=9.86, height=5.43)
      }
    }
  } else if (type == "AC"){
    for(test in 1:length(tests)){
      data <- NULL
      for(value in 1:length(values)){
        armaResults <- read.csv(paste0(inputDir, "/algo-trend/trend", type,
                                       values[value], "-", tests[test], "-",
                                       error, ".csv"),
                                header=TRUE, sep=";")
        
        for(windex in 1:length(windows)){
          data <- rbind(data, cbind(values[value], windows[windex],
                                    sum(armaResults[,windex]) / N))
        }
      }
      
      data <- data.table(data[,1],
                         as.integer(as.character(data[,2])),
                         as.numeric(as.character(data[,3]))-distortion)
      names(data) <- c("m", "window", "value")
      
      pl <- ggplot(data[which(window <= 4096)], aes(x=as.numeric(as.character(window)),
                                                    y=as.numeric(as.character(value))*100,
                                                    group=as.factor(as.numeric(as.character(m))))) +
        geom_line(aes(color=as.factor(as.numeric(as.character(m)))),
                  size=2) +
        xlab("Window Size") +
        ggtitle(substitute(paste(t, "     ", p),
                           list(t=tests[test], p=values[value]))) +
        scale_y_continuous(name="% Detected Stationarity",
                           breaks=c(0, 20, 40, 60, 80, 100),
                           labels=c("0%", "20%", "40%", "60%", "80%", "100%"),
                           limits=c(-10, 100)) +
        scale_color_discrete(name="Trend\nMean") +
        guides(color=guide_legend(override.aes=list(size=2))) +
        theme(axis.title.x=element_text(color='black', size=28, face='bold'),
              axis.title.y=element_text(color='black', size=28, face='bold'),
              axis.text.x=element_text(color='black', size=18, face='bold'),
              axis.text.y=element_text(color='black', size=18, face='bold'),
              axis.line.x=element_line(color='black', size=1.5, linetype='solid'),
              axis.line.y=element_line(color='black', size=1.5, linetype='solid'),
              panel.background=element_rect(fill="transparent", color=NA),
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              plot.title=element_text(color="black", size=24, face="bold"),
              plot.margin=unit(c(1, 20, 0.5, 0.5), "lines"),
              legend.position="bottom",
              legend.justification=c(0, 1),
              legend.title=element_text(color="black", size=18, face="bold"),
              legend.text=element_text(color="black", size=18, face="bold"),
              legend.key=element_rect(fill="white"),
              legend.key.height=unit(2,"line")
        )
      
      ggsave(paste0(outputDir, "/algo-trend/trend", type, "-",
                    tests[test], ".png"),
             pl, width=9.86, height=5.43)
    }
  }
}


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

for(type in types){
	if (type == "M"){
		values <- c(0, 0.5, 1, 2, 3, 4)
		
		distortion <- c(rep(0, length(windows)),
				rep(0.01, length(windows)),
				rep(0.02, length(windows)),
				rep(0.03, length(windows)),
				rep(0.04, length(windows)))
	} else if (type == "V"){
		values <- c(0, 0.5, 1, 2, 3, 4)
		
		distortion <- c(rep(0, length(windows)),
				rep(0.01, length(windows)),
				rep(0.02, length(windows)),
				rep(0.03, length(windows)),
				rep(0.04, length(windows)),
				rep(0.05, length(windows)))
	} else if (type == "AC"){
		values <- c(0.9, 0.45, 0, -0.45, -0.9)
		
		distortion <- c(rep(0, length(windows)),
				rep(0.01, length(windows)),
				rep(0.02, length(windows)),
				rep(0.03, length(windows)),
				rep(0.04, length(windows)))
	}
	
	if ((type == "M") || (type == "V")){
		for(tindex in 1:length(tests)){
			for(phi in 1:length(phis)){
				data <- NULL
				for(value in 1:length(values)){
					armaResults <- read.csv(paste0(inputDir, "/algo-break/break", type,
									values[value], "-phi", phis[phi], "-",
									tests[tindex], "-", error, ".csv"),
							header=TRUE, sep=";")
					
					for(windex in 1:length(windows)){
						data <- rbind(data, cbind(values[value], windows[windex],
										sum(armaResults[,windex]) / N))
					}
				}
				
				data <- data.table(data[,1],
						as.integer(as.character(data[,2])),
						as.numeric(as.character(data[,3]))-distortion)
				names(data) <- c("m", "window", "value")
				
				pl <- ggplot(data, aes(x=as.numeric(as.character(window)),
										y=as.numeric(as.character(value)) * 100,
										group=as.factor(as.numeric(as.character(m))))) +
						geom_line(aes(color=as.factor(as.numeric(as.character(m)))),
								size=2) +
						xlab("Window Size") +
						ggtitle(substitute(paste(t, "     ", phi[1], " = ", p),
										list(t=tests[tindex], p=phis[phi]))) +
						scale_y_continuous(name="% Detected Stationarity",
								breaks=c(0, 20, 40, 60, 80, 100),
								labels=c("0%", "20%", "40%", "60%", "80%", "100%")) +
						scale_color_discrete(name="Structural\nBreak") +
						guides(color=guide_legend(override.aes=list(size=2))) +
						theme(axis.title.x=element_text(color='black', size=28, face='bold'),
								axis.title.y=element_text(color='black', size=28, face='bold'),
								axis.text.x=element_text(color='black', size=18, face='bold'),
								axis.text.y=element_text(color='black', size=18, face='bold'),
								axis.line.x=element_line(color='black', size=1.5, linetype='solid'),
								axis.line.y=element_line(color='black', size=1.5, linetype='solid'),
								panel.background=element_rect(fill="transparent", color=NA),
								panel.grid.minor=element_blank(),
								panel.grid.major=element_blank(),
								plot.title=element_text(color="black", size=24, face="bold"),
								plot.margin=unit(c(1, 8, 0.5, 0.5), "lines"),
								legend.position="bottom",
								legend.justification=c(0, 1),
								legend.title=element_text(color="black", size=18, face="bold"),
								legend.text=element_text(color="black", size=18, face="bold"),
								legend.key=element_rect(fill="white"),
								legend.key.height=unit(2,"line")
						)
				
				ggsave(paste0(outputDir, "/algo-break/break", type, "-phi", phis[phi], "-",
								tests[tindex], ".png"), pl, width=9.86, height=5.43)
			}
		}
	} else if (type == "AC"){
		for(tindex in 1:length(tests)){
			items <- combs(values, 2)
			
			data <- NULL
			for(item in 1:nrow(items)){
				trendAC1 <- items[item, 1]
				trendAC2 <- items[item, 2]
				
				armaResults <- read.csv(paste0(inputDir, "/algo-break/break", type,
								trendAC1, "-", trendAC2, "-",
								tests[tindex], "-", error, ".csv"),
						header=TRUE, sep=";")
				
				for(windex in 1:length(windows)){
					data <- rbind(data, cbind(paste0(trendAC1,"-",trendAC2), windows[windex],
									sum(armaResults[,windex]) / N))
				}
			}
			
			data <- data.table(data[,1],
					as.integer(as.character(data[,2])),
					as.numeric(as.character(data[,3]))-distortion)
			names(data) <- c("m", "window", "value")
			
			pl <- ggplot(data, aes(x=as.numeric(as.character(window)),
									y=as.numeric(as.character(value)) * 100,
									group=as.factor(as.character(m)))) +
					geom_line(aes(color=as.factor(as.character(m))),
							size=2) +
					xlab("Window Size") +
					ggtitle(substitute(paste(t),
									list(t=tests[tindex]))) +
					scale_y_continuous(name="% Detected Stationarity",
							breaks=c(0, 20, 40, 60, 80, 100),
							labels=c("0%", "20%", "40%", "60%", "80%", "100%")) +
					scale_color_discrete(name="Structural\nBreak") +
					guides(color=guide_legend(override.aes=list(size=2))) +
					theme(axis.title.x=element_text(color='black', size=28, face='bold'),
							axis.title.y=element_text(color='black', size=28, face='bold'),
							axis.text.x=element_text(color='black', size=18, face='bold'),
							axis.text.y=element_text(color='black', size=18, face='bold'),
							axis.line.x=element_line(color='black', size=1.5, linetype='solid'),
							axis.line.y=element_line(color='black', size=1.5, linetype='solid'),
							panel.background=element_rect(fill="transparent", color=NA),
							panel.grid.minor=element_blank(),
							panel.grid.major=element_blank(),
							plot.title=element_text(color="black", size=24, face="bold"),
							plot.margin=unit(c(1, 8, 0.5, 0.5), "lines"),
							legend.position="bottom",
							legend.justification=c(0, 1),
							legend.title=element_text(color="black", size=18, face="bold"),
							legend.text=element_text(color="black", size=18, face="bold"),
							legend.key=element_rect(fill="white"),
							legend.key.height=unit(2,"line")
					)
			
			ggsave(paste0(outputDir, "/algo-break/break", type, "-", tests[tindex], ".png"),
					pl, width=9.86, height=5.43)
		}
	}
}
