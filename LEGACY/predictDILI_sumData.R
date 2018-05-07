## Predict DILI compounds

## Settings
	options(stringsAsFactors = FALSE)
	
	library(data.table)
	library(preprocessCore)

	library(parallel)

	library(neuralnet)

	library(ggplot2)

	inputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'
	inputFiles <- list.files(inputDir, pattern = '.Rdata')	

## Functions
	calculateFCvsDMSO <- function(sumData) {
  		if(FALSE) sumData <- df <- CHOP[[1]]
  		
  		# Internal parser function
		.parser <- function(df, catProgress = FALSE, ...) {
			# Treatment names
			if(catProgress) cat('Parseing treatments...\n')
			df$treatment <- gsub('_', '', toupper(gsub(' ', '', df$treatment)))
				df$treatment[df$treatment == 'BUTHIONINESULFOXAMINE'] <- 'BUTHIONINESULFOXIMINE'
				df$treatment[df$treatment == 'STAUROSPORIN'] <- 'STAUSPORIN'
				df$treatment[df$treatment == 'SODIUMARSENITE'] <- 'SODIUM-ARSENITE'
			if(catProgress) cat('1 of 5 steps done.\n')

			# Standardize plateID in sumData
			if(catProgress) cat('Parseing plateIDs...\n')
			tmp <- do.call(rbind, strsplit(toupper(gsub('_wells', '', df$plateID)), '_'))
				repInd  <- grep('REP', tmp[1, ])
				cmaxInd <- grep('CMAX', tmp[1, ])
				timeInd <- grep('H', tmp[1, ])
			df$plateID <- paste(tmp[,repInd], tmp[, cmaxInd], tmp[, timeInd], sep = '_') 
			if(catProgress) cat('2 of 5 steps done.\n')

			# Create uniqueID in sumData
			if(catProgress) cat('Parseing uniqueIDs...\n')
			df$uniqueID <- paste(df$treatment, df$timeID, df$plateID, sep = '_')
			if(catProgress) cat('3 of 5 steps done.\n')

			# Split uniqueID and create CMAX columns
			if(catProgress) cat('Parseing CMAX column...\n')
			uniqueID_splitted <- do.call(rbind, strsplit(df$uniqueID, '_'))
			df$CMAX <- uniqueID_splitted[ ,grep('CMAX', toupper(uniqueID_splitted)[1, ])]
			df$CMAX_numeric <- as.numeric(gsub('CMAX', '', df$CMAX))
			if(catProgress) cat('4 of 5 steps done.\n')

			# Create columns to match against
			if(catProgress) cat('Parseing matching columns...\n')
			df$matchDose <- paste(df$treatment, df$dose_uM, df$CMAX, sep = '_')
			df$matchDoseTime <- paste(df$treatment, df$dose_uM, df$CMAX, df$timeID,sep = '_')
			if(catProgress) cat('5 of 5 steps done.\n\n\n')

			# return
			return(df)
		}

		# do it to it!
		sumData <- as.data.frame(.parser(sumData))
  			colnames(sumData)[grep('GFP', colnames(sumData))] <- 'meanIntensity_GFP'
  		
  		plateIDs <- as.character(unique(sumData$plateID))
  		meanDMSO <- sapply(plateIDs, function(.pID) mean(sumData$meanIntensity_GFP[sumData$plateID == .pID & sumData$treatment == 'DMSO']))
  		
  		sumData$meanPlateDMSO <- NA
  		for(.pID in plateIDs) sumData$meanPlateDMSO[sumData$plateID == .pID] <- meanDMSO[.pID]

  		sumData$GFP_FC <- sumData$meanIntensity_GFP / sumData$meanPlateDMSO

  		rslt <- sumData[sumData$treatment != 'DMSO', ]
  		rslt <- rslt[order(rslt$uniqueID), ]

  		rslt
  	}		


## Load data
  # DILI annotation
  	DILI <- read.delim('/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/DILI_annotation.txt')

  # BIP
	load(paste0(inputDir, inputFiles[1]))	
	BIP <- lapply(outputList, function(x) x$sumData)
  # CHOP
  	load(paste0(inputDir, inputFiles[2]))	
	CHOP <- lapply(outputList, function(x) x$sumData)
  # P21	
  	load(paste0(inputDir, inputFiles[5]))	
	P21 <- lapply(outputList, function(x) x$sumData)

	rm(outputList)

## Parse data
  # DILI annotation
  	DILI$binaryDILI <- DILI$doseDILI <- as.numeric(grepl('Most', DILI$updatedDILI) | grepl('Less', DILI$updatedDILI))
  	DILI$doseDILI[grepl('Most', DILI$updatedDILI)] <- 2 

  	DILI$treatment <- toupper(DILI$treatment)
  	  DILI$treatment[DILI$treatment == "BUTHIONINE SULFOXAMINE"] <- "BUTHIONINESULFOXIMINE"
  	  DILI$treatment[DILI$treatment == "CARBAMAZEPINE"] <- "CARBAMAZAPINE"
  	  DILI$treatment[DILI$treatment == "CYCLOHEXIMIDE"] <- "CYCLOHEXAMINE"
  	  DILI$treatment[DILI$treatment == "CYCLOSPORIN A"] <- "CYCLOSPORINA"
  	  DILI$treatment[DILI$treatment == "ETOPOSIDE"] <- "ETOPOSIDE "
  	  DILI$treatment[DILI$treatment == "SODIUM ARSENITE"] <- "SODIUM-ARSENITE"
  	  DILI$treatment[DILI$treatment == "STAUROSPORIN"] <- "STAUSPORIN"

  	  DILI$treatment <- gsub(' ', '', DILI$treatment)
	
  # sumData
	pBIP  <- lapply(BIP, calculateFCvsDMSO)	
	pCHOP <- lapply(CHOP, calculateFCvsDMSO)	
	pP21  <- lapply(P21, calculateFCvsDMSO)	

  # Create 24hr, 48hr & 72hr short-format data.frames
  	overlappingIDs <- table(c(unique(sapply(pBIP, function(x) x$uniqueID)), unique(sapply(pCHOP, function(x) x$uniqueID)), unique(sapply(pP21, function(x) x$uniqueID))))
  	overlappingIDs <- names(which(overlappingIDs == 3))

  	treatments <- unique(do.call(rbind, strsplit(overlappingIDs, '_'))[,1])

  	mergeReplicates <- function(parsedData, parsedDataName, treatment, timeID, exceptions = c('ASPIRIN', 'DMEM', 'METFORMIN')) {
  		if(FALSE) {
  			parsedData = pCHOP
  			timeID = 24
  			treatment = 'METFORMIN'
        parsedDataName = 'CHOP'
  		}

  		if(!treatment %in% DILI$treatment) return(rep(NA, 6))

  		tmp <- lapply(parsedData, function(x) x[x$treatment == treatment &x$timeID == timeID, ])
  		tmp <- lapply(tmp, function(x) x[order(x$plateID), ])

      
      if(!(length(unique(sapply(tmp, nrow))) == 1)) {
  			cat(treatment,'\n')
  			stop('Bleh.\n')
  		}

  		merged <- rowMeans(sapply(tmp, function(x) x$GFP_FC))
  		names(merged) <- paste0(parsedDataName, '_', tmp[[1]]$CMAX)

      if(treatment %in% exceptions) {
        merged <- merged[seq(1, length(merged), 2)]
      }

  		merged
  	}


  # Reporter & timepoint specific
    assemble_Rprtr_24 <- function() {
      datBIP_24 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pBIP, parsedDataName = 'BIP', treatment = trt, timeID = 24)))
      datCHOP_24 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pCHOP, parsedDataName = 'CHOP', treatment = trt, timeID = 24)))
      datP21_24 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pP21, parsedDataName = 'P21', treatment = trt, timeID = 24)))
        rownames(datBIP_24) <- rownames(datCHOP_24) <- rownames(datP21_24) <- treatments

      datBIP_24 <- data.frame(datBIP_24[which(rowSums(is.na(datBIP_24)) == 0), ])
      datCHOP_24 <- data.frame(datCHOP_24[which(rowSums(is.na(datCHOP_24)) == 0), ])
      datP21_24 <- data.frame(datP21_24[which(rowSums(is.na(datP21_24)) == 0), ])

      datBIP_24$DILI <- sapply(rownames(datBIP_24), function(x) DILI$binaryDILI[DILI$treatment == x])
      datCHOP_24$DILI <- sapply(rownames(datCHOP_24), function(x) DILI$binaryDILI[DILI$treatment == x])
      datP21_24$DILI <- sapply(rownames(datP21_24), function(x) DILI$binaryDILI[DILI$treatment == x])

      list(datBIP_24 = datBIP_24, datCHOP_24 = datCHOP_24, datP21_24 = datP21_24)
    }

    assemble_Rprtr_48 <- function() {
      datBIP_48 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pBIP, parsedDataName = 'BIP', treatment = trt, timeID = 48)))
      datCHOP_48 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pCHOP, parsedDataName = 'CHOP', treatment = trt, timeID = 48)))
      datP21_48 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pP21, parsedDataName = 'P21', treatment = trt, timeID = 48)))
        rownames(datBIP_48) <- rownames(datCHOP_48) <- rownames(datP21_48) <- treatments

      datBIP_48 <- data.frame(datBIP_48[which(rowSums(is.na(datBIP_48)) == 0), ])
      datCHOP_48 <- data.frame(datCHOP_48[which(rowSums(is.na(datCHOP_48)) == 0), ])
      datP21_48 <- data.frame(datP21_48[which(rowSums(is.na(datP21_48)) == 0), ])

      datBIP_48$DILI <- sapply(rownames(datBIP_48), function(x) DILI$binaryDILI[DILI$treatment == x])
      datCHOP_48$DILI <- sapply(rownames(datCHOP_48), function(x) DILI$binaryDILI[DILI$treatment == x])
      datP21_48$DILI <- sapply(rownames(datP21_48), function(x) DILI$binaryDILI[DILI$treatment == x])


      list(datBIP_48 = datBIP_48, datCHOP_48 = datCHOP_48, datP21_48 = datP21_48)
    }

    assemble_Rprtr_72 <- function() {
      datBIP_72 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pBIP, parsedDataName = 'BIP', treatment = trt, timeID = 72)))
      datCHOP_72 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pCHOP, parsedDataName = 'CHOP', treatment = trt, timeID = 72)))
      datP21_72 <- do.call(rbind, lapply(treatments, function(trt) mergeReplicates(parsedData = pP21, parsedDataName = 'P21', treatment = trt, timeID = 72)))
        rownames(datBIP_72) <- rownames(datCHOP_72) <- rownames(datP21_72) <- treatments

      datBIP_72 <- data.frame(datBIP_72[which(rowSums(is.na(datBIP_72)) == 0), ])
      datCHOP_72 <- data.frame(datCHOP_72[which(rowSums(is.na(datCHOP_72)) == 0), ])
      datP21_72 <- data.frame(datP21_72[which(rowSums(is.na(datP21_72)) == 0), ])

      datBIP_72$DILI <- sapply(rownames(datBIP_72), function(x) DILI$binaryDILI[DILI$treatment == x])
      datCHOP_72$DILI <- sapply(rownames(datCHOP_72), function(x) DILI$binaryDILI[DILI$treatment == x])
      datP21_72$DILI <- sapply(rownames(datP21_72), function(x) DILI$binaryDILI[DILI$treatment == x])


      list(datBIP_72 = datBIP_72, datCHOP_72 = datCHOP_72, datP21_72 = datP21_72)
    }
   

  # Timepoints combined
    assemble_24 <- function() {
        mBIP_24 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pBIP, parsedDataName = 'BIP', treatment = trt, timeID = 24))
        mCHOP_24 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pCHOP, parsedDataName = 'CHOP', treatment = trt, timeID = 24))
        mP21_24 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pP21, parsedDataName = 'P21', treatment = trt, timeID = 24))
          names(mBIP_24) <- names(mCHOP_24) <- names(mP21_24) <- treatments

        dat24 <- data.frame(do.call(rbind, mBIP_24), do.call(rbind, mCHOP_24), do.call(rbind, mP21_24))
        dat24 <- dat24[which(rowSums(is.na(dat24)) == 0), ] 

        dat24$DILI <- sapply(rownames(dat24), function(x) DILI$binaryDILI[DILI$treatment == x]) 

        dat24
    }
  
  	assemble_48 <- function() {
        mBIP_48 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pBIP, parsedDataName = 'BIP', treatment = trt, timeID = 48))
        mCHOP_48 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pCHOP, parsedDataName = 'CHOP', treatment = trt, timeID = 48))
        mP21_48 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pP21, parsedDataName = 'P21', treatment = trt, timeID = 48))
          names(mBIP_48) <- names(mCHOP_48) <- names(mP21_48) <- treatments

        dat48 <- data.frame(do.call(rbind, mBIP_48), do.call(rbind, mCHOP_48), do.call(rbind, mP21_48))
        dat48 <- dat48[which(rowSums(is.na(dat48)) == 0), ]  

        dat48$DILI <- sapply(rownames(dat48), function(x) DILI$binaryDILI[DILI$treatment == x]) 

        dat48
    }

    assemble_72 <- function() {
        mBIP_72 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pBIP, parsedDataName = 'BIP', treatment = trt, timeID = 72))
        mCHOP_72 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pCHOP, parsedDataName = 'CHOP', treatment = trt, timeID = 72))
        mP21_72 <- lapply(treatments, function(trt) mergeReplicates(parsedData = pP21, parsedDataName = 'P21', treatment = trt, timeID = 72))
          names(mBIP_72) <- names(mCHOP_72) <- names(mP21_72) <- treatments

        dat72 <- data.frame(do.call(rbind, mBIP_72), do.call(rbind, mCHOP_72), do.call(rbind, mP21_72))
        dat72 <- dat72[which(rowSums(is.na(dat72)) == 0), ]  

        dat72$DILI <- sapply(rownames(dat72), function(x) DILI$binaryDILI[DILI$treatment == x]) 


        dat72
    }


  # Reporters combined
    assemble_BIP <- function() {
      datBIP  <- data.frame(do.call(rbind, mBIP_24), do.call(rbind, mBIP_48), do.call(rbind, mBIP_72))
      datBIP <- datBIP[which(rowSums(is.na(datBIP)) == 0), ]  
        colnames(datBIP) <- paste0(sapply(strsplit(colnames(datBIP), '\\.'), function(x) x[[1]]), '_', c(rep('24', 6), rep('48', 6), rep('72', 6)))

      datBIP$DILI <- sapply(rownames(datBIP), function(x) DILI$binaryDILI[DILI$treatment == x]) 

      datBIP
    }

    assemble_CHOP <- function() {
      datCHOP  <- data.frame(do.call(rbind, mCHOP_24), do.call(rbind, mCHOP_48), do.call(rbind, mCHOP_72))
      datCHOP <- datCHOP[which(rowSums(is.na(datCHOP)) == 0), ]  
        colnames(datCHOP) <- paste0(sapply(strsplit(colnames(datCHOP), '\\.'), function(x) x[[1]]), '_', c(rep('24', 6), rep('48', 6), rep('72', 6)))

      datCHOP$DILI <- sapply(rownames(datCHOP), function(x) DILI$binaryDILI[DILI$treatment == x]) 

      datCHOP
    }
   
    assemble_P21 <- function() {
      datP21  <- data.frame(do.call(rbind, mP21_24), do.call(rbind, mP21_48), do.call(rbind, mP21_72))
      datP21 <- datP21[which(rowSums(is.na(datP21)) == 0), ]  
        colnames(datP21) <- paste0(sapply(strsplit(colnames(datP21), '\\.'), function(x) x[[1]]), '_', c(rep('24', 6), rep('48', 6), rep('72', 6)))

      datP21$DILI <- sapply(rownames(datP21), function(x) DILI$binaryDILI[DILI$treatment == x]) 
    }
    
   
  # All Combined
    assemeble_All <- function() {
      datAll <- data.frame(datBIP, datCHOP, datP21)
      datAll$DILI <- sapply(rownames(datAll), function(x) DILI$binaryDILI[DILI$treatment == x]) 

      datAll
    }
    

## Fit 'n' Predict!
  # Functions
  	.setTrainingIndex <- function(data) {
  		neg <- which(data$DILI == 0)
  		pos <- which(data$DILI == 1)

  		trainingNegative <- sample(neg, round(length(neg) * .7) )
  		trainingPositive <- sample(pos, round(length(pos) * .7) )

  		c(trainingNegative, trainingPositive)
  	}

  	fitPred <- function(data, type = c('glm', 'nn'), index = NULL, onlyIndex = FALSE, nnAccuracy = 0.05, nnNodes = 3) {	
      if(FALSE) {
        data <- assemble_24()
      }

  		maxs <- apply(data[, !grepl('DILI', colnames(data))], 2, max) 
  		mins <- apply(data[, !grepl('DILI', colnames(data))], 2, min)

  		scaled <- as.data.frame(scale(data[, !grepl('DILI', colnames(data))], center = mins, scale = maxs - mins))
  		data[, !grepl('DILI', colnames(data))] <- scaled

      if(is.null(index)) index <- sort(.setTrainingIndex(data))

      if(onlyIndex) return(index)

  		train <- data[index, ]
  		test  <- data[-index, ]

  		f <- as.formula(paste0('DILI ~ ', paste0(grep('CMAX', colnames(data), value = T), collapse = ' + ')))
    		
    		if(type == 'glm') {
    	    fit <- glm(f, data = train, family = binomial(logit), maxit = 1000)
  			  pr <- predict(fit, test)
  			  MSE <- sum((pr - test$DILI)^2)/nrow(test)
  		
  			  missClassifiedProportion <- 1 - (sum((round(pr) + test$DILI) == 2 | (round(pr) + test$DILI) == 0)/nrow(test))

    		}

    		if(type == 'nn') {
    			n <- names(train)
    			fit <- neuralnet(f, data = train, hidden = nnNodes, linear.output = !(classify), threshold = accuracy)

    			pr <- compute(fit, test[, grep('CMAX', colnames(test))])
          MSE <- sum((test$DILI - pr$net.result)^2)/nrow(test)

  			  missClassifiedProportion <- 1 - (sum((round(pr$net.result) + test$DILI) == 2 | (round(pr$net.result) + test$DILI) == 0)/nrow(test))
    		}

      return(list(index = index, fit = fit, pr = pr, mse = MSE, miss = missClassifiedProportion))  
  	}

  # Here we go
    set.seed(123)
    dat <- assemble_Rprtr_24()
    Index <- fitPred(dat$datBIP_24, 'glm', onlyIndex = TRUE)
   
    tmp <- list(fitPred(dat$datBIP_24, 'glm', index = Index), fitPred(dat$datBIP_24, 'nn', index = Index),
                fitPred(dat$datCHOP_24, 'glm', index = Index), fitPred(dat$datCHOP_24, 'nn', index = Index),
                fitPred(dat$datP21_24, 'glm', index = Index), fitPred(dat$datP21_24, 'nn', index = Index))

    Index <- fitPred(dat$datBIP_24, 'glm', onlyIndex = TRUE)
    tmp2 <- list(fitPred(dat$datBIP_24, 'glm', index = Index), fitPred(dat$datBIP_24, 'nn', index = Index),
                fitPred(dat$datCHOP_24, 'glm', index = Index), fitPred(dat$datCHOP_24, 'nn', index = Index),
                fitPred(dat$datP21_24, 'glm', index = Index), fitPred(dat$datP21_24, 'nn', index = Index))



























