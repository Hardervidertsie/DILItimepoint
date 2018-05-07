## Neural network for fun
## TODO
#	more reporters
#	integrate both glm as well as nn fits
#	check difference when round(pr)

## Settings
	options(stringsAsFactors = FALSE)
	
	library(data.table)
	library(preprocessCore)

	library(parallel)

	library(neuralnet)

	library(ggplot2)

	inputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'
	inputFiles <- list.files(inputDir, pattern = '.Rdata')	

## Data
  # DILI annotation
  	DILI <- read.delim('/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/DILI_annotation.txt')
  	  DILI$treatment <- toupper(DILI$treatment)
  	  DILI$binaryDILI <- DILI$doseDILI <- as.numeric(grepl('Most', DILI$updatedDILI) | grepl('Less', DILI$updatedDILI))
  	  DILI$doseDILI[grepl('Most', DILI$updatedDILI)] <- 2 

  	  DILI$treatment[DILI$treatment == "BUTHIONINE SULFOXAMINE"] <- "BUTHIONINE SULFOXIMINE"
  	  DILI$treatment[DILI$treatment == "CARBAMAZEPINE"] <- "CARBAMAZAPINE"
  	  DILI$treatment[DILI$treatment == "CYCLOHEXIMIDE"] <- "CYCLOHEXAMINE"
  	  DILI$treatment[DILI$treatment == "CYCLOSPORIN A"] <- "CYCLOSPORINA"
  	  DILI$treatment[DILI$treatment == "ETOPOSIDE"] <- "ETOPOSIDE "
  	  DILI$treatment[DILI$treatment == "SODIUM ARSENITE"] <- "SODIUM-ARSENITE"
  	  DILI$treatment[DILI$treatment == "STAUROSPORIN"] <- "STAUSPORIN"

  # BIP
	if(!exists('outputList')) load(paste0(inputDir, inputFiles[1]))	
	BIP <- lapply(outputList, function(x) x$sumData)

  # Parse: very clumsy MacGuyer work, but it's just a first!
  	calculateFCvsDMSO <- function(sumData) {
  		if(FALSE) sumData <- BIP[[1]]
  		
  		sumData$treatment <- toupper(sumData$treatment)
  		sumData <- data.frame(sumData)
  		sumData$cmax <- toupper(do.call(rbind, strsplit(as.character(sumData$plateID), '_'))[,2])
  		sumData$ID <- paste(sumData$treatment, paste0(sumData$timeID, 'hr'), sumData$cmax, sep = '_')

  		sumData <- sumData[order(sumData$plateID), ]
  		plateIDs <- as.character(unique(sumData$plateID))
  		
  		meanDMSO <- sapply(plateIDs, function(.pID) mean(sumData$Cytoplasm_Intensity_MeanIntensity_image_GFP[sumData$plateID == .pID & sumData$treatment == 'DMSO']))
  		
  		sumData$meanPlateDMSO <- NA
  		for(.pID in plateIDs) sumData$meanPlateDMSO[sumData$plateID == .pID] <- meanDMSO[.pID]

  		sumData$GFP_FC <- sumData$Cytoplasm_Intensity_MeanIntensity_image_GFP / sumData$meanPlateDMSO

  		rslt <- sumData[sumData$treatment != 'DMSO', ]
  		rslt <- rslt[order(rslt$ID), ]

  		rslt
  	}		

  	pBIP <- lapply(BIP, calculateFCvsDMSO)	
  	dat <- data.frame(pBIP[[1]][, c('treatment', 'timeID', 'cmax')],
  						rep1_meanIntensity = pBIP[[1]]$Cytoplasm_Intensity_MeanIntensity_image_GFP, 
  						rep1_meanPlateDMSO = pBIP[[1]]$meanPlateDMSO,
  						rep1_GFP_FC 	   = pBIP[[1]]$GFP_FC,
  						rep2_meanIntensity = pBIP[[2]]$Cytoplasm_Intensity_MeanIntensity_image_GFP, 
  						rep2_meanPlateDMSO = pBIP[[2]]$meanPlateDMSO,
  						rep2_GFP_FC 	   = pBIP[[2]]$GFP_FC,
  						rep3_meanIntensity = pBIP[[3]]$Cytoplasm_Intensity_MeanIntensity_image_GFP, 
  						rep3_meanPlateDMSO = pBIP[[3]]$meanPlateDMSO,
  						rep3_GFP_FC 	   = pBIP[[3]]$GFP_FC
  						)
  	
  	treatmentsToDrop <- unique(dat[which(!(dat$treatment %in% DILI$treatment)), 'treatment'])
  	dat <- dat[which(!dat$treatment %in% treatmentsToDrop), ]
  	dat$binaryDILI <- dat$doseDILI <- NA

  	for(trt in unique(dat$treatment)) dat[dat$treatment == trt, c('binaryDILI', 'doseDILI')] <- DILI[DILI$treatment == trt, c('binaryDILI', 'doseDILI')]

  	# short formats; replicates merged
  	dat24 <- data.frame(treatment = unique(dat$treatment), 
  						CMAX1 = NA, CMAX5 = NA, CMAX10 = NA, CMAX25 = NA, CMAX50 = NA, CMAX100 = NA)

  	dat24$CMAX1 <- sapply(dat24$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 24 & dat$cmax == '1CMAX'), grep('GFP_FC', colnames(dat))])))
	dat24$CMAX5 <- sapply(dat24$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 24 & dat$cmax == '5CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat24$CMAX10 <- sapply(dat24$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 24 & dat$cmax == '10CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat24$CMAX25 <- sapply(dat24$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 24 & dat$cmax == '25CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat24$CMAX50 <- sapply(dat24$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 24 & dat$cmax == '50CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat24$CMAX100 <- sapply(dat24$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 24 & dat$cmax == '100CMAX'), grep('GFP_FC', colnames(dat))])))
  	  	
  	dat48 <- data.frame(treatment = unique(dat$treatment), 
  						CMAX1 = NA, CMAX5 = NA, CMAX10 = NA, CMAX25 = NA, CMAX50 = NA, CMAX100 = NA)

  	dat48$CMAX1 <- sapply(dat48$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 48 & dat$cmax == '1CMAX'), grep('GFP_FC', colnames(dat))])))
	dat48$CMAX5 <- sapply(dat48$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 48 & dat$cmax == '5CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat48$CMAX10 <- sapply(dat48$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 48 & dat$cmax == '10CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat48$CMAX25 <- sapply(dat48$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 48 & dat$cmax == '25CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat48$CMAX50 <- sapply(dat48$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 48 & dat$cmax == '50CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat48$CMAX100 <- sapply(dat48$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 48 & dat$cmax == '100CMAX'), grep('GFP_FC', colnames(dat))])))
  	  	
  	dat72 <- data.frame(treatment = unique(dat$treatment), 
  						CMAX1 = NA, CMAX5 = NA, CMAX10 = NA, CMAX25 = NA, CMAX50 = NA, CMAX100 = NA)

  	dat72$CMAX1 <- sapply(dat72$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 72 & dat$cmax == '1CMAX'), grep('GFP_FC', colnames(dat))])))
	dat72$CMAX5 <- sapply(dat72$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 72 & dat$cmax == '5CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat72$CMAX10 <- sapply(dat72$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 72 & dat$cmax == '10CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat72$CMAX25 <- sapply(dat72$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 72 & dat$cmax == '25CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat72$CMAX50 <- sapply(dat72$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 72 & dat$cmax == '50CMAX'), grep('GFP_FC', colnames(dat))])))
  	dat72$CMAX100 <- sapply(dat72$treatment, function(trt) mean(unlist(dat[which(dat$treatment == trt & dat$timeID == 72 & dat$cmax == '100CMAX'), grep('GFP_FC', colnames(dat))])))
  	  	

  	dat24$DILI <- dat48$DILI <- dat72$DILI <- sapply(dat24$treatment, function(trt) DILI$binaryDILI[DILI$treatment == trt])

  	rownames(dat24) <- rownames(dat48) <- rownames(dat72) <-   dat24$treatment
  	dat24 <- dat24[, !grepl('treat', colnames(dat24))]
  	dat48 <- dat48[, !grepl('treat', colnames(dat48))]
  	dat72 <- dat72[, !grepl('treat', colnames(dat72))]

  	datAll <- data.frame(dat24, dat48[, grep('CMAX', colnames(dat48))], dat72[, grep('CMAX', colnames(dat72))])

## predictions

  # Functions
  	.setTrainingIndex <- function(data) {
  		trainingNegative <- sample(which(data$DILI == 0), round(sum(data$DILI == 0)*.7))
  		trainingPositive <- sample(which(data$DILI == 1), round(sum(data$DILI == 1)*.7))

  		c(trainingNegative, trainingNegative)
  	}

	fitGLM <- function(data, returnGLM = TRUE) {
		if(FALSE) {
			data <- dat24
			predictFactor <- TRUE
			set.seed(123)
		}

		maxs <- apply(data[, !grepl('DILI', colnames(data))], 2, max) 
		mins <- apply(data[, !grepl('DILI', colnames(data))], 2, min)

		scaled <- as.data.frame(scale(data[, !grepl('DILI', colnames(data))], center = mins, scale = maxs - mins))
		data[, !grepl('DILI', colnames(data))] <- scaled

		index <- .setTrainingIndex(data)
		train <- data[index,]
		test  <- data[-index,]

		f <- as.formula(paste0('DILI ~ ', paste0(grep('CMAX', colnames(data), value = T), collapse = ' + ')))
  		
		lm.fit <- glm(f, data=train, family = binomial(logit))
		summary(lm.fit)
		pr.lm <- predict(lm.fit,test)
		MSE.lm <- sum((pr.lm - test$DILI)^2)/nrow(test)
		
		missClassifiedProportion <- 1 - (sum((round(pr.lm) + test$DILI) == 2 | (round(pr.lm) + test$DILI) == 0)/nrow(test))

		if(returnGLM) return(list(glm = lm.fit, pr = pr.lm, 
								  mse = MSE.lm, 
								  missClassifiedProportion = missClassifiedProportion,
								  index = index))
		if(!returnGLM) return(MSE.lm)
	}

	fitNN <- function(data, returnNN = TRUE, nodes = 3, accuracy = .05, classify = TRUE) {
		maxs <- apply(data[, !grepl('DILI', colnames(data))], 2, max) 
		mins <- apply(data[, !grepl('DILI', colnames(data))], 2, min)

		scaled <- as.data.frame(scale(data[, !grepl('DILI', colnames(data))], center = mins, scale = maxs - mins))
		data[, !grepl('DILI', colnames(data))] <- scaled
		
		index <- .setTrainingIndex(data)
		train <- data[index,]
		test  <- data[-index,]

		n <- names(train)
		f <- as.formula(paste0('DILI ~ ', paste0(grep('CMAX', colnames(data), value = T), collapse = ' + ')))
  		nn <- neuralnet(f, data = train, hidden = nodes, linear.output = !(classify), threshold = accuracy)

  		pr.nn <- compute(nn, test[, grep('CMAX', colnames(test))])

		MSE.nn <- sum((test$DILI - pr.nn$net.result)^2)/nrow(test)

		missClassifiedProportion <- 1 - (sum((round(pr.nn$net.result) + test$DILI) == 2 | (round(pr.nn$net.result) + test$DILI) == 0)/nrow(test))

		
		if(returnNN) return(list(nn = nn, pr = pr.nn, mse = MSE.nn,  missClassifiedProportion = missClassifiedProportion))
		if(!returnNN) return(MSE.nn)
  	}

  	getPropMissCl <- function(fitOutput) {
  		sapply(fitOutput, function(x) x$missClassifiedProportion)
  	}

  set.seed(123)

  # glm
	glm_dat24.MSE <- lapply(1:100, function(x) fitGLM(dat24))
	glm_dat48.MSE <- lapply(1:100, function(x) fitGLM(dat48))
	glm_dat72.MSE <- lapply(1:100, function(x) fitGLM(dat72))
	glm_datAll.MSE <- lapply(1:100, function(x) fitGLM(datAll))

  # neural
  	nn_dat24.MSE <- lapply(1:100, function(x) fitNN(dat24))
	nn_dat48.MSE <- lapply(1:100, function(x) fitNN(dat48))
	nn_dat72.MSE <- lapply(1:100, function(x) fitNN(dat72))
	nn_datAll.MSE <- lapply(1:100, function(x) fitNN(datAll))

	toPlot <- cbind(getPropMissCl(glm_dat24.MSE), getPropMissCl(nn_dat24.MSE),
		  getPropMissCl(glm_dat48.MSE), getPropMissCl(nn_dat48.MSE),
		  getPropMissCl(glm_dat72.MSE), getPropMissCl(nn_dat72.MSE),
		  getPropMissCl(glm_datAll.MSE), getPropMissCl(nn_datAll.MSE))
	

	boxplot(toPlot, title = 'Misclassified propoertion',
  		ylim = c(0, 0.6))











