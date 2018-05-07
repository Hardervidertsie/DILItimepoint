## Settings
	options(stringsAsFactors = FALSE)
	
	library(neuralnet)
	library(nnetpredint)
	library(NeuralNetTools)

	library(data.table)
	library(preprocessCore)
	library(reshape2)
	library(parallel)

	inputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'


## Functions
	mergeReplicates <- function(parsedList, .debug = FALSE) {
		if(.debug) {
			parsedList <- SRXN1
			parsedList <- BIP
		}

		# Parseing
		sumDatas <- lapply(parsedList, function(rep) {
			tmp <- as.data.frame(rep$sumData)
			tmp$treatmentID <- paste0(tmp$treatment, '_', tmp$CMAX, '_', tmp$timeID)
			tmp$dose_uM <- as.numeric(as.character(tmp$dose_uM))

			tmp
		})


		# remove duplicated treatments; select highest dose for duplicates
		# duplicatedTreatments <- sort(names(which(rowSums(sapply(sumDatas, function(tmp) table(tmp$treatmentID) ) > 1) > 0)))
		duplicatedTreatments <- sort(names(which(table(unlist(sapply(sumDatas, function(tmp) names(which(table(tmp$treatmentID) > 1))))) == 3)))
			indices2drop <- lapply(duplicatedTreatments, function(trtID) {
				inds <- sapply(sumDatas, function(tmp) {
					maxDose <- max(tmp[tmp$treatmentID == trtID, 'dose_uM'])
					doses2drop <- which(tmp$treatmentID == trtID & tmp$dose_uM != maxDose)
					
					doses2drop
				})

				if(class(inds) == 'matrix') {
					colnames(inds) <- c('rep1', 'rep2', 'rep3')
				} else {
					names(inds) <- c('rep1', 'rep2', 'rep3')
				}

				inds
			})

		rep1_indices2drop <- unlist(sapply(indices2drop, function(x) {
			if(class(x) == 'matrix') x[, 'rep1'] else x['rep1']
		}))

		rep2_indices2drop <- unlist(sapply(indices2drop, function(x) {
			if(class(x) == 'matrix') x[, 'rep2'] else x['rep2']
		}))

		rep3_indices2drop <- unlist(sapply(indices2drop, function(x) {
			if(class(x) == 'matrix') x[, 'rep3'] else x['rep3']
		}))


		sumDatas[[1]] <- sumDatas[[1]][-rep1_indices2drop, ]
		sumDatas[[2]] <- sumDatas[[2]][-rep2_indices2drop, ]
		sumDatas[[3]] <- sumDatas[[3]][-rep3_indices2drop, ]
		

		experimentsInAllReplicates <- sort(names(which(table(unlist(sapply(sumDatas, function(x) x$treatmentID))) == 3)))
		sumDatas <- lapply(sumDatas, function(rep) {
			tmp <- rep[which(rep$treatmentID %in% experimentsInAllReplicates), ]
			rownames(tmp) <- tmp$treatmentID

			tmp
		})

		sumDatas <- lapply(sumDatas, function(rep) rep[sort(rownames(sumDatas[[1]])), ])


		# Merging
		features2merge <- c(grep('Intensity_MeanIntensity_', colnames(sumDatas[[1]]), value = TRUE),
							'numberOfObjects', 'plateMeanDMSO', 'plateMedianDMSO', 
							grep('_Above_', colnames(sumDatas[[1]]), value = TRUE),
							grep('AnV_', colnames(sumDatas[[1]]), value = TRUE),
							grep('PI_' , colnames(sumDatas[[1]]), value = TRUE))

		featureMeans <- data.frame(lapply(features2merge, function(ftr) {
			cat(ftr, '\n')
			rowMeans(sapply(sumDatas, function(rep) {
				tmp <- rep[, ftr]

				if(class(tmp) == 'factor')    stop(paste(ftr, 'issue.\n'))
				if(class(tmp) == 'character') tmp <- as.numeric(tmp)

				tmp
			}), na.rm = TRUE)
		}))
		
		rownames(featureMeans) <- experimentsInAllReplicates
		colnames(featureMeans) <- features2merge

		# return
		rslt <- data.frame(sumDatas[[1]][experimentsInAllReplicates, c('treatment', 'timeID', 'dose_uM', 'CMAX', 'CMAX_numeric', 'matchDose', 'matchDoseTime')],
						   featureMeans[experimentsInAllReplicates, ])

		return(rslt)
	}

	reshapeMergedReplicates <- function(meanReporter, reporterName, features = c('CMAX_numeric', 'Fraction_Above_1_medianIntregated_plateDMSO', 'AnV_fractionNonNA', 'PI_fractionNonNA'), timePoint, .debug = FALSE) {
		if(.debug) {
			meanReporter <- meanBIP
			features <- c('CMAX_numeric', 'Fraction_Above_2_meanIntregated_plateDMSO', 'AnV_fractionNonNA', 'PI_fractionNonNA')
			timePoint <- 24
		}

		rslt <- DILI[, c('treatment', 'binaryDILI')]
		rslt[, paste0(c(1, 5, 10, 25, 50, 100), 'CMAX_FractionAboveDMSO')] <- NA
		rslt[, paste0(c(1, 5, 10, 25, 50, 100), 'CMAX_nonNA_fractionAnV')] <- NA
		rslt[, paste0(c(1, 5, 10, 25, 50, 100), 'CMAX_nonNA_fractionPI')]  <- NA
			rownames(rslt) <- rslt$treatment
		
		for(trtID in rownames(rslt)) {
			if(length(which(meanReporter$treatment == trtID)) == 0) next
			
			for(ftr in grep('CMAX', colnames(rslt), value = TRUE)) {
				cmax     <- as.numeric(unlist(strsplit(ftr, 'CMAX'))[1])
				if(grepl('DMSO', ftr)) fraction <- 'DMSO'
				if(grepl('PI', ftr))   fraction <- 'PI'
				if(grepl('AnV', ftr))  fraction <- 'AnV'

				.ind <- which(meanReporter$treatment == trtID & meanReporter$timeID == timePoint & meanReporter$CMAX_numeric == cmax)
				.ftr <- grep(fraction, features, value = TRUE)  

				if(length(.ind) == 0) next

				rslt[trtID, ftr] <- meanReporter[.ind, .ftr]

			}
		}

		rslt <- rslt[,colSums(is.na(rslt)) != nrow(rslt)]
		rslt <- rslt[which(rowSums(is.na(rslt)) == 0), ]

		rslt$reporterLine <- reporterName

		rslt
	}

	setTrainingIndex <- function(data, class = c('binaryDILI', 'doseDILI')) {
		if(length(class) != 1) stop('Only supply a single column name as dependent variable.\n')
  		
  		neg <- which(data[, class] == 0)
  		pos <- which(data[, class] == 1)

  		trainingNegative <- sample(neg, round(length(neg) * .7) )
  		trainingPositive <- sample(pos, round(length(pos) * .7) )

  		c(trainingNegative, trainingPositive)
  	}


## Screening data
	if(FALSE) {
		load(paste0(inputDir, 'parsedBIP.RData'))
		meanBIP <- mergeReplicates(BIP)
		meanBIP$treatment <- gsub('SODIUM-ARSENITE', 'SODIUMARSENITE', gsub('STAUSPORIN', 'STAUROSPORIN', gsub('BUTHIONINESULFOXIMINE', 'BUTHIONINESULFOXAMINE', meanBIP$treatment)))
			save(meanBIP, file = paste0(inputDir, 'mergedReplicates_sumDat_BIP.RData'))
			# rm(BIP)

		load(paste0(inputDir, 'parsedSRXN1.RData'))
		meanSRXN1 <- mergeReplicates(SRXN1)	
			save(meanSRXN1, file = paste0(inputDir, 'mergedReplicates_sumDat_SRXN1.RData'))
			# rm(SRXN1)

		load(paste0(inputDir, 'parsedP21.RData'))
		meanP21 <- mergeReplicates(P21)	
			save(meanP21, file = paste0(inputDir, 'mergedReplicates_sumDat_P21.RData'))
			# rm(P21)

		load(paste0(inputDir, 'parsedCHOP.RData'))
		meanCHOP <- mergeReplicates(CHOP)	
			save(meanCHOP, file = paste0(inputDir, 'mergedReplicates_sumDat_CHOP.RData'))
			# rm(CHOP)

		load(paste0(inputDir, 'parsedICAM1.RData'))
		meanICAM1 <- mergeReplicates(ICAM1)	
			save(meanICAM1, file = paste0(inputDir, 'mergedReplicates_sumDat_ICAM1.RData'))
			# rm(ICAM1)

		load(paste0(inputDir, 'parsedHSPA1B.RData'))
		meanHSPA1B <- mergeReplicates(HSPA1B)
		meanHSPA1B$treatment <- gsub('SODIUM-ARSENITE', 'SODIUMARSENITE', gsub('STAUSPORIN', 'STAUROSPORIN', gsub('BUTHIONINESULFOXIMINE', 'BUTHIONINESULFOXAMINE', meanHSPA1B$treatment)))
			save(meanHSPA1B, file = paste0(inputDir, 'mergedReplicates_sumDat_HSPA1B.RData'))
			# rm(HSPA1B)	
	}

	load(paste0(inputDir, 'mergedReplicates_sumDat_BIP.RData'))
	load(paste0(inputDir, 'mergedReplicates_sumDat_SRXN1.RData'))
	load(paste0(inputDir, 'mergedReplicates_sumDat_P21.RData'))
	load(paste0(inputDir, 'mergedReplicates_sumDat_CHOP.RData'))
	load(paste0(inputDir, 'mergedReplicates_sumDat_ICAM1.RData'))
	load(paste0(inputDir, 'mergedReplicates_sumDat_HSPA1B.RData'))
	

## DILI annotation data
	DILI <- read.delim('/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/DILI_annotation.txt')

	DILI$binaryDILI <- DILI$doseDILI <- as.numeric(grepl('Most', DILI$updatedDILI) | grepl('Less', DILI$updatedDILI))
  	DILI$doseDILI[grepl('Most', DILI$updatedDILI)] <- 2 

  	DILI$treatment <- toupper(DILI$treatment)
  	  DILI$treatment[DILI$treatment == "BUTHIONINESULFOXIMINE"] <- "BUTHIONINESULFOXAMINE"
  	  DILI$treatment[DILI$treatment == "CARBAMAZEPINE"]         <- "CARBAMAZAPINE"
  	  DILI$treatment[DILI$treatment == "CYCLOHEXIMIDE"]         <- "CYCLOHEXAMINE"
  	  DILI$treatment[DILI$treatment == "CYCLOSPORIN A"]         <- "CYCLOSPORINA"
  	  DILI$treatment[DILI$treatment == "SODIUM ARSENITE"]       <- "SODIUMARSENITE"
  	
  	  DILI$treatment <- gsub(' ', '', DILI$treatment)

  	  rownames(DILI) <- DILI$treatment


## Reshape reporter data according to DILI annotation
	annotatedBIP <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates(meanBIP,    ,reporterName = 'BIP',    timePoint = .timePoint))
	annotatedSRX <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates(meanSRXN1,  ,reporterName = 'SRXN1',  timePoint = .timePoint))
	annotatedP21 <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates(meanP21,    ,reporterName = 'P21',    timePoint = .timePoint))
	annotatedCHP <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates(meanCHOP,   ,reporterName = 'CHOP',   timePoint = .timePoint))
	annotatedHSP <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates(meanHSPA1B, ,reporterName = 'HSPA1B', timePoint = .timePoint))
	annotatedICM <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates(meanICAM1,  ,reporterName = 'ICAM1',  timePoint = .timePoint))	# misses 24 hour data
		names(annotatedBIP) <- names(annotatedSRX) <- names(annotatedP21) <- names(annotatedCHP) <- names(annotatedHSP) <- names(annotatedICM) <- c('h24', 'h48', 'h72') 

	annotatedAll <- list(BIP    = annotatedBIP, 
						 SRXN1  = annotatedSRX, 
						 P21    = annotatedP21, 
						 CHOP   = annotatedCHP, 
						 HSPA1B = annotatedHSP, 
						 ICAM1  = annotatedICM)

	
## Create objects to input into classifiers
  ## Combine reporters per timepoint
  # 24 hour:	# 154 treatments in all but ICAM1: 108 DILI / 46 nonDILI
  	treatments_24h <- names(which(table(unlist(lapply(annotatedAll[-which(names(annotatedAll) == 'ICAM1')], function(x) { rownames(x$h24) } ))) == 5))
  	classifierInput_24 <- data.frame(binaryDILI = DILI[treatments_24h, 'binaryDILI'],
  									 doseDILI = DILI[treatments_24h, 'doseDILI'],
									 lapply(annotatedAll[-which(names(annotatedAll) == 'ICAM1')], function(x) {
									   tmp <- x$h24[treatments_24h, 3:ncol(x$h24)]
									   tmp[, -grep('reporterLine', colnames(tmp))]
							  		 }))

  # 48 hour:	# 171 treatments: 112 DILI / 59 nonDILI
  	treatments_48h <- names(which(table(unlist(lapply(annotatedAll, function(x) { rownames(x$h48) } ))) == 6))
  	classifierInput_48 <- data.frame(binaryDILI = DILI[treatments_48h, 'binaryDILI'],
  									 doseDILI = DILI[treatments_48h, 'doseDILI'],
									 lapply(annotatedAll, function(x) {
									   tmp <- x$h48[treatments_48h, 3:ncol(x$h48)]
									   tmp[, -grep('reporterLine', colnames(tmp))]
							  		 }))

  # 72 hour:	# 177 treatments in all but ICAM1: 117 DILI / 60 nonDILI
  	treatments_72h <- names(which(table(unlist(lapply(annotatedAll, function(x) { rownames(x$h72) } ))) == 6))
  	classifierInput_72 <- data.frame(binaryDILI = DILI[treatments_72h, 'binaryDILI'],
  									 doseDILI = DILI[treatments_72h, 'doseDILI'],
									 lapply(annotatedAll, function(x) {
									   tmp <- x$h72[treatments_72h, 3:ncol(x$h72)]
									   tmp[, -grep('reporterLine', colnames(tmp))]
							  		 }))


  ## Combine timepoints per reporter
  # BIP:		# 177 treatments: 117 DILI / 60 nonDILI
    treatments_BIP <- names(which(table(unlist(lapply(annotatedBIP, rownames))) == 3))
    classifierInput_BIP <- data.frame(binaryDILI = DILI[treatments_BIP, 'binaryDILI'],
    								  doseDILI = DILI[treatments_BIP, 'doseDILI'],
    								  lapply(annotatedBIP, function(x) {
    								  	x[treatments_BIP, 3:(ncol(x)-1)]
    								  }))

  # SRXN1 		# 149 treatments: 104 DILI / 45 nonDILI
    treatments_SRXN1 <- names(which(table(unlist(lapply(annotatedSRX, rownames))) == 3))
    classifierInput_SRXN1 <- data.frame(binaryDILI = DILI[treatments_SRXN1, 'binaryDILI'],
    									doseDILI = DILI[treatments_SRXN1, 'doseDILI'],
    									lapply(annotatedSRX, function(x) {
    									  x[treatments_SRXN1, 3:(ncol(x)-1)]
    									}))

  # P21 		# 178 treatments: 118 DILI / 60 nonDILI
    treatments_P21 <- names(which(table(unlist(lapply(annotatedP21, rownames))) == 3))
    classifierInput_P21 <- data.frame(binaryDILI = DILI[treatments_P21, 'binaryDILI'],
    									doseDILI = DILI[treatments_P21, 'doseDILI'],
    									lapply(annotatedP21, function(x) {
    									  x[treatments_P21, 3:(ncol(x)-1)]
    									}))

  # CHOP 		# 178 treatments: 118 DILI / 60 nonDILI
    treatments_CHOP <- names(which(table(unlist(lapply(annotatedCHP, rownames))) == 3))
    classifierInput_CHOP <- data.frame(binaryDILI = DILI[treatments_CHOP, 'binaryDILI'],
    									doseDILI = DILI[treatments_CHOP, 'doseDILI'],
    									lapply(annotatedCHP, function(x) {
    									  x[treatments_CHOP, 3:(ncol(x)-1)]
    									}))

  # HSPA1B 		# 177 treatments: 117 DILI / 60 nonDILI
    treatments_HSPA1B <- names(which(table(unlist(lapply(annotatedHSP, rownames))) == 3))
    classifierInput_HSPA1B <- data.frame(binaryDILI = DILI[treatments_HSPA1B, 'binaryDILI'],
    									doseDILI = DILI[treatments_HSPA1B, 'doseDILI'],
    									lapply(annotatedHSP, function(x) {
    									  x[treatments_HSPA1B, 3:(ncol(x)-1)]
    									}))

  # ICAM1 		# 178 treatments: 118 DILI / 60 nonDILI
    treatments_ICAM1 <- names(which(table(unlist(lapply(annotatedICM, rownames))) == 3))
    classifierInput_ICAM1 <- data.frame(binaryDILI = DILI[treatments_ICAM1, 'binaryDILI'],
    									doseDILI = DILI[treatments_ICAM1, 'doseDILI'],
    									lapply(annotatedICM, function(x) {
    									  x[treatments_ICAM1, 3:(ncol(x)-1)]
    									}))


## Neuralnetting
  # 24 hour dataset
  	set.seed(1337)

  	INDEX <- list()
	for(i in 1:10000) INDEX[[i]] <- setTrainingIndex(classifierInput_24, 'binaryDILI')

	fitNN <- function(dat, NodesLayers = c(10, 3), splitBy = c('binaryDILI', 'doseDILI'), index = NULL, .debug = FALSE) {
	  # TO DO
	  	# 1) CMAX level as weight + linear.output boolean
	  	# 2) independent variable selection
	  	# 3) nnetpredint implementation for CIs
	  	# 4) linear.output/err.fct dependent on binaryDILI/doseDILI
	  	# 5) Compute confidence intervals using SEM

	  # Debugging	
		if(.debug) {
			NodesLayers <- c(2, 8, 4)
			splitBy 	<- 'doseDILI'
			index       <- INDEX[[1]]

			dat         <- classifierInput_BIP[, c(1,2, grep('DMSO', colnames(classifierInput_BIP)))]
			dat$doseDILI <- as.numeric(dat$doseDILI == 2)
			colnames(dat) <- gsub('_FractionAboveDMSO', '', colnames(dat))
		}

	  # Parseing
	  	# dat[, class] <- factor(dat[, class], order = TRUE)	

	  # Scale
		maxs <- apply(dat[, grepl('CMAX', colnames(dat))], 2, max) 
  		mins <- apply(dat[, grepl('CMAX', colnames(dat))], 2, min)

  		scaled <- as.data.frame(scale(dat[, grepl('CMAX', colnames(dat))], center = mins, scale = maxs - mins))
  		  if(identical(rownames(dat), rownames(scaled))) scaled <- data.frame(dat[, 1:2], scaled) else stop('Something went wrong.\n')

  	  # Data split
  	  	if(is.null(index)) index <- sort(setTrainingIndex(dat), splitBy)
	
  	  	train <- scaled[index,  ]
  	  	test  <- scaled[-index, ]

  	  # Training neural network
  	  	f <- as.formula(paste0('binaryDILI + doseDILI ~ ', paste0(grep('CMAX', colnames(scaled), value = T), collapse = ' + ')))
    	fit <- neuralnet(formula = f, 
    					 data = train, 
    					 hidden = NodesLayers, 
    					 linear.output = FALSE, 
    					 threshold = 0.05, 
    					 err.fct = 'ce',
    					 rep = 5,
    					 likelihood = TRUE)

      # Evaluation
      	pr <- list(binaryDILI = sapply(1:length(fit$weights), function(i) compute(fit, test[, grep('CMAX', colnames(test))], rep = i)$net.result[, 1]), 
      			   doseDILI   = sapply(1:length(fit$weights), function(i) compute(fit, test[, grep('CMAX', colnames(test))], rep = i)$net.result[, 2]))

      	prSummary <- lapply(pr, function(x) {
      		.mean <- rowMeans(x)
      		.sem  <- apply(x, 1, sd) / sqrt(ncol(x))

      		.lower <- .mean - (1.96 * .sem)
      		.upper <- .mean + (1.96 * .sem)

      		data.frame(mean = .mean, lower = .lower, upper = .upper)
      	})

      	plotnet(fit, alpha = 0.5, pos_col = 'steelblue3', neg_col = 'orangered3', cex = 0.8)











      	predictedDILI <- lapply(1:length(fit$weights), function(i) {
      		net.result <- 
      	


      		colnames(net.result) <- c(c('pred_binaryDILI', 'pred_doseDILI'))

      		net.result
      	})

      	predBinaryDILI <- 

      	apply(sapply(predictedDILI, function(x) x[, 1]), 1, sd)/sqrt(length(predictedDILI))*1.96



      	pr <- compute(fit, test[, grep('CMAX', colnames(test))])






      # Testing neural network
      	ci <- nnetPredInt(object = fit,
      				xTrain = scaled[index, grep('CMAX', colnames(scaled))],
      				yTrain = scaled[index, 1:2],
      				newData = scaled[-index, grep('CMAX', colnames(scaled))])

      	a <- cbind(test$binaryDILI, ci)




      	#
    	pr <- compute(fit, test[, grep('CMAX', colnames(test))])$net.result
    	mse <- sum((test[, class] - pr)^2) / nrow(test)


	}








		
	wrapperFitterAUC <- function(timeID = 24, iter = 500) {
		.bip <- reshapeMergedReplicates(meanBIP, timePoint =s timeID)
		.srx <- reshapeMergedReplicates(meanSRXN1, timePoint = timeID)
		.p21 <- reshapeMergedReplicates(meanP21, timePoint = timeID)
		.chp <- reshapeMergedReplicates(meanCHOP, timePoint = timeID)
		.hsp <- reshapeMergedReplicates(meanHSPA1B, timePoint = timeID)

		compounds2check <- names(which(table(c(rownames(.bip), rownames(.srx), rownames(.p21), rownames(.chp), rownames(.icm), rownames(.hsp))) == 6))
		
		.bip <- .bip[compounds2check, ]
		.srx <- .srx[compounds2check, ]
		.p21 <- .p21[compounds2check, ]
		.chp <- .chp[compounds2check, ]
		.hsp <- .hsp[compounds2check, ]

		predictionsReporter <- data.frame(bip = NA, srx = NA, p21 = NA, chp = NA, hsp = NA)

		for(i in 1:iter) {
			predictionsReporter[i, 'bip'] <- fitPredDMSO_AUC(.bip, index = INDEX[[i]])
			predictionsReporter[i, 'srx'] <- fitPredDMSO_AUC(.srx, index = INDEX[[i]])
			predictionsReporter[i, 'p21'] <- fitPredDMSO_AUC(.p21, index = INDEX[[i]])
			predictionsReporter[i, 'chp'] <- fitPredDMSO_AUC(.chp, index = INDEX[[i]])
			predictionsReporter[i, 'hsp'] <- fitPredDMSO_AUC(.hsp, index = INDEX[[i]])
		}

		predictionsReporter
	}

	fitPredDMSO_AUC <- function(dat, index = NULL) {
		if(FALSE) dat <- .bip
		if(FALSE) index <- INDEX[[1]]

		# maxs <- apply(dat[, grepl('CMAX', colnames(dat))], 2, max) 
  		# mins <- apply(dat[, grepl('CMAX', colnames(dat))], 2, min)

  		# scaled <- as.data.frame(scale(dat[, grepl('CMAX', colnames(dat))], center = mins, scale = maxs - mins))
  		# dat[, grepl('CMAX', colnames(dat))] <- scaled

  		colnames(dat) <- paste0('a_', colnames(dat))
  		dat <- dat[, c('a_binaryDILI', grep('DMSO', colnames(dat), value = TRUE))]

  		if(is.null(index)) index <- sort(.setTrainingIndex(dat))

  		train <- dat[index, ]
  		test  <- dat[-index, ]

  		f <- as.formula(paste0('a_binaryDILI ~ ', paste0(grep('CMAX', colnames(dat), value = T), collapse = ' + ')))

  		fit <- glm(f, data = train, family = binomial(logit), maxit = 1000)
  		pr  <- predict(fit, test)
  		
  		# from r-bloggers
  		# http://stats.stackexchange.com/questions/145566/how-to-calculate-area-under-the-curve-auc-or-the-c-statistic-by-hand
  		truestat <- test$a_binaryDILI
  		testres  <- round(pr)

  			tab <- as.matrix(table(truestat, testres))
  			tot <- colSums(tab)
  			truepos  <- unname(rev(cumsum(rev(tab[2,]))))
  			falsepos <- unname(rev(cumsum(rev(tab[1,]))))
  			totpos <- sum(tab[2, ])
  			totneg <- sum(tab[1, ])

  			sens <- truepos/totpos
  			omspec <- falsepos/totneg

  			sens   <- c(sens, 0)
  			omspec <- c(omspec, 0)

  			height <- (sens[-1]+sens[-length(sens)])/2
			width  <- -diff(omspec) 
		
		sum(height*width)	# AUC

	}

	INDEX <- list()
	for(i in 1:10000) INDEX[[i]] <- .setTrainingIndex(.hsp)

	preds24_AUC <- wrapperFitterAUC(timeID = 24, iter = 500)
	preds48_AUC <- wrapperFitterAUC(timeID = 48, iter = 500)
	preds72_AUC <- wrapperFitterAUC(timeID = 72, iter = 500)

	dev.new()
	par(mfrow = c(1, 3))
	boxplot(preds24_AUC, ylim = c(0.2, 0.8))
		abline(a = .5, b = 0)
	boxplot(preds48_AUC, ylim = c(0.2, 0.8))
		abline(a = .5, b = 0)
	boxplot(preds72_AUC, ylim = c(0.2, 0.8))
		abline(a = .5, b = 0)

	library(neuralnet)
	fitNNDMSO_AUC <- function(dat, nnNodes, errorMeasure = c('mse', 'auc'), index = NULL) {
		if(FALSE) dat <- .bip
		if(FALSE) index <- INDEX[[1]]

		 maxs <- apply(dat[, grepl('CMAX', colnames(dat))], 2, max) 
  		 mins <- apply(dat[, grepl('CMAX', colnames(dat))], 2, min)

  		 scaled <- as.data.frame(scale(dat[, grepl('CMAX', colnames(dat))], center = mins, scale = maxs - mins))
  		 dat[, grepl('CMAX', colnames(dat))] <- scaled

  		colnames(dat) <- paste0('a_', colnames(dat))
  		 # dat <- dat[, c('a_binaryDILI', grep('DMSO', colnames(dat), value = TRUE))]

  		if(is.null(index)) index <- sort(.setTrainingIndex(dat))

  		train <- dat[index, ]
  		test  <- dat[-index, ]

  		f <- as.formula(paste0('a_binaryDILI ~ ', paste0(grep('CMAX', colnames(dat), value = T), collapse = ' + ')))

  			n <- names(train)
    		fit <- neuralnet(f, data = train, hidden = nnNodes, linear.output = FALSE,  threshold = 0.05, rep = 10)

    		pr <- compute(fit, test[, grep('CMAX', colnames(test))])$net.result
    		mse <- sum((test$a_binaryDILI - pr)^2) / nrow(test)
    		# table(rowSums(cbind(as.numeric(pr - median(pr) > 0), test$a_binaryDILI)))
  		
  		# from r-bloggers
  		# http://stats.stackexchange.com/questions/145566/how-to-calculate-area-under-the-curve-auc-or-the-c-statistic-by-hand
  		truestat <- test$a_binaryDILI
  		testres  <- round(pr)

  			tab <- as.matrix(table(truestat, testres))
  			tot <- colSums(tab)
  			truepos  <- unname(rev(cumsum(rev(tab[2,]))))
  			falsepos <- unname(rev(cumsum(rev(tab[1,]))))
  			totpos <- sum(tab[2, ])
  			totneg <- sum(tab[1, ])

  			sens <- truepos/totpos
  			omspec <- falsepos/totneg

  			sens   <- c(sens, 0)
  			omspec <- c(omspec, 0)

  			height <- (sens[-1]+sens[-length(sens)])/2
			width  <- -diff(omspec) 
		
		sum(height*width)	# AUC
	}

	wrapperNNAUC <- function(timeID = 24, iter = 500, ...) {
		.bip <- reshapeMergedReplicates(meanBIP, timePoint = timeID)
		.srx <- reshapeMergedReplicates(meanSRXN1, timePoint = timeID)
		.p21 <- reshapeMergedReplicates(meanP21, timePoint = timeID)
		.chp <- reshapeMergedReplicates(meanCHOP, timePoint = timeID)
		.hsp <- reshapeMergedReplicates(meanHSPA1B, timePoint = timeID)

		compounds2check <- names(which(table(c(rownames(.bip), rownames(.srx), rownames(.p21), rownames(.chp), rownames(.icm), rownames(.hsp))) == 6))
		
		.bip <- .bip[compounds2check, ]
		.srx <- .srx[compounds2check, ]
		.p21 <- .p21[compounds2check, ]
		.chp <- .chp[compounds2check, ]
		.hsp <- .hsp[compounds2check, ]

		

		predictionsReporter <- as.data.frame(do.call(rbind, mclapply(1:iter, function(i) {
			cat(i, '\n')
			#predictionsReporter[i, 'bip'] <- fitNNDMSO_AUC(.bip, index = INDEX[[i]])
			#predictionsReporter[i, 'srx'] <- fitNNDMSO_AUC(.srx, index = INDEX[[i]])
			#predictionsReporter[i, 'p21'] <- fitNNDMSO_AUC(.p21, index = INDEX[[i]])
			#predictionsReporter[i, 'chp'] <- fitNNDMSO_AUC(.chp, index = INDEX[[i]])
			#predictionsReporter[i, 'hsp'] <- fitNNDMSO_AUC(.hsp, index = INDEX[[i]])

			c(fitNNDMSO_AUC(.bip, index = INDEX[[i]], ...), 
			  fitNNDMSO_AUC(.srx, index = INDEX[[i]], ...),
			  fitNNDMSO_AUC(.p21, index = INDEX[[i]], ...),
			  fitNNDMSO_AUC(.chp, index = INDEX[[i]], ...),
			  fitNNDMSO_AUC(.hsp, index = INDEX[[i]], ...))

		}, mc.cores = 5)))

		colnames(predictionsReporter) <- c('bip', 'srx', 'p21', 'chp', 'hsp')

		predictionsReporter
	}

	# 1 layer, 12 neurons
	nn24_12 <- wrapperNNAUC(timeID = 24, iter = 50, nnNodes = 12)
	nn48_12 <- wrapperNNAUC(timeID = 48, iter = 50, nnNodes = 12)
	nn72_12 <- wrapperNNAUC(timeID = 72, iter = 50, nnNodes = 12)

	# 1 layer, 6 neurons
	nn24_6 <- wrapperNNAUC(timeID = 24, iter = 50, nnNodes = 6)
	nn48_6 <- wrapperNNAUC(timeID = 48, iter = 50, nnNodes = 6)
	nn72_6 <- wrapperNNAUC(timeID = 72, iter = 50, nnNodes = 6)

	# 2 layers, 12 and 4 neurons
	nn24_12.4 <- wrapperNNAUC(timeID = 24, iter = 50, nnNodes = c(12, 4))
	nn48_12.4 <- wrapperNNAUC(timeID = 48, iter = 50, nnNodes = c(12, 4))
	nn72_12.4 <- wrapperNNAUC(timeID = 72, iter = 50, nnNodes = c(12, 4))

	# 2 layers, 12 and 8 neurons
	nn24_12.8 <- wrapperNNAUC(timeID = 24, iter = 50, nnNodes = c(12, 8))
	nn48_12.8 <- wrapperNNAUC(timeID = 48, iter = 50, nnNodes = c(12, 8))
	nn72_12.8 <- wrapperNNAUC(timeID = 72, iter = 50, nnNodes = c(12, 8))

	# 3 layers, 12, 8 and 4 neurons
	nn24_12.8.4 <- wrapperNNAUC(timeID = 24, iter = 50, nnNodes = c(12, 8, 4))
	nn48_12.8.4 <- wrapperNNAUC(timeID = 48, iter = 50, nnNodes = c(12, 8, 4))
	nn72_12.8.4 <- wrapperNNAUC(timeID = 72, iter = 50, nnNodes = c(12, 8, 4))


	## With all reporters in there
	nnFitAll <- function(timeID = 24, nnNodes = c(60, 30, 12, 4), iter = 500, ...) {
		.bip <- reshapeMergedReplicates(meanBIP, timePoint = timeID)
		.srx <- reshapeMergedReplicates(meanSRXN1, timePoint = timeID)
		.p21 <- reshapeMergedReplicates(meanP21, timePoint = timeID)
		.chp <- reshapeMergedReplicates(meanCHOP, timePoint = timeID)
		.hsp <- reshapeMergedReplicates(meanHSPA1B, timePoint = timeID)

		compounds2check <- names(which(table(c(rownames(.bip), rownames(.srx), rownames(.p21), rownames(.chp), rownames(.icm), rownames(.hsp))) == 6))
		
		.bip <- .bip[compounds2check, ]
		.srx <- .srx[compounds2check, ]
		.p21 <- .p21[compounds2check, ]
		.chp <- .chp[compounds2check, ]
		.hsp <- .hsp[compounds2check, ]

		.all <- data.frame(treatment = .bip$treatment, binaryDILI = .bip$binaryDILI,
							.bip[, grep('CMAX', colnames(.bip))], 
							.srx[, grep('CMAX', colnames(.srx))],
							.p21[, grep('CMAX', colnames(.p21))],
							.chp[, grep('CMAX', colnames(.chp))],
							.hsp[, grep('CMAX', colnames(.hsp))])

		 maxs <- apply(.all[, grepl('CMAX', colnames(.all))], 2, max) 
  		 mins <- apply(.all[, grepl('CMAX', colnames(.all))], 2, min)

  		 scaled <- as.data.frame(scale(.all[, grepl('CMAX', colnames(.all))], center = mins, scale = maxs - mins))
  		 .all[, grepl('CMAX', colnames(.all))] <- scaled

  		 
  		 
  		 AUCs <- mclapply(1:iter, function(i) {
  		 	index <- INDEX[[i]]

  		 	train <- .all[index, ]
  		 	test  <- .all[-index, ]

  		 	f <- as.formula(paste0('binaryDILI ~ ', paste0(grep('CMAX', colnames(.all), value = T), collapse = ' + ')))

  			n <- names(train)
    		fit <- neuralnet(f, data = train, hidden = nnNodes, linear.output = FALSE,  threshold = 0.05, rep = 50)

    		pr <- compute(fit, test[, grep('CMAX', colnames(test))])$net.result
    		mse <- sum((test$binaryDILI - pr)^2) / nrow(test)
    		
    		# table(rowSums(cbind(as.numeric(pr - median(pr) > 0), test$a_binaryDILI)))
	  		# from r-bloggers
	  		# http://stats.stackexchange.com/questions/145566/how-to-calculate-area-under-the-curve-auc-or-the-c-statistic-by-hand
	  		truestat <- test$binaryDILI
	  		testres  <- pr

	  			tab <- as.matrix(table(truestat, testres))
	  			tot <- colSums(tab)
	  			truepos  <- unname(rev(cumsum(rev(tab[2,]))))
	  			falsepos <- unname(rev(cumsum(rev(tab[1,]))))
	  			totpos <- sum(tab[2, ])
	  			totneg <- sum(tab[1, ])

	  			sens <- truepos/totpos
	  			omspec <- falsepos/totneg

	  			sens   <- c(sens, 0)
	  			omspec <- c(omspec, 0)

	  			height <- (sens[-1]+sens[-length(sens)])/2
				width  <- -diff(omspec) 
			
			sum(height*width)	# AUC
  		 }, mc.cores = 5)
  		
  	}









## Prelim Run
	.bip <- reshapeMergedReplicates(meanBIP, timePoint = 24)
	.srx <- reshapeMergedReplicates(meanSRXN1, timePoint = 24)
	.p21 <- reshapeMergedReplicates(meanP21, timePoint = 24)
	.chp <- reshapeMergedReplicates(meanCHOP, timePoint = 24)
	.hsp <- reshapeMergedReplicates(meanHSPA1B, timePoint = 24)

	compounds2check <- names(which(table(c(rownames(.bip), rownames(.srx), rownames(.p21), rownames(.chp), rownames(.hsp))) == 5))
		
	dat24 <- data.frame(.bip[compounds2check, ], .srx[compounds2check, 3:20], .p21[compounds2check, 3:20], .chp[compounds2check, 3:20], .hsp[compounds2check, 3:17])

	maxs <- apply(dat24[, grepl('CMAX', colnames(dat24))], 2, max) 
  	mins <- apply(dat24[, grepl('CMAX', colnames(dat24))], 2, min)

  	scaled <- as.data.frame(scale(dat24[, grepl('CMAX', colnames(dat24))], center = mins, scale = maxs - mins))
  	dat24[, grepl('CMAX', colnames(dat24))] <- scaled

  	set.seed(123)
	INDEX <- list()
	for(i in 1:10000) INDEX[[i]] <- .setTrainingIndex(dat24)

	index <- INDEX[[1]]

	train <- dat24[index, ]
  	test  <- dat24[-index, ]

  	f <- as.formula(paste0('binaryDILI ~ ', paste0(grep('CMAX', colnames(dat24), value = T), collapse = ' + ')))

  	n <- names(train)
    fit <- neuralnet(f, data = train, hidden = c(10, 3), linear.output = FALSE,  threshold = 0.05, rep = 10)

    pr <- compute(fit, test[, grep('CMAX', colnames(test))])$net.result
    mse <- sum((test$a_binaryDILI - pr)^2) / nrow(test)
















































	preds24nn_12_4_AUC <- wrapperNNAUC(timeID = 24, iter = 50)




	preds24nn_AUC_3nodes <- wrapperNNAUC(timeID = 24, iter = 100)
	preds48nn_AUC_3nodes <- wrapperNNAUC(timeID = 48, iter = 100)
	preds72nn_AUC_3nodes <- wrapperNNAUC(timeID = 72, iter = 100)

	dev.new()
	par(mfrow = c(1, 3))
	boxplot(preds24nn_AUC_3nodes, ylim = c(0, 1))
		abline(a = .5, b = 0)
	boxplot(preds48nn_AUC_3nodes, ylim = c(0.2, 0.8))
		abline(a = .5, b = 0)
	boxplot(preds72nn_AUC_3nodes, ylim = c(0.2, 0.8))
		abline(a = .5, b = 0)


	
	preds24nn_AUC_3nodes_3layers <- wrapperNNAUC(timeID = 24, iter = 100)
	preds48nn_AUC_3nodes_3layers <- wrapperNNAUC(timeID = 48, iter = 100)
	preds72nn_AUC_3nodes_3layers <- wrapperNNAUC(timeID = 72, iter = 100)

	dev.new()
	par(mfrow = c(1, 3))
	boxplot(preds24nn_AUC_3nodes_3layers, ylim = c(0, 1))
		abline(a = .5, b = 0)
	boxplot(preds48nn_AUC_3nodes_3layers, ylim = c(0.2, 0.8))
		abline(a = .5, b = 0)
	boxplot(preds72nn_AUC_3nodes_3layers, ylim = c(0.2, 0.8))
		abline(a = .5, b = 0)	











	## TO DO
	# 1a) Run prediction using the same index for the 5 reporterlines
	# 1b) Generate composite score whether an iteration predicted a test compound correct
	#	  This composite score should take along what reporterline predicted the compound as DILI (also timepoint to take along)
	# 2) Increase iterations to 10000


	compositeFitter <- function(timeID = 24, iter = 500) {
		.bip <- reshapeMergedReplicates(meanBIP, timePoint = 24)
		.srx <- reshapeMergedReplicates(meanSRXN1, timePoint = 24)
		.p21 <- reshapeMergedReplicates(meanP21, timePoint = 24)
		.chp <- reshapeMergedReplicates(meanCHOP, timePoint = 24)
		.hsp <- reshapeMergedReplicates(meanHSPA1B, timePoint = 24)

		compounds2check <- names(which(table(c(rownames(.bip), rownames(.srx), rownames(.p21), rownames(.chp), rownames(.icm), rownames(.hsp))) == 6))
		
		.bip <- .bip[compounds2check, ]
		.srx <- .srx[compounds2check, ]
		.p21 <- .p21[compounds2check, ]
		.chp <- .chp[compounds2check, ]
		.hsp <- .hsp[compounds2check, ]

		compositePrediction <- NA

		for(i in 1:iter) {
			predBIP <- fitPred(.bip, index = INDEX[[i]], returnIndivCompoundPrediction = TRUE)
			predSRX <- fitPred(.srx, index = INDEX[[i]], returnIndivCompoundPrediction = TRUE)
			predP21 <- fitPred(.p21, index = INDEX[[i]], returnIndivCompoundPrediction = TRUE)
			predCHP <- fitPred(.chp, index = INDEX[[i]], returnIndivCompoundPrediction = TRUE)
			predHSP <- fitPred(.hsp, index = INDEX[[i]], returnIndivCompoundPrediction = TRUE)
		
			rslt <- data.frame(DILI = as.logical(predBIP$class), 
						predBIP$pred, 
						predSRX$pred, 
						predP21$pred, 
						predCHP$pred, 
						predHSP$pred)
			rslt$compositeScore <- rowSums(rslt[, grep('pred', colnames(rslt))]) >= 4
			sideBySide <- rslt[, c('DILI', 'compositeScore')]

			compositePrediction[i] <- length(which(rowSums(sideBySide) == 0 | rowSums(sideBySide) == 2)) / nrow(sideBySide)
		
			
		}

		compositePrediction
	}


	fitPred <- function(dat, returnIndivCompoundPrediction = FALSE, index = NULL) {
		if(FALSE) dat <- .bip
		if(FALSE) index <- INDEX[[1]]

		maxs <- apply(dat[, grepl('CMAX', colnames(dat))], 2, max) 
  		mins <- apply(dat[, grepl('CMAX', colnames(dat))], 2, min)

  		scaled <- as.data.frame(scale(dat[, grepl('CMAX', colnames(dat))], center = mins, scale = maxs - mins))
  		dat[, grepl('CMAX', colnames(dat))] <- scaled

  		colnames(dat) <- paste0('a_', colnames(dat))

  		if(is.null(index)) index <- sort(.setTrainingIndex(dat))

  		train <- dat[index, ]
  		test  <- dat[-index, ]

  	
  		f <- as.formula(paste0('a_binaryDILI ~ ', paste0(grep('CMAX', colnames(dat), value = T), collapse = ' + ')))

  		fit <- glm(f, data = train, family = binomial(logit), maxit = 1000)
  		pr  <- predict(fit, test)
  		MSE <- sum((pr - test$a_binaryDILI)^2)/nrow(test)

  			sideBySide <- data.frame(pred = round(pr) > 0, class = test$a_binaryDILI)
  		
  		fractionProperlyPredicted <- length(which(rowSums(sideBySide) == 0 | rowSums(sideBySide) == 2)) / nrow(sideBySide)

  		if(returnIndivCompoundPrediction) return(sideBySide)
  		
  		return(fractionProperlyPredicted)
	}

	wrapperFitter <- function(timeID = 24, iter = 500) {
		.bip <- reshapeMergedReplicates(meanBIP, timePoint = 24)
		.srx <- reshapeMergedReplicates(meanSRXN1, timePoint = 24)
		.p21 <- reshapeMergedReplicates(meanP21, timePoint = 24)
		.chp <- reshapeMergedReplicates(meanCHOP, timePoint = 24)
		.hsp <- reshapeMergedReplicates(meanHSPA1B, timePoint = 24)

		compounds2check <- names(which(table(c(rownames(.bip), rownames(.srx), rownames(.p21), rownames(.chp), rownames(.icm), rownames(.hsp))) == 6))
		
		.bip <- .bip[compounds2check, ]
		.srx <- .srx[compounds2check, ]
		.p21 <- .p21[compounds2check, ]
		.chp <- .chp[compounds2check, ]
		.hsp <- .hsp[compounds2check, ]

		predictionsReporter <- data.frame(bip = NA, srx = NA, p21 = NA, chp = NA, hsp = NA)

		for(i in 1:iter) {
			predictionsReporter[i, 'bip'] <- fitPred(.bip, index = INDEX[[i]])
			predictionsReporter[i, 'srx'] <- fitPred(.srx, index = INDEX[[i]])
			predictionsReporter[i, 'p21'] <- fitPred(.p21, index = INDEX[[i]])
			predictionsReporter[i, 'chp'] <- fitPred(.chp, index = INDEX[[i]])
			predictionsReporter[i, 'hsp'] <- fitPred(.hsp, index = INDEX[[i]])
		}

		predictionsReporter
	}



	preds24Composite <- compositeFitter()


	predsBIP <- wrapperFitter(meanBIP)
	predsP21 <- wrapperFitter(meanP21)
	predsCHOP <- wrapperFitter(meanCHOP)
	predsHSPA1B <- wrapperFitter(meanHSPA1B)
	predsSRXN1 <- wrapperFitter(meanSRXN1)






	