## TO DO ##
  	# 1) CMAX concentration as input parameter
  	# 2) independent variable selection
  	# 5) GET AIC/NIC/BIC WORKING FFS! (model evaluation)
  	# 6) ICAM1 can go up AND down compared to DMSO (due to addition of TNF at 24H)

## Settings
	options(stringsAsFactors = FALSE)
	
	library(neuralnet)
	library(nnetpredint)
	library(NeuralNetTools)
	library(caret)

	library(ggplot2)

	library(data.table)
	library(preprocessCore)
	library(reshape2)
	library(parallel)

	inputDir <- '/Users/Wouter/stack/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'


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
		features2merge <- c(grep('_Intensity_', colnames(sumDatas[[1]]), value = TRUE),
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

	mergeReplicates_fractionDown <- function(parsedList, .debug = FALSE) {
		if(.debug) {
			parsedList <- ICAM1
		}

		# Parseing
		sumDatas <- lapply(parsedList, function(rep) {
			tmp <- as.data.frame(rep$sumData_rev)
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
		features2merge <- c(grep('_Intensity_', colnames(sumDatas[[1]]), value = TRUE),
							'numberOfObjects', 'plateMeanDMSO', 'plateMedianDMSO', 
							grep('_Below_', colnames(sumDatas[[1]]), value = TRUE),
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

		return(rslt)		# ICAM1
	}


	reshapeMergedReplicates_DMSO <- function(meanReporter, reporterName, features = c('CMAX_numeric', 'Fraction_Above_1_medianIntregated_plateDMSO', 'Fraction_Above_2_medianIntregated_plateDMSO', 'Fraction_Above_3_medianIntregated_plateDMSO', 'Fraction_Above_5_medianIntregated_plateDMSO'), timePoint, .debug = FALSE) {
		if(.debug) {
			meanReporter <- meanICAM1
			features <- c('CMAX_numeric', 'Fraction_Above_1_medianIntregated_plateDMSO', 'Fraction_Above_2_medianIntregated_plateDMSO', 'Fraction_Above_3_medianIntregated_plateDMSO', 'Fraction_Above_5_medianIntregated_plateDMSO')
			timePoint <- 48
		}

		rslt <- DILI[, c('treatment', 'binaryDILI')]
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_Above_1_medianIntregated_plateDMSO')] <- NA 
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_Above_2_medianIntregated_plateDMSO')] <- NA
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_Above_3_medianIntregated_plateDMSO')] <- NA
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_Above_5_medianIntregated_plateDMSO')] <- NA
				rownames(rslt) <- rslt$treatment

		if(reporterName == 'ICAM1') {
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_Below_1_medianIntregated_plateDMSO')] <- NA
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_Below_2_medianIntregated_plateDMSO')] <- NA
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_Below_3_medianIntregated_plateDMSO')] <- NA
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_Below_5_medianIntregated_plateDMSO')] <- NA

			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_netchange_1_medianIntegrated_plateDMSO')] <- NA
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_netchange_2_medianIntegrated_plateDMSO')] <- NA
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_netchange_3_medianIntegrated_plateDMSO')] <- NA
			rslt[, paste0('CMAX', c(1, 5, 10, 25, 50, 100), '_Fraction_netchange_5_medianIntegrated_plateDMSO')] <- NA
		}		
		
		for(trtID in rownames(rslt)) {
			if(length(which(meanReporter$treatment == trtID)) == 0) next
			
			for(ftr in grep('CMAX', colnames(rslt), value = TRUE)) {
				cmax     <- as.numeric(gsub('CMAX', '', unlist(strsplit(ftr, '_'))[1]))
				
				.ind <- which(meanReporter$treatment == trtID & meanReporter$timeID == timePoint & meanReporter$CMAX_numeric == cmax)
				.ftr <- paste(unlist(strsplit(ftr, '_'))[-1], collapse = '_')

				if(length(.ind) == 0) next

				rslt[trtID, ftr] <- meanReporter[.ind, .ftr]

			}
		}

		rslt <- rslt[,colSums(is.na(rslt)) != nrow(rslt)]
		rslt <- rslt[which(rowSums(is.na(rslt)) == 0), ]

		rslt$reporterLine <- reporterName

		rslt
	}

	setTrainingIndex <- function(data, class = c('binaryDILI', 'doseDILI'), negN = 40, posN = 80) {
		if(length(class) != 1) stop('Only supply a single column name as dependent variable.\n')
  		
  		neg <- which(data[, class] == 0)
  		pos <- which(data[, class] == 1)

  		trainingNegative <- sample(neg, negN)
  		trainingPositive <- sample(pos, posN)

  		c(trainingNegative, trainingPositive)
  	}


## Screening data
	if(TRUE) {
		load(paste0(inputDir, 'parsedBIP.RData'))
		meanBIP <- mergeReplicates(BIP)
		meanBIP$treatment <- gsub('SODIUM-ARSENITE', 'SODIUMARSENITE', gsub('STAUSPORIN', 'STAUROSPORIN', gsub('BUTHIONINESULFOXIMINE', 'BUTHIONINESULFOXAMINE', meanBIP$treatment)))
			save(meanBIP, file = paste0(inputDir, 'mergedReplicates_sumDat_BIP.RData'))
			 rm(BIP)

		load(paste0(inputDir, 'parsedSRXN1.RData'))
		meanSRXN1 <- mergeReplicates(SRXN1)	
			save(meanSRXN1, file = paste0(inputDir, 'mergedReplicates_sumDat_SRXN1.RData'))
			 rm(SRXN1)

		load(paste0(inputDir, 'parsedP21.RData'))
		meanP21 <- mergeReplicates(P21)	
			save(meanP21, file = paste0(inputDir, 'mergedReplicates_sumDat_P21.RData'))
			 rm(P21)

		load(paste0(inputDir, 'parsedCHOP.RData'))
		meanCHOP <- mergeReplicates(CHOP)	
			save(meanCHOP, file = paste0(inputDir, 'mergedReplicates_sumDat_CHOP.RData'))
			 rm(CHOP)

		load(paste0(inputDir, 'parsedICAM1.RData'))
		meanICAM1 <- mergeReplicates(ICAM1)	
		meanICAM1_down <- mergeReplicates_fractionDown(ICAM1)	
			if(identical(rownames(meanICAM1), rownames(meanICAM1_down))) meanICAM1 <- cbind(meanICAM1, meanICAM1_down[, grep('Below', colnames(meanICAM1_down), value = TRUE)])
			for(i in c(1, 2, 3, 5)) {
				.netChange_median <- meanICAM1[, grep(paste0('Fraction_Above_', i, '_median'), colnames(meanICAM1), value = TRUE)] - meanICAM1[, grep(paste0('Fraction_Below_', i, '_median'), colnames(meanICAM1), value = TRUE)]
				.netChange_mean   <- meanICAM1[, grep(paste0('Fraction_Above_', i, '_mean'), colnames(meanICAM1), value = TRUE)] - meanICAM1[, grep(paste0('Fraction_Below_', i, '_mean'), colnames(meanICAM1), value = TRUE)]
			
				meanICAM1[, paste0('Fraction_netchange_', i, '_medianIntegrated_plateDMSO')] <- .netChange_median
				meanICAM1[, paste0('Fraction_netchange_', i, '_meanIntegrated_plateDMSO')]   <- .netChange_mean
			}

			save(meanICAM1, file = paste0(inputDir, 'mergedReplicates_sumDat_ICAM1.RData'))
			 rm(ICAM1)

		load(paste0(inputDir, 'parsedHSPA1B.RData'))
		meanHSPA1B <- mergeReplicates(HSPA1B)
		meanHSPA1B$treatment <- gsub('SODIUM-ARSENITE', 'SODIUMARSENITE', gsub('STAUSPORIN', 'STAUROSPORIN', gsub('BUTHIONINESULFOXIMINE', 'BUTHIONINESULFOXAMINE', meanHSPA1B$treatment)))
			save(meanHSPA1B, file = paste0(inputDir, 'mergedReplicates_sumDat_HSPA1B.RData'))
			 rm(HSPA1B)	

		load(paste0(inputDir, 'parsedBTG2.RData'))
		meanBTG2 <- mergeReplicates(BTG2)	
			save(meanBTG2, file = paste0(inputDir, 'mergedReplicates_sumDat_BTG2.RData'))
			 rm(BTG2)	

		load(paste0(inputDir, 'parsedHMOX1.RData'))
		meanHMOX1 <- mergeReplicates(HMOX1)	
			save(meanHMOX1, file = paste0(inputDir, 'mergedReplicates_sumDat_HMOX1.RData'))
			 rm(HMOX1)				


	}

	#load(paste0(inputDir, 'mergedReplicates_sumDat_BIP.RData'))
	#load(paste0(inputDir, 'mergedReplicates_sumDat_SRXN1.RData'))
	#load(paste0(inputDir, 'mergedReplicates_sumDat_P21.RData'))
	#load(paste0(inputDir, 'mergedReplicates_sumDat_CHOP.RData'))
	#load(paste0(inputDir, 'mergedReplicates_sumDat_ICAM1.RData'))
	#load(paste0(inputDir, 'mergedReplicates_sumDat_HSPA1B.RData'))
	#load(paste0(inputDir, 'mergedReplicates_sumDat_BTG2.RData'))
	#load(paste0(inputDir, 'mergedReplicates_sumDat_HMOX1.RData'))
	
	if(TRUE) {		# Save summData's to tab delim for Lhasa data deposit (mail Steven Hiemstra)
		save2text <- function(df, cols = c('treatment', 'timeID', 'dose_uM', 'CMAX_numeric', 'Fraction_Above_2_meanIntregated_plateDMSO', 'PI_fractionNonNA', grep('_Intensity_', colnames(df), value = TRUE)), file) {
			df2save <- df[, cols]
				colnames(df2save)[2] <- 'timepoint'
				colnames(df2save)[4] <- 'CMAX'
				colnames(df2save)[5] <- 'fraction_positiveGFPcells_2xDMSO'
				colnames(df2save)[6] <- 'fraction_positivePIcells'
				colnames(df2save)[7] <- 'meanIntensity_GFP'

			write.table(df2save, file = file, sep = '\t', row.names = FALSE, quote = FALSE)
		}

		save2text(df = meanBIP,    file = paste0(inputDir, 'sumDat_BIP.txt'))
		save2text(df = meanSRXN1,  file = paste0(inputDir, 'sumDat_SRXN1.txt'))
		save2text(df = meanP21,    file = paste0(inputDir, 'sumDat_P21.txt'))
		save2text(df = meanCHOP,   file = paste0(inputDir, 'sumDat_CHOP.txt'))
		save2text(df = meanICAM1,  file = paste0(inputDir, 'sumDat_ICAM1.txt'))
		save2text(df = meanHSPA1B, file = paste0(inputDir, 'sumDat_HSPA1B.txt'))
		save2text(df = meanBTG2,   file = paste0(inputDir, 'sumDat_BTG2.txt'))
		save2text(df = meanHMOX1,  file = paste0(inputDir, 'sumDat_HMOX1.txt'))
	
	}

	# TOX PPT Monday 20 Feb 2017
	if(FALSE) {
		# Timeline plot for experiment
		library(ggplot2)

		expDesign <- data.frame(timepoint = c(0, 24, 48, 72), 
									   event = c('Compound', 'Imaging', 'Imaging', 'Imaging'))

		ggplot(expDesign, aes(x = timepoint, y = 0)) + # + ylim(-1, 1)
			geom_point(size = 5) +
			scale_x_continuous(breaks = c(0, 24, 48, 72), limits = c(-5, 77)) +
			geom_segment(aes(x = 0, y = -0.1, xend = 0, yend = 0.1), size = 1.5) +
			geom_segment(aes(x = 24, y = -0.1, xend = 24, yend = 0.1), size = 1.5) +
			geom_segment(aes(x = 48, y = -0.1, xend = 48, yend = 0.1), size = 1.5) +
			geom_segment(aes(x = 72, y = -0.1, xend = 72, yend = 0.1), size = 1.5) +
			
			geom_segment(aes(x = -2, y = 0, xend = 74, yend = 0), size = 1) +

			ylim(-0.5, 0.5) 


			geom_segment(aes(y=0,yend=y,xend=year))

	}


## DILI annotation data
	if(FALSE) {	## Old annotation by Steven Hiemstra
		DILI <- read.delim('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/DILIscreen/DILI_annotation.txt')

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
	}


## newDILI annotation data
	DILI <- read.delim('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/DILIscreen/DILI_annotation_dissolved_10012017.txt')  
	DILI <- DILI[which(DILI$treatment != ''), ]
	
		DILI$treatment <- gsub(' ', '', toupper(DILI$treatment))
	  	  DILI$treatment[DILI$treatment == "CARBAMAZEPINE"]         <- "CARBAMAZAPINE"
	  	  DILI$treatment[DILI$treatment == "G.F.DMEM"]         		<- "DMEM"

	  	DILI[nrow(DILI) + 1, ] <- DILI[which(DILI$treatment == 'TUNICAMYCIN'), ]
	  	DILI[nrow(DILI), 'treatment'] <- 'TUNICAMYCINC'

	  	DILI[nrow(DILI) + 1, ] <- DILI[which(DILI$treatment == 'ETOPOSIDE'), ]
	  	DILI[nrow(DILI), 'treatment'] <- 'ETOPOSIDEC'

	DILI$binaryDILI <- DILI$doseDILI <- as.numeric(grepl('Most', DILI$DILIConcern) | grepl('Less', DILI$DILIConcern))
	DILI$doseDILI[grepl('Most', DILI$DILIConcern)] <- 2 

	rownames(DILI) <- DILI$treatment


## Reshape reporter data according to DILI annotation
	annotatedBIP <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates_DMSO(meanBIP,    ,reporterName = 'BIP',    timePoint = .timePoint))
	annotatedSRX <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates_DMSO(meanSRXN1,  ,reporterName = 'SRXN1',  timePoint = .timePoint))
	annotatedP21 <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates_DMSO(meanP21,    ,reporterName = 'P21',    timePoint = .timePoint))
	annotatedCHP <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates_DMSO(meanCHOP,   ,reporterName = 'CHOP',   timePoint = .timePoint))
	annotatedHSP <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates_DMSO(meanHSPA1B, ,reporterName = 'HSPA1B', timePoint = .timePoint))
	annotatedICM <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates_DMSO(meanICAM1,  ,reporterName = 'ICAM1',  timePoint = .timePoint))	# misses 24 hour data
	annotatedBTG <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates_DMSO(meanBTG2,   ,reporterName = 'BTG2',   timePoint = .timePoint))
	annotatedHMX <- lapply(c(24, 48, 72), function(.timePoint) reshapeMergedReplicates_DMSO(meanHMOX1,  ,reporterName = 'HMOX1',  timePoint = .timePoint))
		names(annotatedBIP) <- names(annotatedSRX) <- names(annotatedP21) <- names(annotatedCHP) <- names(annotatedHSP) <- names(annotatedICM) <- names(annotatedBTG) <- names(annotatedHMX) <- c('h24', 'h48', 'h72') 

	annotatedAll <- list(BIP    = annotatedBIP, 
						 SRXN1  = annotatedSRX, 
						 P21    = annotatedP21, 
						 CHOP   = annotatedCHP, 
						 HSPA1B = annotatedHSP, 
						 ICAM1  = annotatedICM,
						 BTG2   = annotatedBTG,
						 HMOX1  = annotatedHMX) # β ε


## Create objects to input into classifiers
  ## Combine reporters per timepoint
  # 24 hour:	# 154 treatments in all but ICAM1: 108 DILI / 46 nonDILI
  	treatments_24h <- names(which(table(unlist(lapply(annotatedAll[-which(names(annotatedAll) == 'ICAM1')], function(x) { rownames(x$h24) } ))) == 7))
  	classifierInput_24 <- data.frame(binaryDILI = DILI[treatments_24h, 'binaryDILI'],
  									 doseDILI = DILI[treatments_24h, 'doseDILI'],
									 lapply(annotatedAll[-which(names(annotatedAll) == 'ICAM1')], function(x) {
									   tmp <- x$h24[treatments_24h, 3:ncol(x$h24)]
									   tmp[, -grep('reporterLine', colnames(tmp))]
							  		 }))

  # 48 hour:	# 171 treatments: 112 DILI / 59 nonDILI
  	treatments_48h <- names(which(table(unlist(lapply(annotatedAll, function(x) { rownames(x$h48) } ))) == 8))
  	classifierInput_48 <- data.frame(binaryDILI = DILI[treatments_48h, 'binaryDILI'],
  									 doseDILI = DILI[treatments_48h, 'doseDILI'],
									 lapply(annotatedAll, function(x) {
									   tmp <- x$h48[treatments_48h, 3:ncol(x$h48)]
									   tmp[, -grep('reporterLine', colnames(tmp))]
							  		 }))

  # 72 hour:	# 177 treatments in all but ICAM1: 117 DILI / 60 nonDILI
  	treatments_72h <- names(which(table(unlist(lapply(annotatedAll, function(x) { rownames(x$h72) } ))) == 8))
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

  # BTG2 		# 178 treatments: 118 DILI / 60 nonDILI
    treatments_BTG2 <- names(which(table(unlist(lapply(annotatedBTG, rownames))) == 3))
    classifierInput_BTG2 <- data.frame(binaryDILI = DILI[treatments_BTG2, 'binaryDILI'],
    									doseDILI = DILI[treatments_BTG2, 'doseDILI'],
    									lapply(annotatedBTG, function(x) {
    									  x[treatments_BTG2, 3:(ncol(x)-1)]
    									}))

  # HMOX1 		# 178 treatments: 118 DILI / 60 nonDILI
    treatments_HMOX1 <- names(which(table(unlist(lapply(annotatedHMX, rownames))) == 3))
    classifierInput_HMOX1 <- data.frame(binaryDILI = DILI[treatments_HMOX1, 'binaryDILI'],
    									doseDILI = DILI[treatments_HMOX1, 'doseDILI'],
    									lapply(annotatedHMX, function(x) {
    									  x[treatments_HMOX1, 3:(ncol(x)-1)]
    									}))


## Neuralnetting
 # 24 hour dataset
  	set.seed(1337)

  	INDEX <- list()
	for(i in 1:1000) INDEX[[i]] <- setTrainingIndex(classifierInput_24, 'binaryDILI')


	fitNN <- function(dat, NodesLayers = c(10, 3), splitBy = c('binaryDILI', 'doseDILI'), inputParameterGrep = c('CMAX', '_1', '_2', '_3', '_5'), index = NULL, nCores, .debug = FALSE) {
	  # TO DO
	  	# 3) nnetpredint implementation for CIs
	  	# 4) Compute confidence intervals using SEM
	  	# 5) GET AIC/NIC/BIC WORKING FFS! (model evaluation)
	  	#		in that regard: histogram of specificity/sensitivity only contains single values (not continious)
	  	#		Not sure if this is problem; might be due to verdeling dili/non-dili compounds in training and test set OR INDEX[1]
	  	#		

	  # Debugging	
		if(.debug) {
			NodesLayers <- c(10, 3)
			splitBy 	<- 'binaryDILI'
			index       <- INDEX[[1]]

			#dat           <- classifierInput_BIP[, c(1,2, grep('DMSO', colnames(classifierInput_BIP)))]
			dat 		  <- classifierInput_ALL[compoundBob, ]
			#dat$doseDILI  <- as.numeric(dat$doseDILI == 2)
			colnames(dat) <- gsub('_medianIntregated_plateDMSO', '', gsub('Fraction_Above_', '', colnames(dat)))
		}


		# Data split
  	  	if(is.null(index)) index <- sort(setTrainingIndex(dat, splitBy))
	
  	  	train <- dat[index,  ]
  	  	test  <- dat[-index, ]

  	  # Training neural network
  	  	f <- as.formula(paste0('binaryDILI ~ ', paste0(grep(inputParameterGrep, colnames(dat), value = T), collapse = ' + ')))
    	fitList <- mclapply(1:10, function(i) { 
    			cat(as.character(Sys.time()), '    Initiating thread', i, '\n')
    				subThreatFit <- neuralnet(formula = f, 
    				data = train, 
    				hidden = NodesLayers, 
    				stepmax = 5e+05,
    				linear.output = FALSE, 
    				threshold = 0.001, 
    				err.fct = 'ce',
    				rep = 1,
    				likelihood = TRUE) 
    			cat(as.character(Sys.time()), '    Finished thread', i, '\n')

    			return(subThreatFit)
    		   } , mc.cores = nCores)

    	.combineSubThreatFits <- function(fitList) {
    		rslt <- fitList[[1]]

    		rslt$net.result          <- do.call(c, lapply(fitList,     function(subThreatFit) { subThreatFit$net.result          } ))
    		rslt$weights             <- do.call(c, lapply(fitList,     function(subThreatFit) { subThreatFit$weights             } ))
    		rslt$startweights        <- do.call(c, lapply(fitList,     function(subThreatFit) { subThreatFit$startweights        } ))
    		rslt$generalized.weights <- do.call(c, lapply(fitList,     function(subThreatFit) { subThreatFit$generalized.weights } ))
    		rslt$result.matrix       <- do.call(cbind, lapply(fitList, function(subThreatFit) { subThreatFit$result.matrix       } ))

    		rslt
    	}

    	fit <- .combineSubThreatFits(fitList)

    	if(FALSE) plotnet(fitList[[1]], alpha = 0.25, pos_col = 'steelblue3', neg_col = 'orangered3', cex = 0.75)

      # Evaulating the trained networks
    	pr <- lapply(1:length(fit$net.result), function(.rep) { compute(fit, test[, grep(inputParameterGrep, colnames(test))], rep = .rep)$net.result  } )

    	.calcAUC <- function(truestat, testres) {
    		# http://stats.stackexchange.com/questions/145566/how-to-calculate-area-under-the-curve-auc-or-the-c-statistic-by-hand

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
		
			rslt <- sum(height*width)	# AUC
    	}

    	AUCs <- sapply(pr, function(predDILI) .calcAUC(truestat = round(test$binaryDILI), testres = predDILI))
    	spec <- sapply(1:length(pr), function(i) specificity(as.factor(round(pr[[i]])[,1]), as.factor(test$binaryDILI)))
    	sens <- sapply(1:length(pr), function(i) sensitivity(as.factor(round(pr[[i]])[,1]), as.factor(test$binaryDILI)))
    	
      # return values
      	rslt <- list(train = train, test = test, call = f, trainingFit = fit, testPred = pr, AUCs = AUCs, spec = spec, sens = sens)

      	return(rslt)
	}

  ## Save environment 
	save.image(file = paste0(inputDir, 'nnInput_28Feb2017.RData'))	


## heatmap for Bob
	load(paste0(inputDir, 'nnInput_28Feb2017.RData'))

	library(pheatmap)

	createPheatmapData <- function(classifierInput, fractionAbove = c(1, 2, 3, 5), .debug = FALSE) {
		if(.debug) {
			classifierInput <- classifierInput_48
			fractionAbove   <- 3
		}

		classifierInput[, c(1, 2, grep(paste0('_', fractionAbove), colnames(classifierInput)))]
	}

	plotPheatmap <- function(pheatmapData, title = '', treatmentAnnot = c('binary', 'dose', 'severe')) {
		d <- pheatmapData[, 3:ncol(pheatmapData)]
		  colsICAM1 <-  colnames(d)[grepl('ICAM1', colnames(d)) & grepl('netchange', colnames(d))]
		  colsOthers <- colnames(d)[!grepl('ICAM1', colnames(d))]

		d <- d[, c(colsOthers, colsICAM1)]

		if(treatmentAnnot == 'binary') {
			annotRow <- pheatmapData[, 'binaryDILI', drop = FALSE]
			colnames(annotRow) <- 'DILI'
		}

		if(treatmentAnnot == 'dose') {
			annotRow <- data.frame(DILI = pheatmapData$doseDILI)
			rownames(annotRow) <- rownames(d)
		}

		if(treatmentAnnot == 'severe') {
			annotRow <- data.frame(DILI = as.numeric(pheatmapData$doseDILI == 2))
			rownames(annotRow) <- rownames(d)
		}

		annotCol <- data.frame(CMAX = as.numeric(factor(as.numeric(gsub('CMAX', '', do.call(rbind, strsplit(do.call(rbind, strsplit(colnames(d), '\\.'))[,2], '_'))[,1])), ordered = TRUE)))
			rownames(annotCol) <- colnames(d)

		breakMax <- max(abs(d))
		myBreaks <- c(seq((-1 * breakMax), 0, length.out = ceiling(100/2) + 1), 
		          		  seq(breakMax/100, breakMax, length.out = floor(100/2)))

		pheatmap(d, 
			 clustering_distance_rows = 'correlation',
			 cluster_cols = FALSE, cluster_rows = TRUE,
			 breaks = myBreaks, 
			 legend = FALSE, 
			 border_color = NA, 
			 annotation_row = annotRow,
			 annotation_col = annotCol,
			 annotation_legend = FALSE,
			 labels_col = '', 
			 fontsize_row = 3,
			 main = title)
	}

	compoundBob <- rownames(DILI[DILI$DILIConcern != 'N/A', ])


	# Timepoint datasets
	hm_24H <- lapply(c(1, 2, 3, 5), function(i) { createPheatmapData(classifierInput = classifierInput_24, fractionAbove = i) } )
	hm_48H <- lapply(c(1, 2, 3, 5), function(i) { createPheatmapData(classifierInput = classifierInput_48, fractionAbove = i) } )
	hm_72H <- lapply(c(1, 2, 3, 5), function(i) { createPheatmapData(classifierInput = classifierInput_72, fractionAbove = i) } )
		names(hm_24H) <- names(hm_48H) <- names(hm_72H) <- c('DMSO_1', 'DMSO_2', 'DMSO_3', 'DMSO_5')

	pdf(paste0(inputDir, '_heatmap_24H_doseAnnotated.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_24H[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_24H[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_24H[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_24H[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'dose')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_48H_doseAnnotated.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_48H[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_48H[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_48H[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_48H[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'dose')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_72H_doseAnnotated.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_72H[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_72H[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_72H[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'dose')
		plotPheatmap(na.omit(hm_72H[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'dose')
	dev.off()


	# Reporter datasets
	hm_BIP    <- lapply(c(1, 2, 3, 5), function(dmsoCutoff) classifierInput_BIP[, c(1, 2, grep(paste0('Above_', dmsoCutoff), colnames(classifierInput_BIP)))] )
	hm_SRXN1  <- lapply(c(1, 2, 3, 5), function(dmsoCutoff) classifierInput_SRXN1[, c(1, 2, grep(paste0('Above_', dmsoCutoff), colnames(classifierInput_SRXN1)))] )
	hm_P21    <- lapply(c(1, 2, 3, 5), function(dmsoCutoff) classifierInput_P21[, c(1, 2, grep(paste0('Above_', dmsoCutoff), colnames(classifierInput_P21)))] ) 
	hm_CHOP   <- lapply(c(1, 2, 3, 5), function(dmsoCutoff) classifierInput_CHOP[, c(1, 2, grep(paste0('Above_', dmsoCutoff), colnames(classifierInput_CHOP)))] )
	hm_HSPA1B <- lapply(c(1, 2, 3, 5), function(dmsoCutoff) classifierInput_HSPA1B[, c(1, 2, grep(paste0('Above_', dmsoCutoff), colnames(classifierInput_HSPA1B)))] )
	hm_ICAM1  <- lapply(c(1, 2, 3, 5), function(dmsoCutoff) classifierInput_ICAM1[, c(1, 2, grep(paste0('Above_', dmsoCutoff), colnames(classifierInput_ICAM1)))] )
	hm_BTG2   <- lapply(c(1, 2, 3, 5), function(dmsoCutoff) classifierInput_BTG2[, c(1, 2, grep(paste0('Above_', dmsoCutoff), colnames(classifierInput_BTG2)))] )
	hm_HMOX1  <- lapply(c(1, 2, 3, 5), function(dmsoCutoff) classifierInput_HMOX1[, c(1, 2, grep(paste0('Above_', dmsoCutoff), colnames(classifierInput_HMOX1)))] )
	
	pdf(paste0(inputDir, '_heatmap_BIP.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_BIP[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_BIP[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_BIP[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_BIP[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_SRXN1.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_SRXN1[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_SRXN1[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_SRXN1[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_SRXN1[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_P21.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_P21[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_P21[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_P21[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_P21[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_CHOP.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_CHOP[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_CHOP[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_CHOP[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_CHOP[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_HSPA1B.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_HSPA1B[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_HSPA1B[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_HSPA1B[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_HSPA1B[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_ICAM1.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_ICAM1[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_ICAM1[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_ICAM1[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_ICAM1[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_BTG2.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_BTG2[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_BTG2[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_BTG2[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_BTG2[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, '_heatmap_HMOX1.pdf'), width = 5, height = 10)
		plotPheatmap(na.omit(hm_HMOX1[[1]][compoundBob, ]), title = 'GFP positive fraction of cells >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_HMOX1[[2]][compoundBob, ]), title = 'GFP positive fraction of cells >2x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_HMOX1[[3]][compoundBob, ]), title = 'GFP positive fraction of cells >3x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_HMOX1[[4]][compoundBob, ]), title = 'GFP positive fraction of cells >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	# Hypergeometric enrichment
	hgEnrichment <- function(hm_data, dmsoCutoff = c('DMSO_1', 'DMSO_2', 'DMSO_3', 'DMSO_5'), clusters = 5) {
		# hm_data <- hm_48H
		# dmsoCutoff <- 'DMSO_2'
		# annot <- 'severe'

		dat <- hm_data[[dmsoCutoff]]

		d   <- dat[, 3:ncol(dat)]
		annotRow <- data.frame(DILI = as.numeric(dat$doseDILI == 2))
			rownames(annotRow) <- rownames(d)

		ph <- pheatmap(d, clustering_distance_rows = 'correlation', cluster_cols = FALSE, cluster_rows = TRUE)
		annotRow$cl <- paste0('cl', cutree(ph$tree_row, clusters))


		# Hypergeometric enrichment test
		totalTreat   <- rownames(annotRow)
		totalTreat.L <- length(totalTreat)

		diliTreat   <- totalTreat[as.logical(annotRow$DILI)]
		diliTreat.L <- length(diliTreat)

		clDILI   <- lapply(unique(annotRow$cl), function(.cl) { rownames(annotRow[which(annotRow$cl == .cl), ]) } )
		clDILI.L <- sapply(clDILI, length)

		overlap_diliTreat_clDILI   <- lapply(clDILI, function(x) x[which(x %in% diliTreat)])
		overlap_diliTreat_clDILI.L <- sapply(overlap_diliTreat_clDILI, length)

		hg_test <- sapply(1:clusters, function(i) {
		phyper(overlap_diliTreat_clDILI.L[i],
			   diliTreat.L,
			   totalTreat.L - diliTreat.L,
			   clDILI.L[i],
			   lower.tail = FALSE, log.p = FALSE)
		})
		names(hg_test) <- paste0('cl', 1:clusters)

		rslt <- list(ph, annotRow, hg_test)
		names(rslt) <- c('ph', 'annot', 'hg_test')

		rslt
		
	}

	hg_hm_24H <- hgEnrichment(hm_24H, dmsoCutoff = 'DMSO_3', clusters = 8)
	hg_hm_24H$clOrder <- hg_hm_24H$ph$tree_row$labels[hg_hm_24H$ph$tree_row$order]
	hg_hm_24H$hg_test_ordered <- data.frame(hg_hm_24H$annot[hg_hm_24H$clOrder, ], hg_hm_24H$hg_test[hg_hm_24H$annot[hg_hm_24H$clOrder, 'cl']])

	hg_hm_48H <- hgEnrichment(hm_48H, 'DMSO_2')
	hg_hm_48H$clOrder <- hg_hm_48H$ph$tree_row$labels[hg_hm_48H$ph$tree_row$order]
	hg_hm_48H$hg_test_ordered <- data.frame(hg_hm_48H$annot[hg_hm_48H$clOrder, ], hg_hm_48H$hg_test[hg_hm_48H$annot[hg_hm_48H$clOrder, 'cl']])

	hg_hm_72H <- hgEnrichment(hm_72H, 'DMSO_3')
	hg_hm_72H$clOrder <- hg_hm_72H$ph$tree_row$labels[hg_hm_72H$ph$tree_row$order]
	hg_hm_72H$hg_test_ordered <- data.frame(hg_hm_72H$annot[hg_hm_72H$clOrder, ], hg_hm_72H$hg_test[hg_hm_72H$annot[hg_hm_72H$clOrder, 'cl']])

## heatmap for DILI bookchapter
	createPheatmapData <- function(classifierInput, fractionAbove = c(1, 2, 3, 5), .debug = FALSE) {
		if(.debug) {
			classifierInput <- classifierInput_48
			fractionAbove   <- 3
		}

		rslt <- classifierInput[, c(1, 2, grep(paste0('_', fractionAbove), colnames(classifierInput)))]
			allColumns_ButICAM1   <- grep('ICAM1', colnames(rslt), value = TRUE, invert = TRUE)
			ICAM1Columns_selected <- grep('netchange', grep('ICAM1', colnames(rslt), value = TRUE), value = TRUE)
	
		rslt[, c(allColumns_ButICAM1, ICAM1Columns_selected)]
	}

	hm_24H <- lapply(c(1, 2, 3, 5), function(i) { createPheatmapData(classifierInput = classifierInput_24, fractionAbove = i) } )
	hm_48H <- lapply(c(1, 2, 3, 5), function(i) { createPheatmapData(classifierInput = classifierInput_48, fractionAbove = i) } )
	hm_72H <- lapply(c(1, 2, 3, 5), function(i) { createPheatmapData(classifierInput = classifierInput_72, fractionAbove = i) } )
		names(hm_24H) <- names(hm_48H) <- names(hm_72H) <- c('DMSO_1', 'DMSO_2', 'DMSO_3', 'DMSO_5')

	plotPheatmap <- function(pheatmapData, title = '', treatmentAnnot = c('binary', 'dose', 'severe')) {
		d <- pheatmapData[, 3:ncol(pheatmapData)]

		if(treatmentAnnot == 'binary') {
			annotRow <- pheatmapData[, 'binaryDILI', drop = FALSE]
			colnames(annotRow) <- 'DILI'
		}

		if(treatmentAnnot == 'dose') {
			annotRow <- data.frame(DILI = pheatmapData$doseDILI)
			rownames(annotRow) <- rownames(d)
		}

		if(treatmentAnnot == 'severe') {
			annotRow <- data.frame(DILI = as.numeric(pheatmapData$doseDILI == 2))
			rownames(annotRow) <- rownames(d)
		}

		annotCol <- data.frame(CMAX = as.numeric(factor(as.numeric(gsub('CMAX', '', do.call(rbind, strsplit(do.call(rbind, strsplit(colnames(d), '\\.'))[,2], '_'))[,1])), ordered = TRUE)))
			rownames(annotCol) <- colnames(d)

		# breakMax <- max(abs(d))
		breakMax <- 1
		myBreaks <- c(seq((-1 * breakMax), 0, length.out = ceiling(100/2) + 1), 
		          		  seq(breakMax/100, breakMax, length.out = floor(100/2)))

		pheatmap(d, 
			 clustering_distance_rows = 'correlation',
			 cluster_cols = FALSE, cluster_rows = TRUE,
			 breaks = myBreaks, 
			 legend = TRUE, 
			 border_color = NA, 
			 annotation_row = annotRow,
			 annotation_col = annotCol,
			 annotation_legend = FALSE,
			 labels_col = '', 
			 fontsize_row = 8,
			 main = title)
	}

	compoundBob <- rownames(DILI[DILI$DILIConcern != 'N/A', ])

	pdf(paste0(inputDir, 'chapterDILI_heatmap_24H_severeAnnotated.pdf'), width = 10, height = 20)
		# plotPheatmap(na.omit(hm_24H[[1]][compoundBob, ]), title = 'GFP fraction >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_24H[[2]][compoundBob, ]), title = 'GFP fraction >2x DMSO', treatmentAnnot = 'severe')
		# plotPheatmap(na.omit(hm_24H[[3]][compoundBob, ]), title = 'GFP fraction >3x DMSO', treatmentAnnot = 'severe')
		# plotPheatmap(na.omit(hm_24H[[4]][compoundBob, ]), title = 'GFP fraction >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, 'chapterDILI_heatmap_48H_severeAnnotated.pdf'), width = 10, height = 20)
		# plotPheatmap(na.omit(hm_48H[[1]][compoundBob, ]), title = 'GFP fraction >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_48H[[2]][compoundBob, ]), title = 'GFP fraction >2x DMSO', treatmentAnnot = 'severe')
		# plotPheatmap(na.omit(hm_48H[[3]][compoundBob, ]), title = 'GFP fraction >3x DMSO', treatmentAnnot = 'severe')
		# plotPheatmap(na.omit(hm_48H[[4]][compoundBob, ]), title = 'GFP fraction >5x DMSO', treatmentAnnot = 'severe')
	dev.off()

	pdf(paste0(inputDir, 'chapterDILI_heatmap_72H_severeAnnotated.pdf'), width = 10, height = 20)
		# plotPheatmap(na.omit(hm_72H[[1]][compoundBob, ]), title = 'GFP fraction >1x DMSO', treatmentAnnot = 'severe')
		plotPheatmap(na.omit(hm_72H[[2]][compoundBob, ]), title = 'GFP fraction >2x DMSO', treatmentAnnot = 'severe')
		# plotPheatmap(na.omit(hm_72H[[3]][compoundBob, ]), title = 'GFP fraction >3x DMSO', treatmentAnnot = 'severe')
		# plotPheatmap(na.omit(hm_72H[[4]][compoundBob, ]), title = 'GFP fraction >5x DMSO', treatmentAnnot = 'severe')
	dev.off()



	## neuralnet example plots for SOT 2017 poster Baltimore
			
		pdf(paste0(inputDir, 'posterSOT_nn3_seperateNets.pdf'), width = 6, height = 6)
			nets2plot <- sapply(unique(do.call(rbind, strsplit(names(nn_3), '_'))[,1]), function(x) sample(grep(x, names(nn_3)), 1))

			for(i in nets2plot) {
				cat(i, '\n')
				tmp <- nn_3[[i]]$trainingFit

				plotnet(tmp, 
					alpha = 0.25, 
					pos_col = 'steelblue3', 
					neg_col = 'orangered3', 
					cex_val = 0.45,
					circle_cex = 2,
					max_sp = FALSE,
					x_names =  gsub('\\.', ' - ', toupper(do.call(rbind, strsplit(tmp$model.list$variables, '_'))[, 1])),
					y_names = 'DILI') 
			}
		dev.off()


		pdf(paste0(inputDir, 'posterSOT_nn3_combinedNets.pdf'), width = 12, height = 12)
			load(paste0(inputDir, 'nnOutput_50i_above3DMSO.RData'))

			nets2plot <- sample(length(nn_ALL3), 5)

			#for(i in nets2plot) {
			#	cat(i, '\n')
			#	tmp <- nn_ALL3[[i]]$trainingFit

				plotnet(tmp, 
					alpha = 0.25, 
					pos_col = 'steelblue3', 
					neg_col = 'orangered3', 
					cex_val = 0.45,
					circle_cex = 2,
					max_sp = FALSE,
					x_names =  gsub('\\.', ' - ', toupper(do.call(rbind, strsplit(tmp$model.list$variables, '_'))[, 1])),
					y_names = 'DILI') 
			# }
		dev.off()


		pdf(paste0(inputDir, 'posterSOT_HMOX_SRXN1_CHOP_25_10_10_combinedNets.pdf'), width = 12, height = 12)
		   #load(paste0(inputDir, 'nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_10_10.RData'))
			#nets2plot <- sample(length(nn_ALL3_HMOX_SRXN1_CHOP), 5)

			#for(i in nets2plot) {
			#	cat(i, '\n')
			#	tmp <- nn_ALL3[[i]]$trainingFit

				plotnet(nn_ALL3_HMOX_SRXN1_CHOP[[4]]$trainingFit, 
					alpha = 0.25, 
					pos_col = 'steelblue3', 
					neg_col = 'orangered3', 
					cex_val = 0.45,
					circle_cex = 2,
					max_sp = FALSE,
					x_names =  gsub('\\.', ' - ', toupper(do.call(rbind, strsplit(tmp$model.list$variables, '_'))[, 1])),
					y_names = 'DILI') 
			# }
		dev.off()








## PCA for MIP-DILI Jan 2017
	plotPCA <- function(pcaData, title = '', .debug = FALSE) {
		if(.debug) {
			pcaData <- hm_48H[[2]]
		}

		d <- pcaData[, 3:ncol(pcaData)]
		  colsICAM1 <-  colnames(d)[grepl('ICAM1', colnames(d)) & grepl('netchange', colnames(d))]
		  colsOthers <- colnames(d)[!grepl('ICAM1', colnames(d))]

		d <- d[, c(colsOthers, colsICAM1)]

		dDist <- dist(d)
		pca <- cmdscale(dDist, eig = TRUE, k = 5)$points
			colnames(pca) <- paste0('PC', 1:ncol(pca))

		tmp <- melt(data.frame(pcaData[, 1:2], pca), id.var = c('binaryDILI', 'doseDILI'))
		tmp <- do.call(rbind, lapply(1:5, function(i) {
  			tmp$x <- tmp$value[tmp$variable == paste0('PC', i)]
  			tmp$xPC <- paste0('PC', i)
  			tmp$yPC <- tmp$variable

  			tmp
  		} ))

  		# tmp$dili <- c('nonDILI', 'DILI')[tmp$binaryDILI + 1]
  		tmp$dili <- c('non-DILI', 'less-DILI', 'severe-DILI')[1 + tmp$doseDILI]

  		ggplot(tmp, aes(x = x, y = value)) +
	  		geom_point(aes(fill = dili), shape = 21, size = 2, alpha = 0.9) +
	  		facet_grid(yPC ~ xPC) +
	  		xlab('') + ylab('') + ggtitle(title)
	}

	pdf(paste0(inputDir, 'dose_pca_24H.pdf'), width = 12, height = 10)
		plotPCA(hm_24H[[1]], title = 'GFP positive fraction of cells >1x DMSO')
		plotPCA(hm_24H[[2]], title = 'GFP positive fraction of cells >2x DMSO')
		plotPCA(hm_24H[[3]], title = 'GFP positive fraction of cells >3x DMSO')
		plotPCA(hm_24H[[4]], title = 'GFP positive fraction of cells >5x DMSO')
	dev.off()

	pdf(paste0(inputDir, 'dose_pca_48H.pdf'), width = 12, height = 10)
		plotPCA(hm_48H[[1]], title = 'GFP positive fraction of cells >1x DMSO')
		plotPCA(hm_48H[[2]], title = 'GFP positive fraction of cells >2x DMSO')
		plotPCA(hm_48H[[3]], title = 'GFP positive fraction of cells >3x DMSO')
		plotPCA(hm_48H[[4]], title = 'GFP positive fraction of cells >5x DMSO')
	dev.off()

	pdf(paste0(inputDir, 'dose_pca_72H.pdf'), width = 12, height = 10)
		plotPCA(hm_72H[[1]], title = 'GFP positive fraction of cells >1x DMSO')
		plotPCA(hm_72H[[2]], title = 'GFP positive fraction of cells >2x DMSO')
		plotPCA(hm_72H[[3]], title = 'GFP positive fraction of cells >3x DMSO')
		plotPCA(hm_72H[[4]], title = 'GFP positive fraction of cells >5x DMSO')
	dev.off()


## LDA
	library(MASS)


	plotLDA <- function(d) {
		f <- as.formula(paste0('binaryDILI ~ ', paste(grep('CMAX', colnames(d), value = TRUE), collapse = ' + ')))
		ldaFit <- lda(f, d, CV = FALSE)


		.getCoefLDA <- function(x){
		  if (!is.null(Terms <- x$terms)) {
		    data <- model.frame(x)
		    X <- model.matrix(delete.response(Terms), data)
		    g <- model.response(data)
		    xint <- match("(Intercept)", colnames(X), nomatch = 0L)
		    if (xint > 0L) 
		      X <- X[, -xint, drop = FALSE]
		  }
		  means <- colMeans(x$means)
		  X <- scale(X, center = means, scale = FALSE) %*% x$scaling
		  rtrn <- as.data.frame(cbind(X,labels=as.character(g)))
		  rtrn <- data.frame(X,labels=as.character(g))
		  return(rtrn)
		}

		toPlot <- .getCoefLDA(ldaFit)
		toPlot$labels <- as.character(as.numeric(!(as.logical(as.numeric(toPlot$labels)))))

		ggplot(toPlot, aes(x = LD1, colour = NULL)) +
			geom_histogram(alpha = 0.4, aes(fill = labels), position = 'identity', binwidth = 0.15) +
			theme(legend.position="none")

		# ct <- table(d$binaryDILI, ldaFit$class)
		# diag(prop.table(ct, 1))
		# # total percent correct
		# sum(diag(prop.table(ct)))	
	}

	pdf(paste0(inputDir, '_lda_24H.pdf'), width = 6, height = 2)
		plotLDA(hm_24H[[1]][compoundBob, ])
		plotLDA(hm_24H[[2]][compoundBob, ])
		plotLDA(hm_24H[[3]][compoundBob, ])
		plotLDA(hm_24H[[4]][compoundBob, ])
	dev.off()

	pdf(paste0(inputDir, '_lda_48H.pdf'), width = 6, height = 2)
		plotLDA(hm_48H[[1]][compoundBob, ])
		plotLDA(hm_48H[[2]][compoundBob, ])
		plotLDA(hm_48H[[3]][compoundBob, ])
		plotLDA(hm_48H[[4]][compoundBob, ])
	dev.off()

	pdf(paste0(inputDir, '_lda_72H.pdf'), width = 6, height = 2)
		plotLDA(hm_72H[[1]][compoundBob, ])
		plotLDA(hm_72H[[2]][compoundBob, ])
		plotLDA(hm_72H[[3]][compoundBob, ])
		plotLDA(hm_72H[[4]][compoundBob, ])
	dev.off()

	# train / test splitting - DEPRECATED
	predLDA <- function(d, iterations = 50, .debug = FALSE) {
		if(.debug) {
			d <- hm_48H[[3]]
			iterations <- 50
		}

		.calcAUC <- function(truestat, testres) {
	    		# http://stats.stackexchange.com/questions/145566/how-to-calculate-area-under-the-curve-auc-or-the-c-statistic-by-hand

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
			
				rslt <- sum(height*width)	# AUC
	   	}

		set.seed(1337)
		INDEX <- list()
			for(i in 1:iterations) {
				INDEX[[i]] <- setTrainingIndex(d, 'binaryDILI', 
				negN = floor(2 * table(d$binaryDILI)/3)[1], 
				posN = floor(2 * table(d$binaryDILI)/3)[2])
			}

		pr <- mclapply(1:iterations, function(i) {
			index <- INDEX[[i]]

			train <- d[index,  ]
  	  		test  <- d[-index, ]

  	  		f <- as.formula(paste0('binaryDILI ~ ', paste(grep('CMAX', colnames(train), value = TRUE), collapse = ' + ')))
			ldaFit <- lda(f, train, CV = FALSE)

			predict(ldaFit, test)$class
		})
		
		AUCs <- sapply(pr, function(predDILI) .calcAUC(truestat = round(test$binaryDILI), testres = predDILI))

		# ...
	}


## neuralnetting - server continuation
	options(stringsAsFactors = FALSE)
	set.seed(1337)
	
	library(neuralnet)
	library(nnetpredint)
	library(NeuralNetTools)
	library(caret)

	library(ggplot2)

	library(data.table)
	library(preprocessCore)
	library(reshape2)
	library(parallel)

	load('/data/wouter/DILI/nnInput.RData')

  ## Drop control compounds
	compoundBob <- rownames(DILI[DILI$DILIConcern != 'N/A', ])

  ## Set indici
	setTrainingIndex <- function(data, class = c('binaryDILI', 'doseDILI'), negN = 21, posN = 70) {
		if(length(class) != 1) stop('Only supply a single column name as dependent variable.\n')
  		
  		neg <- which(data[, class] == 0)
  		pos <- which(data[, class] == 1)

  		trainingNegative <- sample(neg, negN)
  		trainingPositive <- sample(pos, posN)

  		c(trainingNegative, trainingPositive)
  	}

  	set.seed(1337)
  	INDEX <- list()
	for(i in 1:1000) INDEX[[i]] <- setTrainingIndex(classifierInput_24, 'binaryDILI')


  ## Examplary plot
  	# nn_48H    <- lapply(1:5, function(i)  fitNN(dat = classifierInput_48, NodesLayers = c(24, 12, 12), inputParameterGrep = '_2', index = INDEX[[i]], nCores = 10))
  	
  	nn_BIP    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BIP[compoundBob, ],    NodesLayers = c(10, 10, 3), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_SRXN1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_SRXN1[compoundBob, ],  NodesLayers = c(10, 10, 3), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_P21    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_P21[compoundBob, ],    NodesLayers = c(10, 10, 3), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_CHOP   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_CHOP[compoundBob, ],   NodesLayers = c(10, 10, 3), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_HSPA1B <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HSPA1B[compoundBob, ], NodesLayers = c(10, 10, 3), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_ICAM1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ICAM1[compoundBob, ],  NodesLayers = c(10, 10, 3), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_BTG2   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BTG2[compoundBob, ],   NodesLayers = c(10, 10, 3), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_HMOX1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1[compoundBob, ],  NodesLayers = c(10, 10, 3), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
		nn_3 <- c(nn_BIP, nn_SRXN1, nn_P21, nn_CHOP, nn_HSPA1B, nn_ICAM1, nn_BTG2, nn_HMOX1)
		names(nn_3) <- c(paste0('BIP_', 1:50), paste0('SRXN1_', 1:50), paste0('P21_', 1:50), paste0('CHOP_', 1:50), paste0('HSPA1B_', 1:50), paste0('ICAM1_', 1:50), paste0('BTG2_', 1:50), paste0('HMOX1_', 1:50))
		save(nn_3, file = '/data/wouter/DILI/nnOutput_50i_above3DMSO__3layers-10-10-3.RData')


	nn_BIP    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BIP[compoundBob, ],    NodesLayers = c(6, 6), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_SRXN1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_SRXN1[compoundBob, ],  NodesLayers = c(6, 6), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_P21    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_P21[compoundBob, ],    NodesLayers = c(6, 6), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_CHOP   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_CHOP[compoundBob, ],   NodesLayers = c(6, 6), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_HSPA1B <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HSPA1B[compoundBob, ], NodesLayers = c(6, 6), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_ICAM1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ICAM1[compoundBob, ],  NodesLayers = c(6, 6), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_BTG2   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BTG2[compoundBob, ],   NodesLayers = c(6, 6), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_HMOX1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1[compoundBob, ],  NodesLayers = c(6, 6), inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
		nn_3 <- c(nn_BIP, nn_SRXN1, nn_P21, nn_CHOP, nn_HSPA1B, nn_ICAM1, nn_BTG2, nn_HMOX1)
		names(nn_3) <- c(paste0('BIP_', 1:50), paste0('SRXN1_', 1:50), paste0('P21_', 1:50), paste0('CHOP_', 1:50), paste0('HSPA1B_', 1:50), paste0('ICAM1_', 1:50), paste0('BTG2_', 1:50), paste0('HMOX1_', 1:50))
		save(nn_3, file = '/data/wouter/DILI/nnOutput_50i_above3DMSO__3layers-6-6.RData')		


  ## TIL
	nn_BIP    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BIP[compoundBob, ],    inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
	nn_SRXN1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_SRXN1[compoundBob, ],  inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
	nn_P21    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_P21[compoundBob, ],    inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
	nn_CHOP   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_CHOP[compoundBob, ],   inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
	nn_HSPA1B <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HSPA1B[compoundBob, ], inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
	nn_ICAM1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ICAM1[compoundBob, ],  inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
	nn_BTG2   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BTG2[compoundBob, ],   inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
	nn_HMOX1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1[compoundBob, ],  inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
		nn_1 <- c(nn_BIP, nn_SRXN1, nn_P21, nn_CHOP, nn_HSPA1B, nn_ICAM1, nn_BTG2, nn_HMOX1)
		names(nn_1) <- c(paste0('BIP_', 1:50), paste0('SRXN1_', 1:50), paste0('P21_', 1:50), paste0('CHOP_', 1:50), paste0('HSPA1B_', 1:50), paste0('ICAM1_', 1:50), paste0('BTG2_', 1:50), paste0('HMOX1_', 1:50))
		save(nn_1, file = '/data/wouter/DILI/nnOutput_50i_above1DMSO.RData')

	nn_BIP    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BIP[compoundBob, ],    inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_SRXN1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_SRXN1[compoundBob, ],  inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_P21    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_P21[compoundBob, ],    inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_CHOP   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_CHOP[compoundBob, ],   inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_HSPA1B <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HSPA1B[compoundBob, ], inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_ICAM1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ICAM1[compoundBob, ],  inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_BTG2   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BTG2[compoundBob, ],   inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_HMOX1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1[compoundBob, ],  inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
		nn_2 <- c(nn_BIP, nn_SRXN1, nn_P21, nn_CHOP, nn_HSPA1B, nn_ICAM1, nn_BTG2, nn_HMOX1)
		names(nn_2) <- c(paste0('BIP_', 1:50), paste0('SRXN1_', 1:50), paste0('P21_', 1:50), paste0('CHOP_', 1:50), paste0('HSPA1B_', 1:50), paste0('ICAM1_', 1:50), paste0('BTG2_', 1:50), paste0('HMOX1_', 1:50))
		save(nn_2, file = '/data/wouter/DILI/nnOutput_50i_above2DMSO.RData')

	nn_24H    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_24[compoundBob, ],    inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_48H    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_48[compoundBob, ],    inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
	nn_72H    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_72[compoundBob, ],    inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
		names(nn_24H) <- paste0('nn_24H_', 1:50)
		names(nn_48H) <- paste0('nn_48H_', 1:50)
		names(nn_72H) <- paste0('nn_72H_', 1:50)

		save(nn_24H, file = '/data/wouter/DILI/nnOutput_50i_24H_above2DMSO.RData')
		save(nn_48H, file = '/data/wouter/DILI/nnOutput_50i_48H_above2DMSO.RData')
		save(nn_72H, file = '/data/wouter/DILI/nnOutput_50i_72H_above2DMSO.RData')

	treatmentsInAll     <- names(which(table(c(rownames(classifierInput_24), rownames(classifierInput_48), rownames(classifierInput_72))) == 3))
	classifierInput_ALL <- data.frame(classifierInput_24[treatmentsInAll, ], 
									  classifierInput_48[treatmentsInAll, 3:ncol(classifierInput_48)], 
									  classifierInput_72[treatmentsInAll, 3:ncol(classifierInput_72)])

	nn_ALL1    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ALL[compoundBob, ],    inputParameterGrep = '_1', splitBy = 'binaryDILI', nCores = 10))
		names(nn_ALL1) <- paste0('nn_above1DMSO_', 1:50)
		save(nn_ALL1, file = '/data/wouter/DILI/nnOutput_50i_above1DMSO_reportersCombined.RData')

	nn_ALL2    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ALL[compoundBob, ],    inputParameterGrep = '_2', splitBy = 'binaryDILI', nCores = 10))
		names(nn_ALL2) <- paste0('nn_above2DMSO_', 1:50)
		save(nn_ALL2, file = '/data/wouter/DILI/nnOutput_50i_above2DMSO_reportersCombined.RData')

	nn_ALL3    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ALL[compoundBob, ],    inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
		names(nn_ALL3) <- paste0('nn_above3DMSO_', 1:50)
		save(nn_ALL3, file = '/data/wouter/DILI/nnOutput_50i_above3DMSO_reportersCombined.RData')

	nn_ALL5    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ALL[compoundBob, ],    inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
		names(nn_ALL5) <- paste0('nn_above5DMSO_', 1:50)
		save(nn_ALL5, file = '/data/wouter/DILI/nnOutput_50i_above5DMSO_reportersCombined.RData')

	nn_ALL_CMAX <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ALL[compoundBob, ],    inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
		names(nn_ALL_CMAX) <- paste0('nn_allCMAXPars_', 1:50)
		save(nn_ALL_CMAX, file = '/data/wouter/DILI/nnOutput_50i_allCMAXPars_reportersCombined.RData')


	## TOX PPT Mon 20 Feb 2017
	ind <- c(grep('HMOX', colnames(classifierInput_ALL)), grep('SRXN', colnames(classifierInput_ALL)), grep('CHOP', colnames(classifierInput_ALL)))
	classifierInput_HMOX1_SRXN1_CHOP <- classifierInput_ALL[, c(1, 2, ind)]

	nn_ALL3_HMOX_SRXN1_CHOP <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1_SRXN1_CHOP[compoundBob, ],    inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 20, NodesLayers = c(50, 10)))
		names(nn_ALL3_HMOX_SRXN1_CHOP) <- paste0('nn_above3DMSO_', 1:50)
		save(nn_ALL3_HMOX_SRXN1_CHOP, file = '/data/wouter/DILI/nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined.RData')

	nn_ALL3_HMOX_SRXN1_CHOP <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1_SRXN1_CHOP[compoundBob, ],    inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 20, NodesLayers = c(25, 10)))
		names(nn_ALL3_HMOX_SRXN1_CHOP) <- paste0('nn_above3DMSO_', 1:50)
		save(nn_ALL3_HMOX_SRXN1_CHOP, file = '/data/wouter/DILI/nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_10.RData')	

	nn_ALL3_HMOX_SRXN1_CHOP <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1_SRXN1_CHOP[compoundBob, ],    inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 20, NodesLayers = c(25, 10, 10)))
		names(nn_ALL3_HMOX_SRXN1_CHOP) <- paste0('nn_above3DMSO_', 1:50)
		save(nn_ALL3_HMOX_SRXN1_CHOP, file = '/data/wouter/DILI/nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_10_10.RData')	

	nn_ALL3_HMOX_SRXN1_CHOP <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1_SRXN1_CHOP[compoundBob, ],    inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 20, NodesLayers = c(25, 50, 50)))
		names(nn_ALL3_HMOX_SRXN1_CHOP) <- paste0('nn_above3DMSO_', 1:50)
		save(nn_ALL3_HMOX_SRXN1_CHOP, file = '/data/wouter/DILI/nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_50_50.RData')	





  ## EMT
  	nn_BIP    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BIP[compoundBob, ],    inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_SRXN1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_SRXN1[compoundBob, ],  inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_P21    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_P21[compoundBob, ],    inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_CHOP   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_CHOP[compoundBob, ],   inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_HSPA1B <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HSPA1B[compoundBob, ], inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_ICAM1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ICAM1[compoundBob, ],  inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_BTG2   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BTG2[compoundBob, ],   inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
	nn_HMOX1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1[compoundBob, ],  inputParameterGrep = '_3', splitBy = 'binaryDILI', nCores = 10))
		nn_3 <- c(nn_BIP, nn_SRXN1, nn_P21, nn_CHOP, nn_HSPA1B, nn_ICAM1, nn_BTG2, nn_HMOX1)
		names(nn_3) <- c(paste0('BIP_', 1:50), paste0('SRXN1_', 1:50), paste0('P21_', 1:50), paste0('CHOP_', 1:50), paste0('HSPA1B_', 1:50), paste0('ICAM1_', 1:50), paste0('BTG2_', 1:50), paste0('HMOX1_', 1:50))
		save(nn_3, file = '/data/wouter/DILI/nnOutput_50i_above3DMSO.RData')
	

	nn_BIP    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BIP[compoundBob, ],    inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
	nn_SRXN1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_SRXN1[compoundBob, ],  inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
	nn_P21    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_P21[compoundBob, ],    inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
	nn_CHOP   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_CHOP[compoundBob, ],   inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
	nn_HSPA1B <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HSPA1B[compoundBob, ], inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
	nn_ICAM1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ICAM1[compoundBob, ],  inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
	nn_BTG2   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BTG2[compoundBob, ],   inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
	nn_HMOX1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1[compoundBob, ],  inputParameterGrep = '_5', splitBy = 'binaryDILI', nCores = 10))
		nn_5 <- c(nn_BIP, nn_SRXN1, nn_P21, nn_CHOP, nn_HSPA1B, nn_ICAM1, nn_BTG2, nn_HMOX1)
		names(nn_5) <- c(paste0('BIP_', 1:50), paste0('SRXN1_', 1:50), paste0('P21_', 1:50), paste0('CHOP_', 1:50), paste0('HSPA1B_', 1:50), paste0('ICAM1_', 1:50), paste0('BTG2_', 1:50), paste0('HMOX1_', 1:50))
		
		save(nn_5, file = '/data/wouter/DILI/nnOutput_50i_above5DMSO.RData')

  ## TOX	
	nn_BIP    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BIP[compoundBob, ],    inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
	nn_SRXN1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_SRXN1[compoundBob, ],  inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
	nn_P21    <- lapply(1:50, function(i)  fitNN(dat = classifierInput_P21[compoundBob, ],    inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
	nn_CHOP   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_CHOP[compoundBob, ],   inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
	nn_HSPA1B <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HSPA1B[compoundBob, ], inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
	nn_ICAM1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_ICAM1[compoundBob, ],  inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
	nn_BTG2   <- lapply(1:50, function(i)  fitNN(dat = classifierInput_BTG2[compoundBob, ],   inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
	nn_HMOX1  <- lapply(1:50, function(i)  fitNN(dat = classifierInput_HMOX1[compoundBob, ],  inputParameterGrep = 'CMAX', splitBy = 'binaryDILI', nCores = 10))
		nn_CMAX <- c(nn_BIP, nn_SRXN1, nn_P21, nn_CHOP, nn_HSPA1B, nn_ICAM1, nn_BTG2, nn_HMOX1)
		names(nn_CMAX) <- c(paste0('BIP_', 1:50), paste0('SRXN1_', 1:50), paste0('P21_', 1:50), paste0('CHOP_', 1:50), paste0('HSPA1B_', 1:50), paste0('ICAM1_', 1:50), paste0('BTG2_', 1:50), paste0('HMOX1_', 1:50))
		save(nn_CMAX, file = '/data/wouter/DILI/nnOutput_50i_allParameters.RData')
	

## Neuralnetting - local plotting of AUCs
	plotOutcomeNN <- function(nnList, feature = c('AUCs', 'spec', 'sens'), .debug = FALSE, offset = 0, returnData = FALSE) {
		if(.debug) {
			nnList <- nn_3
		}

		nnList <- nnList[which(sapply(nnList, function(nnFit) length(nnFit[[feature]] )) == 10)]
		df <- data.frame(nnFit = names(nnList), 
					 reporter = do.call(rbind, strsplit(names(nnList), '_'))[,1], 
					 indexID  = do.call(rbind, strsplit(names(nnList), '_'))[,2])
			rownames(df) <- df$nnFit
		
		dfFeat <- data.frame(t(sapply(nnList, function(nnFit) nnFit[[feature]])))
			colnames(dfFeat) <- paste0('auc_N', 1:10)

		if(identical(rownames(df), rownames(dfFeat))) {
			dfFeat <- data.frame(df, dfFeat)
		} else {
			stop('Something went wrong!\n')
		}	

		mdfFeat <- melt(dfFeat, id.var = c('nnFit', 'reporter', 'indexID'))
		mdfFeat$value <- mdfFeat$value + offset

		if(returnData) return(list(dfFeat = dfFeat, mdfFeat = mdfFeat))

		ggplot(mdfFeat, aes(x = reporter, y = value)) +
			geom_hline(yintercept = 0.5) +
			geom_boxplot() +
			ylim(0, 1) + ggtitle(feature) +
			theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
			ylab('') + xlab('')

	}

  # MIP-DILI poster AUC plots
  	pdf(paste0(inputDir, 'POSTER_nn_AUC_perReporter_allTimepoints_above3DMSO.pdf'), width = 3, height = 5)
		# load(paste0(inputDir, 'nnOutput_50i_above3DMSO.RData'))			# nn_3
		plotOutcomeNN(nn_3, 'AUCs', offset = 0)
		# rm(nn_2)
	dev.off()

	pdf(paste0(inputDir, 'POSTER_nn_AUC_combinedReporter_allTimepoints_above2DMSO.pdf'), width = 1.5, height = 5)
		# load(paste0(inputDir, 'nnOutput_50i_72H_above2DMSO.RData'))			# nn_72H
		plotOutcomeNN(nn_72H, 'AUCs', offset = 0)
		# rm(nn_72H)
	dev.off()

  # TOX PPT Mon 20 Feb 2017
  	pdf(paste0('/Users/Wouter/Desktop/', 'TOXPPT_nn_AUC_perReporter_allTimepoints_above3DMSO.pdf'), width = 3, height = 5)
		# load(paste0(inputDir, 'nnOutput_50i_above3DMSO.RData'))			# nn_3
		plotOutcomeNN(nn_3, 'AUCs', offset = 0)
		# rm(nn_2)
	dev.off()

	pdf(paste0('/Users/Wouter/Desktop/', 'TOXPPT_nn_sens_perReporter_allTimepoints_above3DMSO.pdf'), width = 3, height = 5)
		# load(paste0(inputDir, 'nnOutput_50i_above3DMSO.RData'))			# nn_3
		plotOutcomeNN(nn_3, 'spec', offset = 0)
		# rm(nn_2)
	dev.off()

	pdf(paste0('/Users/Wouter/Desktop/', 'TOXPPT_nn_spec_perReporter_allTimepoints_above3DMSO.pdf'), width = 3, height = 5)
		# load(paste0(inputDir, 'nnOutput_50i_above3DMSO.RData'))			# nn_3
		plotOutcomeNN(nn_3, 'sens', offset = 0)
		# rm(nn_2)
	dev.off()


	# HMOX1, SRXN1 & CHOP
	load(paste0(inputDir, 'nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_10.RData'))
	pdf(paste0(inputDir, 'TOXPPT_nn_Accuracy_HMOX1_SRXN1_CHOP_allTimepoints_above3DMSO___25_10.pdf'), width = 1.5, height = 5)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'AUCs', offset = 0)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'spec', offset = 0)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'sens', offset = 0)
	dev.off()


	load(paste0(inputDir, 'nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_10_10.RData'))
	pdf('/Users/Wouter/Desktop/HMOX1_SRXN1_CHOP_exampleNN_25_10_10.pdf', width = 8, height = 5)
		tmp <- nn_ALL3_HMOX_SRXN1_CHOP[[2]]$trainingFit
		plotnet(tmp, 
						alpha = 0.25, 
						pos_col = 'steelblue3', 
						neg_col = 'orangered3', 
						cex_val = 0.3,
						circle_cex = 2,
						max_sp = FALSE,
						x_names =  gsub('\\.', ' - ', toupper(do.call(rbind, strsplit(tmp$model.list$variables, '_'))[, 1])),
						y_names = 'DILI') 
	dev.off()

	load(paste0(inputDir, 'nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_10_10.RData'))
	pdf(paste0(inputDir, 'TOXPPT_nn_Accuracy_HMOX1_SRXN1_CHOP_allTimepoints_above3DMSO___25_10_10.pdf'), width = 1.5, height = 5)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'AUCs', offset = 0)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'spec', offset = 0)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'sens', offset = 0)
	dev.off()

	load(paste0(inputDir, 'nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_50_50.RData'))
	pdf(paste0(inputDir, 'TOXPPT_nn_Accuracy_HMOX1_SRXN1_CHOP_allTimepoints_above3DMSO___25_50_50.pdf'), width = 1.5, height = 5)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'AUCs', offset = 0)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'spec', offset = 0)
		plotOutcomeNN(nn_ALL3_HMOX_SRXN1_CHOP, 'sens', offset = 0)
	dev.off()

	nn_ALL3_HMOX_SRXN1_CHOP_compPredScores <-  mclapply(nn_ALL3_HMOX_SRXN1_CHOP, getCompPred, mc.cores = 6)

	pdf(paste0(inputDir, 'TOXPPT_nn_compPredScores_HMOX1_SRXN1_CHOP_allTimepoints_above3DMSO.pdf'), width = 12, height = 4)
		plotThoseScores(nn_ALL3_HMOX_SRXN1_CHOP_compPredScores)
	dev.off()	

	pdf('/Users/Wouter/Desktop/TOXPPT_nn_HMOX1_SRXN1_CHOP_allTimepoints_above3DMSO.pdf', width = 12, height = 6)
		tmp <- nn_ALL3_HMOX_SRXN1_CHOP[[10]]$trainingFit
			plotnet(tmp, 
							alpha = 0.25, 
							pos_col = 'steelblue3', 
							neg_col = 'orangered3', 
							cex_val = 0.3,
							circle_cex = 2,
							max_sp = FALSE,
							x_names =  gsub('\\.', ' - ', toupper(do.call(rbind, strsplit(tmp$model.list$variables, '_'))[, 1])),
							y_names = 'DILI') 
	dev.off()	

	highSpec <- data.frame(which(sapply(nn_ALL3_HMOX_SRXN1_CHOP, function(x) x$sens ) > 0.7, arr.ind = TRUE))
	highSpec[, 'iter'] <- names(nn_ALL3_HMOX_SRXN1_CHOP)[highSpec[, 'col']]
	highSpec[, 'uniqueID'] <- paste0(highSpec$iter, '_', highSpec$row)

	highSens <- data.frame(which(sapply(nn_ALL3_HMOX_SRXN1_CHOP, function(x) x$spec ) > 0.7, arr.ind = TRUE))
	highSens[, 'iter'] <- names(nn_ALL3_HMOX_SRXN1_CHOP)[highSens[, 'col']]
	highSens[, 'uniqueID'] <- paste0(highSens$iter, '_', highSens$row)

	highAUC <- data.frame(which(sapply(nn_ALL3_HMOX_SRXN1_CHOP, function(x) x$AUCs ) > 0.7, arr.ind = TRUE))
	highAUC[, 'iter'] <- names(nn_ALL3_HMOX_SRXN1_CHOP)[highAUC[, 'col']]
	highAUC[, 'uniqueID'] <- paste0(highAUC$iter, '_', highAUC$row)

	# goodPerformingNNs <- highSpec[which(highSpec$uniqueID %in% highSens$uniqueID), ]
	goodPerformingNNs <- highAUC

	goodNN_ALL3_HMOX_SRXN1_CHOP <- nn_ALL3_HMOX_SRXN1_CHOP[unique(goodPerformingNNs$iter)]
	for(i in names(goodNN_ALL3_HMOX_SRXN1_CHOP)) {
		.ind <- goodPerformingNNs$row[which(goodPerformingNNs$iter == i)]
		goodNN_ALL3_HMOX_SRXN1_CHOP[[i]]$testPred <- goodNN_ALL3_HMOX_SRXN1_CHOP[[i]]$testPred[.ind]
	}

	tmp <- lapply(goodNN_ALL3_HMOX_SRXN1_CHOP, getCompPred)


	# Acquire compounds that are often predicted wrong / good
	getCompPred <- function(nnOutput, .debug = FALSE) {
		if(.debug) {
			nnOutput <- nn_3[[1]]
		}

		trueTestClass <- nnOutput$test$binaryDILI
		predTestClass <- round(do.call(cbind, nnOutput$testPred))

		rowSums(apply(predTestClass, 2, function(x) x == trueTestClass)) / 10
	}

	plotThoseScores <- function(compPredScores) {
		compPredScores <- lapply(compPredScores, function(x) x[!is.na(x)] )
		compounds <- unique(unlist(lapply(compPredScores, names)))# [1:131]

		compPredScores <- lapply(compounds, function(x) as.numeric(na.omit(sapply(compPredScores, function(y) y[x]  ))))
			names(compPredScores) <- compounds

		compoundAnnotation <- DILI[compounds, 'binaryDILI', drop = FALSE]	

		for(i in names(compPredScores)) {
			compPredScores[[i]] <- data.frame(compound = i, predScore = compPredScores[[i]])
		}

		compPredScores <- compPredScores[names(sort(sapply(compPredScores, function(x) mean(x$predScore))))]
		compPredScoresLong <- do.call(rbind, compPredScores)

		# ggplot(compPredScoresLong, aes(x = compound, y = predScore)) + geom_boxplot()


		meanPredScores <- data.frame(compound = names(compPredScores),
									 meanPredScore = sapply(compPredScores, function(x) mean(x$predScore)),
									 binaryDILI = c('nonDILI', 'DILI')[1 + DILI[names(compPredScores), 'binaryDILI']])
		

		ggplot(meanPredScores, aes(x = compound, y = meanPredScore, fill = binaryDILI)) + 
			 geom_bar(position = position_dodge(), stat = 'identity') +
			 scale_x_discrete(limits = meanPredScores$compound) +
			 theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1), legend.position = 'none') +
			 ylim(0, 1) + xlab('') + ylab('')
	}

	plotThoseScores_SOT <- function(compPredScores, y.textSize = 12, x.textSize = 10, returnData = FALSE) {
		compPredScores <- lapply(compPredScores, function(x) x[!is.na(x)] )
		compounds <- unique(unlist(lapply(compPredScores, names)))# [1:131]

		compPredScores <- lapply(compounds, function(x) as.numeric(na.omit(sapply(compPredScores, function(y) y[x]  ))))
			names(compPredScores) <- compounds

		compoundAnnotation <- DILI[compounds, 'doseDILI', drop = FALSE]	

		for(i in names(compPredScores)) {
			compPredScores[[i]] <- data.frame(compound = i, predScore = compPredScores[[i]])
		}

		compPredScores <- compPredScores[names(sort(sapply(compPredScores, function(x) mean(x$predScore))))]
		compPredScoresLong <- do.call(rbind, compPredScores)

		# ggplot(compPredScoresLong, aes(x = compound, y = predScore)) + geom_boxplot()


		meanPredScores <- data.frame(compound = names(compPredScores),
									 meanPredScore = sapply(compPredScores, function(x) mean(x$predScore)),
									 DILI = c('non-DILI', 'less-DILI', 'severe-DILI')[1 + DILI[names(compPredScores), 'doseDILI']])
		
		if(returnData) return(list(compPredScores = compPredScores, meanPredScores = meanPredScores))

		ggplot(meanPredScores, aes(x = compound, y = meanPredScore, fill = DILI)) + 
			 geom_bar(position = position_dodge(), stat = 'identity') +
			 # scale_x_discrete(limits = meanPredScores$compound) +
			 theme(axis.text.y = element_text(size = y.textSize), 
			 	   axis.text.x = element_text(size = x.textSize, angle = 45, hjust = 1), 
			 	   strip.text.x = element_text(size = 20),
			 	   legend.position = 'none') +
			 ylim(0, 1) + xlab('') + ylab('') +
			 facet_grid(. ~ DILI, scales = 'free_x', space = 'free_x')
	}

	plotThoseScoresPerReporter <- function(compPredScores, y.textSize = 12, x.textSize = 10, returnData = FALSE) {
		compPredScores <- lapply(compPredScores, function(x) x[!is.na(x)] )
		compounds <- unique(unlist(lapply(compPredScores, names)))# [1:131]

		#compPredScores <- lapply(compounds, function(x) as.numeric(na.omit(sapply(compPredScores, function(y) y[x]  ))))
		#		names(compPredScores) <- compounds

		compPredScores <- lapply(compounds, function(x) {
			tmp <- sapply(compPredScores, function(y) y[x])
			tmp[!is.na(tmp)]
		}  )
			names(compPredScores) <- compounds

		compoundAnnotation <- DILI[compounds, 'doseDILI', drop = FALSE]	

		for(i in names(compPredScores)) {
			tmp <- data.frame(compound = i, predScore = compPredScores[[i]])
			tmp$reporter <- do.call(rbind, strsplit(rownames(tmp), '_'))[, 1]
			compPredScores[[i]] <- tmp
		}

		compPredScores <- compPredScores[names(sort(sapply(compPredScores, function(x) mean(x$predScore))))]
		reporters <- sort(unique(unlist(sapply(compPredScores, function(x) x$reporter))))

		meanPredScores <- do.call(rbind, lapply(reporters, function(.rep) {
			do.call(rbind, lapply(compPredScores, function(y) y[which(y$reporter == .rep), ]))
			}))
		# meanPredScores$doseDILI <- c('non-DILI', 'less-DILI', 'severe-DILI')[1 + DILI[names(compPredScores), 'doseDILI']]
		meanPredScores$doseDILI <- c('non-DILI', 'less-DILI', 'severe-DILI')[1 + DILI[meanPredScores$compound, 'doseDILI']]
		
		meanPredScores$UID <- paste0(meanPredScores$compound, '_', meanPredScores$reporter)

		.meanPredScores <- data.frame(UID = unique(meanPredScores$UID), predScore = NA)
			.meanPredScores$compound <- do.call(rbind, strsplit(.meanPredScores$UID, '_'))[, 1]
			.meanPredScores$reporter <- do.call(rbind, strsplit(.meanPredScores$UID, '_'))[, 2]
			.meanPredScores$DILI <- c('non-DILI', 'less-DILI', 'severe-DILI')[1 + DILI[.meanPredScores$compound, 'doseDILI']]

			for(uid in .meanPredScores$UID) .meanPredScores$predScore[.meanPredScores$UID == uid] <- mean(meanPredScores$predScore[meanPredScores$UID == uid])

		if(returnData) return(list(compPredScores = compPredScores, meanPredScores = meanPredScores))		

		ggplot(.meanPredScores, aes(x = compound, y = predScore, fill = DILI)) + 
			 geom_bar(position = position_dodge(), stat = 'identity') +
			 # scale_x_discrete(limits = unique(meanPredScores$compound)) +
			 theme(axis.text.y = element_text(size = y.textSize), 
			 	   axis.text.x = element_text(size = x.textSize, angle = 45, hjust = 1), 
			 	   strip.text.x = element_text(size = 20),
			 	   legend.position = 'none') +
			 	   ylim(0, 1) + xlab('') + ylab('') +
			 facet_grid(reporter ~ DILI, scales = 'free_x', space = 'free_x')
	}


	# nn_1_compPredScores    <- mclapply(nn_1, getCompPred, mc.cores = 6)
	# nn_2_compPredScores    <- mclapply(nn_2, getCompPred, mc.cores = 6)
	# nn_3_compPredScores    <- mclapply(nn_3, getCompPred, mc.cores = 6)
	# nn_5_compPredScores    <- mclapply(nn_5, getCompPred, mc.cores = 6)
	# nn_CMAX_compPredScores <- mclapply(nn_CMAX, getCompPred, mc.cores = 6)

	# SOT POSTER 2017 NEURALNET NN ANN - all timepoints
	  # Seperate reporters - correct prediction rate
		pdf('/Users/Wouter/Desktop/combined_nn_3_compoundPredictabilityScores.pdf', width = 36, height = 6)
			plotThoseScores_SOT(nn_3_compPredScores)
		dev.off()	

		pdf('/Users/Wouter/Desktop/seperate_nn_3_compoundPredictabilityScores.pdf', width = 36, height = 24)
			plotThoseScoresPerReporter(nn_3_compPredScores)
		dev.off()	

		pdf('/Users/Wouter/Desktop/seperate_nn_3_compoundPredictabilityScores_AUC.pdf', width = 3, height = 5)
			plotOutcomeNN(nn_3, 'AUCs', offset = 0)
		dev.off()	

		pdf('/Users/Wouter/Desktop/combined_nn_3_compoundPredictabilityScores_AUC.pdf', width = 3, height = 5)
			plotOutcomeNN(nn_ALL3, 'AUCs', offset = 0)
		dev.off()	

	  # HMOX1, SRXN1 & CHOP combined - correct prediction rate
		load(paste0(inputDir, 'nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_10_10.RData'))
		nn_ALL3_HMOX_SRXN1_CHOP_compPredScores <- mclapply(nn_ALL3_HMOX_SRXN1_CHOP, getCompPred, mc.cores = 6)

		pdf('/Users/Wouter/Desktop/combined_nn_HMOX_SRXN1_CHOP_compoundPredictabilityScores_25_10_10.pdf', width = 12, height = 3)
			plotThoseScores_SOT(nn_ALL3_HMOX_SRXN1_CHOP_compPredScores)
		dev.off()	

		load(paste0(inputDir, 'nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_10.RData'))
		nn_ALL3_HMOX_SRXN1_CHOP_compPredScores <- mclapply(nn_ALL3_HMOX_SRXN1_CHOP, getCompPred, mc.cores = 6)

		pdf('/Users/Wouter/Desktop/combined_nn_HMOX_SRXN1_CHOP_compoundPredictabilityScores_25_10.pdf', width = 12, height = 3)
			plotThoseScores_SOT(nn_ALL3_HMOX_SRXN1_CHOP_compPredScores)
		dev.off()	

		load(paste0(inputDir, 'nnOutput_50i_above3DMSO_HMOX1_SRXN1_CHOP_Combined__layers_25_50_50.RData'))
		nn_ALL3_HMOX_SRXN1_CHOP_compPredScores <- mclapply(nn_ALL3_HMOX_SRXN1_CHOP, getCompPred, mc.cores = 6)

		pdf('/Users/Wouter/Desktop/combined_nn_HMOX_SRXN1_CHOP_compoundPredictabilityScores_25_50_50.pdf', width = 12, height = 3)
			plotThoseScores_SOT(nn_ALL3_HMOX_SRXN1_CHOP_compPredScores)
		dev.off()	

	  # all reporters combined - correct prediction rate
	  	load(paste0(inputDir, 'nnOutput_50i_above3DMSO_reportersCombined.RData'))			# nn_ALL3
	  	nn_ALL3_compPredScores <- mclapply(nn_ALL3, getCompPred, mc.cores = 6)

		pdf('/Users/Wouter/Desktop/combined_nn_ALL3_compoundPredictabilityScores.pdf', width = 12, height = 3)
			plotThoseScores_SOT(nn_ALL3_compPredScores)
		dev.off()	

		

	# plotThoseScores(nn_1_compPredScores)
	# plotThoseScores(nn_2_compPredScores)
	pdf('/Users/Wouter/Desktop/nn_3_compoundPredictabilityScores.pdf', width = 12, height = 4)
		plotThoseScores(nn_3_compPredScores)
	dev.off()	

	pdf('/Users/Wouter/Desktop/nn_3_compoundPredictabilityScores_HMOX1.pdf', width = 12, height = 4)
	 	nn_3_hmox1_compPredScores <- nn_3_compPredScores[grep('HMOX1', names(nn_3_compPredScores))]
	 	plotThoseScores(nn_3_hmox1_compPredScores)
	dev.off() 

	pdf('/Users/Wouter/Desktop/nn_3_compoundPredictabilityScores_SRXN1.pdf', width = 12, height = 4)
	 	nn_3_SRXN1_compPredScores <- nn_3_compPredScores[grep('SRXN1', names(nn_3_compPredScores))]
	 	plotThoseScores(nn_3_SRXN1_compPredScores)
	dev.off() 

	pdf('/Users/Wouter/Desktop/nn_3_compoundPredictabilityScores_CHOP.pdf', width = 12, height = 4)
	 	nn_3_CHOP_compPredScores <- nn_3_compPredScores[grep('CHOP', names(nn_3_compPredScores))]
	 	plotThoseScores(nn_3_CHOP_compPredScores)
	dev.off() 

	pdf('/Users/Wouter/Desktop/nn_3_compoundPredictabilityScores_BIP.pdf', width = 12, height = 4)
	 	nn_3_BIP_compPredScores <- nn_3_compPredScores[grep('BIP', names(nn_3_compPredScores))]
	 	plotThoseScores(nn_3_BIP_compPredScores)
	dev.off() 

	# plotThoseScores(nn_5_compPredScores)
	# plotThoseScores(nn_CMAX_compPredScores)

	load(paste0(inputDir, 'nnOutput_50i_above3DMSO__3layers-6-6.RData'))		## nn_3
	nn_3__6_6 <- nn_3
	rm(nn_3)

	load(paste0(inputDir, 'nnOutput_50i_above3DMSO__3layers-10-10-3.RData'))	## nn_3
	nn_3__10_10_3 <- nn_3
	rm(nn_3)

	nn_3__6_6_compPredScores <- mclapply(nn_3__6_6, getCompPred, mc.cores = 6)
	nn_3__10_10_3_compPredScores <- mclapply(nn_3__10_10_3, getCompPred, mc.cores = 6)

	plotThoseScores(nn_3__6_6_compPredScores)
	plotThoseScores(nn_3__10_10_3_compPredScores)

	
	

	pdf('/Users/Wouter/Desktop/nn_3__10_10_3__CHOP_exampleNN.pdf', width = 8, height = 5)
		tmp <- nn_3__10_10_3[[182]]$trainingFit
		plotnet(tmp, 
						alpha = 0.25, 
						pos_col = 'steelblue3', 
						neg_col = 'orangered3', 
						cex_val = 0.3,
						circle_cex = 2,
						max_sp = FALSE,
						x_names =  gsub('\\.', ' - ', toupper(do.call(rbind, strsplit(tmp$model.list$variables, '_'))[, 1])),
						y_names = 'DILI') 
	dev.off()





  # AUCs per reporter - all time points
	pdf(paste0(inputDir, '_nn_AUC_perReporter_allTimepoints_above1DMSO.pdf'), width = 3, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_above1DMSO.RData'))			# nn_1
		plotOutcomeNN(nn_1, 'AUCs')
		rm(nn_1)
	dev.off()

	pdf(paste0(inputDir, '_nn_AUC_perReporter_allTimepoints_above2DMSO.pdf'), width = 3, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_above2DMSO.RData'))			# nn_2
		plotOutcomeNN(nn_2, 'AUCs')
		rm(nn_2)
	dev.off()

	pdf(paste0(inputDir, '_nn_AUC_perReporter_allTimepoints_above3DMSO.pdf'), width = 3, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_above3DMSO.RData'))			# nn_3
		plotOutcomeNN(nn_3, 'AUCs')
		rm(nn_3)
	dev.off()

	pdf(paste0(inputDir, '_nn_AUC_perReporter_allTimepoints_above5DMSO.pdf'), width = 3, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_above5DMSO.RData'))			# nn_5
		plotOutcomeNN(nn_5, 'AUCs')
		rm(nn_5)
	dev.off()

	pdf(paste0(inputDir, '_nn_AUC_perReporter_allTimepoints_allParameters.pdf'), width = 3, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_allParameters.RData'))			# nn_CMAX
		plotOutcomeNN(nn_CMAX, 'AUCs')
		rm(nn_CMAX)
	dev.off()


  # AUCs - reporters & timepoints combined 
  ## !! TEST SET NOT PROPERLY INDEXED !! ##
	pdf(paste0(inputDir, '_nn_AUC_combinedReporter_allTimepoints_above1DMSO.pdf'), width = 1.5, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_above1DMSO_reportersCombined.RData'))			# nn_ALL1
		plotOutcomeNN(nn_ALL1, 'AUCs')
		rm(nn_ALL1)
	dev.off()

	pdf(paste0(inputDir, '_nn_AUC_combinedReporter_allTimepoints_above2DMSO.pdf'), width = 1.5, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_above2DMSO_reportersCombined.RData'))			# nn_ALL2
		plotOutcomeNN(nn_ALL2, 'AUCs')
		rm(nn_ALL2)
	dev.off()

	pdf(paste0(inputDir, '_nn_AUC_combinedReporter_allTimepoints_above3DMSO.pdf'), width = 1.5, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_above3DMSO_reportersCombined.RData'))			# nn_ALL3
		plotOutcomeNN(nn_ALL3, 'AUCs')
		rm(nn_ALL3)
	dev.off()

	pdf(paste0(inputDir, '_nn_AUC_combinedReporter_allTimepoints_above5DMSO.pdf'), width = 1.5, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_above5DMSO_reportersCombined.RData'))			# nn_ALL5
		plotOutcomeNN(nn_ALL5, 'AUCs')
		rm(nn_ALL5)
	dev.off()

	pdf(paste0(inputDir, '_nn_AUC_combinedReporter_allTimepoints_allParameters.pdf'), width = 1.5, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_allCMAXPars_reportersCombined.RData'))			# nn_ALL_CMAX
		plotOutcomeNN(nn_ALL_CMAX, 'AUCs')
		rm(nn_ALL_CMAX)
	dev.off()



	

  # AUCs - 24H - reporters & timepoints combined 
  ## !! TEST SET NOT PROPERLY INDEXED !! ##
  	pdf(paste0(inputDir, '_nnOutput_AUC_combinedReporter_24H_above2DMSO.pdf'), width = 1.5, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_24H_above2DMSO.RData'))			# nn_24H
		plotOutcomeNN(nn_24H, 'AUCs')
		rm(nn_24H)
	dev.off()

	pdf(paste0(inputDir, '_nnOutput_AUC_combinedReporter_48H_above2DMSO.pdf'), width = 1.5, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_48H_above2DMSO.RData'))			# nn_48H
		plotOutcomeNN(nn_48H, 'AUCs')
		rm(nn_48H)
	dev.off()

	pdf(paste0(inputDir, '_nnOutput_AUC_combinedReporter_72H_above2DMSO.pdf'), width = 1.5, height = 5)
		load(paste0(inputDir, 'nnOutput_50i_72H_above2DMSO.RData'))			# nn_72H
		plotOutcomeNN(nn_72H, 'AUCs')
		rm(nn_72H)
	dev.off()


  # 


## Neuralnetting - local plotting of trained networks
  # Example plot for graphical representation of NNs
	# load(paste0(inputDir, 'nnOutput_50i_above3DMSO.RData'))	
	

	pdf(paste0(inputDir, 'Example_nnFits_45inputNodes.pdf'), width = 14, height = 8)
		for(i in 1:2) {
			tmp <- nn_48H[[i]]$trainingFit	# 45 input nodes
			# tmp <- nn_3[[i]]$trainingFit	# 18 input nodes

			cat(i, '\n')
			plotnet(tmp, 
					alpha = 0.25, 
					pos_col = 'steelblue3', 
					neg_col = 'orangered3', 
					cex_val = 0.3,
					circle_cex = 2,
					max_sp = FALSE,
					x_names =  gsub('\\.', ' - ', toupper(do.call(rbind, strsplit(tmp$model.list$variables, '_'))[, 1])),
					y_names = 'DILI') 
		}
	dev.off()


## TOX PPT Mon 20 Feb 2017 / SOT poster
	## Settings™
	options(stringsAsFactors = FALSE)
	
	library(neuralnet)
	library(nnetpredint)
	library(NeuralNetTools)
	library(caret)

	library(ggplot2)

	library(data.table)
	library(preprocessCore)
	library(reshape2)
	library(parallel)

	inputDir <- '/Users/Wouter/stack/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'
	load(paste0(inputDir, 'nnOutput_50i_above2DMSO.RData'))

	tmp <- nn_2[[1]]
	tmp$call





	## 1) imabalnce
	## oversamplen kleine class cases (maakt gewicht groter van kleine class cases)
	## oversampling

	##- 2) multitask networks
	## 

	## 3) 70/15/15 train/test/validation

	## 4) scaling
	## min-max / zero-mean (variance)

	## 5) loss function
	## evaluate loss during training in test/validation set (per epoch)

	## MCC mathhews correlation coefficient
	## confusion matrix mee nemen 

	## act.fct
	## logistic niet cool (sigmoid beter idee OF _softmax_)

	## softmax
	## 2 output nodes: 1 voor DILI (0-1) en 1 voor nonDILI (0-1) -> samen 1

	## Dropout
	## 


	## Welke optimizer? rprop+
	## stochastic gradient decent with momentum 
	## 

	## learning rate
	## pytion CC LAM: 0.0001


	## eerste layer groter dan input layer -> 2e layer subst. kleiner (gooit onnodige info weg)


	## batch size
	## class imbalance door batch vorming tijdens trainen




























