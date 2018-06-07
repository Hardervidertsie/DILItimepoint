## Parser functions for timepoint DILI screen
## All function accept a single replicate from outputList

## Settings



	options(stringsAsFactors = FALSE)
	
	library(data.table)
	library(preprocessCore)
	library(parallel)

## Functions
	# meta column parser
	.parseMetaColumns <- function(dt, mc, verbose = TRUE, .debug = FALSE,...) {
		if(.debug) {
			verbose <- FALSE
			dt <- outputList[[1]]$sumData
		}

		# Treatment names
			if(verbose) cat('Parseing treatments...\n')
			dt$treatment <- sub('_', '', toupper(sub(' ', '', dt$treatment, fixed = TRUE)), fixed = TRUE)
			if(verbose) cat('1 of 4 steps done.\n')

		# Standardize plateID in sumData
			if(verbose) cat('Parseing plateIDs...\n')
			dt$plateID <- gsub('HSPA1B_', '', dt$plateID)
			dt$plateID <- toupper(sub('_wells', '', dt$plateID, fixed = TRUE))
			if(verbose) cat('2 of 4 steps done.\n')

			
		# Add repID, CMAX & CMAX_numeric
			if(verbose) cat('Parseing repIDs and CMAXs...\n')
			plateID_splitted <- do.call(rbind, strsplit(dt$plateID, '_'))
			plateID_splitted <- plateID_splitted[, !grepl('HMOX1', plateID_splitted[1, ])]		# HMOX1 platelayout mess

				colnames(plateID_splitted) <- c('a', 'b', 'c')
				colnames(plateID_splitted)[grep('REP',  plateID_splitted[1, ])]  <- 'repID'
				colnames(plateID_splitted)[grep('CMAX',  plateID_splitted[1, ])] <- 'CMAX'
				colnames(plateID_splitted)[grep('H',  plateID_splitted[1, ])]    <- 'TIME'


dt$repID <- as.numeric(sub('REP', '', plateID_splitted[, 'repID'], fixed = TRUE))
			dt$CMAX  <- plateID_splitted[, 'CMAX']
			dt$CMAX_numeric <- as.numeric(sub('CMAX', '', dt$CMAX, fixed = TRUE))
			if(verbose) cat('3 of 4 steps done.\n')

		# Add columns to match against
			if(verbose) cat('Creating columns to match against...\n')
			dt$matchDose <- paste(dt$treatment, dt$dose_uM, dt$CMAX, sep = '_')
			dt$matchDoseTime <- paste(dt$matchDose, dt$timeID, sep = '_')
			dt$uniqueID <- paste(dt$matchDoseTime, dt$repID, sep = '_')

				# legacy uniqueID
				# dt$uniqueID <- paste(dt$treatment, dt$timeID, dt$dose_uM, dt$plateID, sep = '_')
			if(verbose) cat('4 of 4 steps done.\n')

		# return
			return(dt)	
	}

	# Add mean DMSO for all plateIDs
	.addSummaryDMSO <- function(myDT, sumData, mc, verbose = FALSE, .debug = FALSE) {
		if(.debug) {
			verbose <- FALSE

			sumData <- pSD
			myDT <- pDT
		}

		# Pre-run checks/parses
			duplicatedEvents <- which(duplicated(sumData$uniqueID))
			if(length(duplicatedEvents) > 0) sumData <- sumData[-duplicatedEvents, ]
			
			if(!identical(unique(myDT$uniqueID), sumData$uniqueID)) stop('myDT and sumData do not seem to match.\n')

			if(!any(colnames(sumData) == 'numberOfObjects')) sumData$numberOfObjects <- table(myDT$uniqueID)[sumData$uniqueID]
		
		# Objects to assess DMSO in myDT
			featureGFP <- grep('GFP', grep('_Intensity_Integrated', colnames(myDT), value = TRUE), value = TRUE)
			plateIDs <- sumData$plateID

			indDMSO <- which(myDT$treatment == 'DMSO')

		# Acquire mean/median DMSO for each plateID
			if(verbose) cat('Acquiring GFP intensities per treatment...')

			if(mc) { 
				DMSOstats <- do.call(rbind, mclapply(plateIDs, function(.plateID) {
					.ind <- intersect(which(myDT$plateID == .plateID), indDMSO)
					c(mean(myDT[.ind, ][[featureGFP]], na.rm = TRUE), median(myDT[.ind, ][[featureGFP]], na.rm = TRUE) )
				}, mc.cores = 1, mc.preschedule = FALSE))
			} else {
				DMSOstats <- do.call(rbind, lapply(plateIDs, function(.plateID) {
					.ind <- intersect(which(myDT$plateID == .plateID), indDMSO)
					c(mean(myDT[.ind, ][[featureGFP]], na.rm = TRUE), median(myDT[.ind, ][[featureGFP]], na.rm = TRUE) )
				} ))
			}	

			rownames(DMSOstats) <- plateIDs
			colnames(DMSOstats) <- c('meanDMSO', 'medianDMSO')

		# Add meanDMSO & medianDMSO to myDT and sumData objects
			sumData$plateMeanDMSO <- DMSOstats[sumData$plateID, 'meanDMSO']	
			sumData$plateMedianDMSO <- DMSOstats[sumData$plateID, 'medianDMSO']	

			myDT$plateMeanDMSO <- DMSOstats[myDT$plateID, 'meanDMSO']	
			myDT$plateMedianDMSO <- DMSOstats[myDT$plateID, 'medianDMSO']	

			if(verbose) cat(' Done.\n')


		# Add booleans to myDT that indicates whether an entry (i.e. cell or nucleus) is above the DMSO plate mean/median.
			if(verbose) cat('Calculating fraction of objects above DMSO thresholds...')

			myDT$Above_1_meanIntregated_plateDMSO   <- myDT[[featureGFP]] > myDT$plateMeanDMSO
			myDT$Above_1_medianIntregated_plateDMSO <- myDT[[featureGFP]] > myDT$plateMedianDMSO

			myDT$Above_2_meanIntregated_plateDMSO   <- myDT[[featureGFP]] > (2 * myDT$plateMeanDMSO)
			myDT$Above_2_medianIntregated_plateDMSO <- myDT[[featureGFP]] > (2 * myDT$plateMedianDMSO)

			myDT$Above_3_meanIntregated_plateDMSO   <- myDT[[featureGFP]] > (3 * myDT$plateMeanDMSO)
			myDT$Above_3_medianIntregated_plateDMSO <- myDT[[featureGFP]] > (3 * myDT$plateMedianDMSO)

			myDT$Above_5_meanIntregated_plateDMSO   <- myDT[[featureGFP]] > (5 * myDT$plateMeanDMSO)
			myDT$Above_5_medianIntregated_plateDMSO <- myDT[[featureGFP]] > (5 * myDT$plateMedianDMSO)


		# Add number of nuclei/cells above DMSO for each treatment to sumData
			if(mc) {
				objectsAboveDMSO <- do.call(rbind, mclapply(sumData$uniqueID, function(.uniqueID) {
				  c(	
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_1_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_1_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_2_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_2_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_3_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_3_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_5_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_5_medianIntregated_plateDMSO']], na.rm = TRUE)) # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']])
				  )	
				}, mc.cores = 1, mc.preschedule = FALSE))	
			} else {
				objectsAboveDMSO <- do.call(rbind, lapply(sumData$uniqueID, function(.uniqueID) {
				  c(	
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_1_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_1_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_2_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_2_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_3_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_3_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_5_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Above_5_medianIntregated_plateDMSO']], na.rm = TRUE)) # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']])
				  )	
				}))	
			}	

			colnames(objectsAboveDMSO) <- paste0('Objects_', grep('Above_', colnames(myDT), value = TRUE))


		# Calclate positive fraction of nuclei/cells for DMSO thresholds
			fractionObjectsAboveDMSO <- apply(objectsAboveDMSO, 2, function(x) x/sumData$numberOfObjects)
				colnames(fractionObjectsAboveDMSO) <- gsub('Objects', 'Fraction', colnames(fractionObjectsAboveDMSO))

			if(verbose) cat(' Done.')
		# Parse
			sumData <- cbind(sumData, objectsAboveDMSO, fractionObjectsAboveDMSO)

		# Return
			return(list(myDT = myDT, sumData = sumData))
	}

	# ICAM1 - allow for lower GFP compared to DMSO-TNFa
	.addSummaryDMSO_rev <- function(myDT, sumData, mc, verbose = FALSE, .debug = FALSE) {
		if(.debug) {
			verbose <- FALSE

			sumData <- pSD
			myDT <- pDT
		}

		# Pre-run checks/parses
			duplicatedEvents <- which(duplicated(sumData$uniqueID))
			if(length(duplicatedEvents) > 0) sumData <- sumData[-duplicatedEvents, ]
			
			if(!identical(unique(myDT$uniqueID), sumData$uniqueID)) stop('myDT and sumData do not seem to match.\n')

			if(!any(colnames(sumData) == 'numberOfObjects')) sumData$numberOfObjects <- table(myDT$uniqueID)[sumData$uniqueID]
		
		# Objects to assess DMSO in myDT
			featureGFP <- grep('GFP', grep('_Intensity_Integrated', colnames(myDT), value = TRUE), value = TRUE)
			plateIDs <- sumData$plateID

			indDMSO <- which(myDT$treatment == 'DMSO')

		# Acquire mean/median DMSO for each plateID
			if(verbose) cat('Acquiring GFP intensities per treatment...')

			if(mc) { 
				DMSOstats <- do.call(rbind, mclapply(plateIDs, function(.plateID) {
					.ind <- intersect(which(myDT$plateID == .plateID), indDMSO)
					c(mean(myDT[.ind, ][[featureGFP]], na.rm = TRUE), median(myDT[.ind, ][[featureGFP]], na.rm = TRUE) )
				}, mc.cores = 1, mc.preschedule = FALSE))
			} else {
				DMSOstats <- do.call(rbind, lapply(plateIDs, function(.plateID) {
					.ind <- intersect(which(myDT$plateID == .plateID), indDMSO)
					c(mean(myDT[.ind, ][[featureGFP]], na.rm = TRUE), median(myDT[.ind, ][[featureGFP]], na.rm = TRUE) )
				} ))
			}	

			rownames(DMSOstats) <- plateIDs
			colnames(DMSOstats) <- c('meanDMSO', 'medianDMSO')

		# Add meanDMSO & medianDMSO to myDT and sumData objects
			sumData$plateMeanDMSO <- DMSOstats[sumData$plateID, 'meanDMSO']	
			sumData$plateMedianDMSO <- DMSOstats[sumData$plateID, 'medianDMSO']	

			myDT$plateMeanDMSO <- DMSOstats[myDT$plateID, 'meanDMSO']	
			myDT$plateMedianDMSO <- DMSOstats[myDT$plateID, 'medianDMSO']	

			if(verbose) cat(' Done.\n')


		# Add booleans to myDT that indicates whether an entry (i.e. cell or nucleus) is BELOW the DMSO plate mean/median.
			if(verbose) cat('Calculating fraction of objects below DMSO thresholds...')

			myDT$Below_1_meanIntregated_plateDMSO   <- myDT[[featureGFP]] < myDT$plateMeanDMSO
			myDT$Below_1_medianIntregated_plateDMSO <- myDT[[featureGFP]] < myDT$plateMedianDMSO

			myDT$Below_2_meanIntregated_plateDMSO   <- myDT[[featureGFP]] < ((1/2) * myDT$plateMeanDMSO)
			myDT$Below_2_medianIntregated_plateDMSO <- myDT[[featureGFP]] < ((1/2) * myDT$plateMedianDMSO)

			myDT$Below_3_meanIntregated_plateDMSO   <- myDT[[featureGFP]] < ((1/3) * myDT$plateMeanDMSO)
			myDT$Below_3_medianIntregated_plateDMSO <- myDT[[featureGFP]] < ((1/3) * myDT$plateMedianDMSO)

			myDT$Below_5_meanIntregated_plateDMSO   <- myDT[[featureGFP]] < ((1/5) * myDT$plateMeanDMSO)
			myDT$Below_5_medianIntregated_plateDMSO <- myDT[[featureGFP]] < ((1/5) * myDT$plateMedianDMSO)


		# Add number of nuclei/cells BELOW DMSO for each treatment to sumData
			if(mc) {
				objectsBelowDMSO <- do.call(rbind, mclapply(sumData$uniqueID, function(.uniqueID) {
				  c(	
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_1_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_1_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_2_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_2_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_3_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_3_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_5_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_5_medianIntregated_plateDMSO']], na.rm = TRUE)) # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']])
				  )	
				}, mc.cores = 1, mc.preschedule = FALSE))	
			} else {
				objectsBelowDMSO <- do.call(rbind, lapply(sumData$uniqueID, function(.uniqueID) {
				  c(	
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_1_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_1_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_2_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_2_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_3_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_3_medianIntregated_plateDMSO']], na.rm = TRUE)), # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
				
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_5_meanIntregated_plateDMSO']], na.rm = TRUE)), #   / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']]),
					(sum(myDT[myDT$uniqueID == .uniqueID, ][['Below_5_medianIntregated_plateDMSO']], na.rm = TRUE)) # / (sumData[sumData$uniqueID == .uniqueID][['numberOfObjects']])
				  )	
				}))	
			}	

			colnames(objectsBelowDMSO) <- paste0('Objects_', grep('Below_', colnames(myDT), value = TRUE))


		# Calclate positive fraction of nuclei/cells for DMSO thresholds
			fractionObjectsBelowDMSO <- apply(objectsBelowDMSO, 2, function(x) x/sumData$numberOfObjects)
				colnames(fractionObjectsBelowDMSO) <- gsub('Objects', 'Fraction', colnames(fractionObjectsBelowDMSO))

			if(verbose) cat(' Done.')
		# Parse
			sumData <- cbind(sumData, objectsBelowDMSO, fractionObjectsBelowDMSO)

		# Return
			return(list(myDT = myDT, sumData = sumData))
	}

	# Add AnV 
	.addAnV <- function(myDT, sumData, mc, .debug = FALSE) {
		if(.debug) {
			verbose <- FALSE

			sumData <- pSD
			myDT <- pDT
		}

		# Pre-run checks
			if(!identical(unique(myDT$uniqueID), sumData$uniqueID)) stop('myDT and sumData do not seem to match.\n')
			if(!any(colnames(sumData) == 'numberOfObjects')) sumData$numberOfObjects <- table(myDT$uniqueID)[sumData$uniqueID]

		# Add total number of objects with non-NA AnV, fraction of non-NA AnV objects and mean/median of non-NA AnV objects
			if(mc) {	
				AnVstats <- do.call(rbind, mclapply(sumData$uniqueID, function(.uniqueID) {
					.AnV <- myDT$AnV_masked_primaryID_AreaShape_Area[which(myDT$uniqueID == .uniqueID)]
					
					objectsAnV  <- sum(!(is.na(.AnV)))
					fractionAnV <- objectsAnV / length(.AnV)

					meanAnV   <- mean(.AnV, na.rm = TRUE)
					medianAnV <- median(.AnV, na.rm = TRUE)

					c(.uniqueID, objectsAnV, fractionAnV, meanAnV, medianAnV)

				}, mc.cores = 1, mc.preschedule = FALSE))
			} else {
				AnVstats <- do.call(rbind, lapply(sumData$uniqueID, function(.uniqueID) {
					.AnV <- myDT$AnV_masked_primaryID_AreaShape_Area[which(myDT$uniqueID == .uniqueID)]
					
					objectsAnV  <- sum(!(is.na(.AnV)))
					fractionAnV <- objectsAnV / length(.AnV)

					meanAnV   <- mean(.AnV, na.rm = TRUE)
					medianAnV <- median(.AnV, na.rm = TRUE)

					c(.uniqueID, objectsAnV, fractionAnV, meanAnV, medianAnV)

				} ))
			}	

			colnames(AnVstats) <- c('uniqueID', 'AnV_numberNonNA', 'AnV_fractionNonNA', 'AnV_meanNonNA', 'AnV_medianNonNA')

		# Parse
			if(identical(AnVstats[, 'uniqueID'], sumData$uniqueID)) sumData <- cbind(sumData, AnVstats)		

			
		# Return
			return(sumData)
	}

	# Add PI
	.addPI <- function(myDT, sumData, mc, .debug = FALSE) {
		if(.debug) {
			verbose <- FALSE

			sumData <- pSD
			myDT <- pDT
		}

		# Pre-run checks
			if(!identical(unique(myDT$uniqueID), sumData$uniqueID)) stop('myDT and sumData do not seem to match.\n')
			if(!any(colnames(sumData) == 'numberOfObjects')) sumData$numberOfObjects <- table(myDT$uniqueID)[sumData$uniqueID]

		# Add total number of objects with non-NA PI, fraction of non-NA PI objects and mean/median of non-NA PI objects
			if(mc) {	
				PIstats <- do.call(rbind, mclapply(sumData$uniqueID, function(.uniqueID) {
					.PI <- myDT$PI_masked_primaryID_AreaShape_Area[which(myDT$uniqueID == .uniqueID)]
					
					objectsPI  <- sum(!(is.na(.PI)))
					fractionPI <- objectsPI / length(.PI)

					meanPI   <- mean(.PI, na.rm = TRUE)
					medianPI <- median(.PI, na.rm = TRUE)

					c(.uniqueID, objectsPI, fractionPI, meanPI, medianPI)

				}, mc.cores = 1, mc.preschedule = FALSE))
			} else {
				PIstats <- do.call(rbind, lapply(sumData$uniqueID, function(.uniqueID) {
					.PI <- myDT$PI_masked_primaryID_AreaShape_Area[which(myDT$uniqueID == .uniqueID)]
					
					objectsPI  <- sum(!(is.na(.PI)))
					fractionPI <- objectsPI / length(.PI)

					meanPI   <- mean(.PI, na.rm = TRUE)
					medianPI <- median(.PI, na.rm = TRUE)

					c(.uniqueID, objectsPI, fractionPI, meanPI, medianPI)

				} ))
			}	

			colnames(PIstats) <- c('uniqueID', 'PI_numberNonNA', 'PI_fractionNonNA', 'PI_meanNonNA', 'PI_medianNonNA')

		# Parse
			if(identical(PIstats[, 'uniqueID'], sumData$uniqueID)) sumData <- cbind(sumData, PIstats)		

		# Return
			return(sumData)
	}

	# Wrapper function
	outputListParser <- function(list2parse, ICAM1 = FALSE, MC = FALSE, subMC = FALSE, mcCores = 1, .debug = FALSE, ...) {
		if(.debug) {
			list2parse <- outputList
			subMC <- TRUE
			mcCores <- 6

			rep <- outputList[[1]]
		}

		# Pre-run checks

		# Add meta data
			cat('Adding meta data... ')
			if(MC) {	
				.outputList <- mclapply(list2parse, function(rep) {
					list(sumData = .parseMetaColumns(rep$sumData, mc = subMC), 
						 myDT    = .parseMetaColumns(rep$myDT, mc = subMC))
				}, mc.cores = 1, mc.preschedule = FALSE)
			} else {
				.outputList <- lapply(list2parse, function(rep) {
					list(sumData = .parseMetaColumns(rep$sumData, mc = subMC), 
						 myDT    = .parseMetaColumns(rep$myDT, mc = subMC))
				} )
			}					
			cat('Done!\n')

		# Add DMSO, PI and AnV summary statistics 
			cat('Adding DMSO, PI and AnV data... ')
				.outputList <- lapply(.outputList, function(rep) {
					cat('\nParseing replicate', rep$sumData$repID[1], '...')

					.dmso <- .addSummaryDMSO(rep$myDT, rep$sumData, mc = subMC)
					.anv  <- .addAnV(.dmso$myDT, .dmso$sumData, mc = subMC)
					.pi   <- .addPI(.dmso$myDT, .anv, mc = subMC)

					if(ICAM1) {
						.dmso_rev <- .addSummaryDMSO_rev(rep$myDT, rep$sumData, mc = subMC)
						.anv_rev  <- .addAnV(.dmso_rev$myDT, .dmso_rev$sumData, mc = subMC)
						.pi_rev   <- .addPI(.dmso_rev$myDT, .anv_rev, mc = subMC)

						list(sumData = .pi, 
						 	 myDT = .dmso$myDT,
						 	 sumData_rev = .pi_rev,
						 	 myDT_rev = .dmso_rev$myDT)

					} else {
						list(sumData = .pi, 
						 	 myDT = .dmso$myDT)
					}

					
				}) 

			cat('\nDone!\n')
		
		# Return
			return(.outputList)	
	}

## Parse reporter screens
	if(FALSE) {
	  inputDir <- 'D:/analysis/DILI timepoint/data/'
		
	  load(paste0(inputDir, 'BIP/', 'outputList BIP.Rdata'))
	  BIP <- outputListParser(list2parse = outputList, verbose = TRUE, subMC = FALSE, MC = FALSE, mcCores = 1)
	  save(BIP, file = paste0(inputDir, 'BIP/', 'parsedBIP.RData'))
rm('BIP')
rm('outputList')
gc(reset=TRUE)
	 
 load(paste0(inputDir, 'SRXN1/', 'outputList110516.Rdata'))
	  SRXN1 <- outputListParser(list2parse = outputList,verbose = TRUE, subMC = FALSE, MC = FALSE, mcCores = 1)
	  save(SRXN1, file = paste0(inputDir, 'SRXN1/', 'parsedSRXN1.RData'))
	  rm('SRXN1')
	  rm('outputList')
	  gc(reset=TRUE)
	  
	  load(paste0(inputDir, 'P21/', 'outputList12052016.Rdata'))
	  P21 <- outputListParser(list2parse = outputList, verbose = TRUE, subMC = FALSE, MC = FALSE, mcCores = 1)
	  save(P21, file = paste0(inputDir,'P21/', 'parsedP21.RData'))	
	  rm('P21')
	  rm('outputList')
	  gc(reset=TRUE)
	  
	  load(paste0(inputDir,'ICAM1/', 'outputListICAM1.Rdata'))
	  ICAM1 <- outputListParser(list2parse = outputList, ICAM1 = TRUE, verbose = TRUE, subMC = FALSE, MC = FALSE, mcCores = 1)
	  save(ICAM1, file = paste0(inputDir,'ICAM1/', 'parsedICAM1.RData'))	
	  rm('ICAM1')
	  rm('outputList')
	  gc(reset=TRUE)
	  
	  load(paste0(inputDir, 'CHOP/', 'outputList110516.Rdata'))
	  CHOP <- outputListParser(list2parse = outputList, verbose = TRUE, subMC = FALSE, MC = FALSE, mcCores = 1)
	  save(CHOP, file = paste0(inputDir,'CHOP/', 'parsedCHOP.RData'))			
	  rm('CHOP')
	  rm('outputList')
	  gc(reset=TRUE)
	  
	  load(paste0(inputDir, 'HSPA1B/', 'outputList HSPA1B.Rdata'))
	  HSPA1B <- outputListParser(list2parse = outputList, verbose = TRUE, subMC = FALSE, MC = FALSE, mcCores = 1)
	  save(HSPA1B, file = paste0(inputDir, 'HSPA1B/', 'parsedHSPA1B.RData'))	
	  rm('HSPA1B')
	  rm('outputList')
	  gc(reset=TRUE)
	  
	  load(paste0(inputDir, 'BTG2/', 'outputList.Rdata'))
	  BTG2 <- outputListParser(list2parse = outputList, verbose = TRUE, subMC = FALSE, MC = FALSE, mcCores = 1)
	  save(BTG2, file = paste0(inputDir, 'BTG2/', 'parsedBTG2.RData'))	
	  rm('BTG2')
	  rm('outputList')
	  gc(reset=TRUE)
	  
	  load(paste0(inputDir, 'HMOX1/', 'outputList.Rdata'))
	  HMOX1 <- outputListParser(list2parse = outputList, verbose = TRUE, subMC = FALSE, MC = FALSE, mcCores = 1)
	  save(HMOX1, file = paste0(inputDir, 'HMOX1/', 'parsedHMOX1.RData'))		

	}

