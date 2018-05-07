## First go at caret package & classifiers
## 08 June 2016 - Wouter den Hollander

## Variables to pounder about:
#	timeID
#	replicate merger (REP1 in CHOP != REP1 in ICAM1)
#	Control GFP substraction
# 	All data in single data.frame or seperate cell lines
# 	CMAX as numeric value

## Settings
	options(stringsAsFactors = FALSE)
	
	library(data.table)
	library(caret)

	inputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'


## Functions
	parseOutputList <- function(outputList, .debug = FALSE) {
		if(.debug) {
			load("/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/outputListCHOP.Rdata")
		}

		.outputList <- outputList
		names(.outputList) <- c('replicate_1', 'replicate_2', 'replicate_3')

		for(i in 1:length(.outputList)) {
			tmp <- .outputList[[i]]

			tmp$myDT$treatment <- gsub('_', '', toupper(gsub(' ', '', tmp$myDT$treatment)))
			
			tmp$sumData$treatment <- gsub('_', '', toupper(gsub(' ', '', tmp$sumData$treatment)))
			tmp$sumData$plateID <- toupper(gsub('_wells', '', tmp$sumData$plateID))
			tmp$sumData$uniqueID <- apply(tmp$sumData, 1, function(x) paste(x[-length(x)], collapse = '_'))

			plateID_splitted <- do.call(rbind, strsplit(tmp$sumData$plateID, '_'))
			cmaxInd <- grep('CMAX', toupper(plateID_splitted)[1, ])

			tmp$sumData$CMAX <- plateID_splitted[, cmaxInd]
			tmp$sumData$CMAX_numeric <- as.numeric(gsub('CMAX', '', tmp$sumData$CMAX))

			tmp$sumData$matchDose <- paste(tmp$sumData$treatment, tmp$sumData$dose_uM, tmp$sumData$CMAX, sep = '_')
			tmp$sumData$matchDoseTime <- paste(tmp$sumData$treatment, tmp$sumData$dose_uM, tmp$sumData$CMAX, tmp$sumData$timeID,sep = '_')

			.outputList[[i]] <- tmp
		}

		return(.outputList)
	}


## Data
	# load
	inputFiles <- list.files(inputDir, pattern = '.Rdata')

	dat <- datMeta <- list()
		for(i in 1:length(inputFiles)) {
			load(paste0(inputDir, inputFiles[[i]]))

			dat[[i]] <- parseOutputList(outputList)
			datMeta[[i]] <- outputList[[3]]$metaCSVData	
				names(dat)[[i]] <- names(datMeta)[[i]] <- gsub('outputList', '', gsub('.Rdata', '', inputFiles[[i]]))

			rm(outputList)
			cat(paste0('Done with ', names(dat)[[i]], '.\n'))
		}	

	# remove wells with too little/no compound pipetted from $sumData
	failedWells <- lapply(datMeta, function(x) x[x$replID > 3])
	failedWells <- lapply(failedWells, function(x) {
		x$treatment <- toupper(x$treatment)
		x
	})

	for(i in names(dat)) {
		.failedWells <- failedWells[[i]]
		.failedWellsInfo <- toupper(do.call(rbind, lapply(strsplit(as.character(.failedWells$plateID), '_'), function(x) x[1:3])))

		repInd <- grep('REP', .failedWellsInfo[1, ])
		cmaxInd <- grep('CMAX', .failedWellsInfo[1, ])

		for(j in 1:nrow(.failedWells)) {
			if(.failedWells$treatment[j] == 'DMEM') next

			.replicate <- gsub('REP', 'replicate_', .failedWellsInfo[j, repInd])
			.sumDat <- dat[[i]][[.replicate]]$sumData

			removeFromSumDat <- which(.sumDat$treatment == .failedWells$treatment[j] & .sumDat$CMAX == .failedWellsInfo[j, cmaxInd] & .sumDat$timeID == .failedWells$timeID[j])

			if(length(removeFromSumDat) != 1) stop('Something went wrong.\n')

		}

	}

## Parseing
	# While most of this could be done quicker using merge()/merge_all()/Reduce(),
	# the sparsity of overlap is somewhat scary and therefor a lot is done manually.

	# Generate merged sumData data.frame() per cell line
	ReporterSumData <- lapply(dat, function(x) lapply(x, function(y) y$sumData))
	ReporterSumData <- lapply(1:length(ReporterSumData), function(i) {
					reporter <- names(ReporterSumData)[i]
					cat(reporter, '\n')
					
					x <- ReporterSumData[[i]]
					df <- data.frame(Reduce(function(...) merge(..., all = TRUE, by = 'matchDoseTime'), x))
					
					# Generate single meta columns
					newMeta <- data.frame(sapply(c('treatment', 'timeID', 'dose_uM', 'CMAX', 'CMAX_numeric', 'matchDose'), function(y) {
						colsToMerge <- c(y, paste0(y, '.x'), paste0(y, '.y'))
						apply(df[, colsToMerge], 1, function(z) unique(na.omit(z)))
					}), matchDoseTime = df$matchDoseTime)

					
					# Sanity check mode: force replicate IDs to be fixed and relate to raw data (i.e. dat$reporter)
					extractedFeature <- strsplit(grep('_GFP', colnames(df), value = TRUE), '\\.')[[1]][1]

					repMergeIndices <- paste0(extractedFeature, gsub('plateID', '', 
										c(rep1 = grep('plateID', colnames(df), value = TRUE)[which(grepl('REP1', df[, grep('plateID', colnames(df))]))],
										  rep2 = grep('plateID', colnames(df), value = TRUE)[which(grepl('REP2', df[, grep('plateID', colnames(df))]))],
										  rep3 = grep('plateID', colnames(df), value = TRUE)[which(grepl('REP3', df[, grep('plateID', colnames(df))]))])
										))

					if(grepl('Nuclei', extractedFeature)) {
							df <- data.frame(newMeta, 
										  Nuclei_Intensity_MeanIntensity_image_GFP_REP1 = df[,repMergeIndices[1]],
										  Nuclei_Intensity_MeanIntensity_image_GFP_REP2 = df[,repMergeIndices[2]],
										  Nuclei_Intensity_MeanIntensity_image_GFP_REP3 = df[,repMergeIndices[3]]
										)
						} else {
							df <- data.frame(newMeta, 
										  Cytoplasm_Intensity_MeanIntensity_image_GFP_REP1 = df[,repMergeIndices[1]],
										  Cytoplasm_Intensity_MeanIntensity_image_GFP_REP2 = df[,repMergeIndices[2]],
										  Cytoplasm_Intensity_MeanIntensity_image_GFP_REP3 = df[,repMergeIndices[3]]
										)
						}

					colnames(df)[grep('_GFP', colnames(df))] <- paste0(reporter, '_', colnames(df)[grep('_GFP', colnames(df))])	
					
					# Take average of DMEM entries and create final data.frame()
					duplicatedDMEMs <- names(which(table(df$matchDoseTime) > 1))
					dfDMEM <- df[which(df$matchDoseTime %in% duplicatedDMEMs),]
					dfNoDMEM <- df[-which(df$matchDoseTime %in% duplicatedDMEMs),]

					dfDMEMAv <- dfDMEM[!duplicated(dfDMEM$matchDoseTime), ]
						for(y in dfDMEMAv$matchDoseTime) {
							dfDMEMAv[dfDMEMAv$matchDoseTime == y, grep('GFP', colnames(dfDMEMAv))] <- colMeans(dfDMEM[which(dfDMEM$matchDoseTime == y), grep('GFP', colnames(dfDMEM))])
						}

					# Create and parse final object
					finalDF <- rbind(dfNoDMEM, dfDMEMAv)
						finalDF$timeID <- as.numeric(finalDF$timeID)
						finalDF$dose_uM <- as.numeric(finalDF$dose_uM)
						finalDF$CMAX_numeric <- as.numeric(finalDF$CMAX_numeric)

					# finalDF[which(finalDF$treatment == 'DMSO' & finalDF$dose_uM > 0 & finalDF$timeID == 48), ]

					return(finalDF)
				})
		names(ReporterSumData) <- names(dat)

	# Generate single data.frame() which contains all data: replicates seperate
	ReporterSumDataAll <- data.frame(matchDoseTime = unique(unlist(sapply(ReporterSumData, function(x) x$matchDoseTime))))
		tmp <- do.call(rbind, strsplit(ReporterSumDataAll$matchDoseTime, '_'))
		ReporterSumDataAll$dose_uM <- as.numeric(tmp[,2])
		ReporterSumDataAll$timeID <- as.numeric(tmp[,4])
		ReporterSumDataAll$CMAX_numeric <- as.numeric(gsub('CMAX', '', tmp[,3]))
		ReporterSumDataAll$treatment <- tmp[,1]

		for(i in 1:nrow(ReporterSumDataAll)) {
			.matchDoseTime <- ReporterSumDataAll$matchDoseTime[i]

			GFP <- lapply(ReporterSumData, function(x) x[x$matchDoseTime == .matchDoseTime, grep('GFP', colnames(x))])
			if(i == 1) ReporterSumDataAll[, sapply(GFP, colnames)] <- NA

			for(j in 1:length(GFP)) {
				.cellData <- unlist(GFP[[j]])
				ReporterSumDataAll[i, names(.cellData)] <- .cellData
			}
		}	
	ReporterSumDataAll <- ReporterSumDataAll[rowSums(is.na(ReporterSumDataAll)) < 9, ]	# removed measurements only present in single cell line



	# Generate single data.frame() which contains all data: replicates merged
	CHOP_Nuclei_Intensity_MeanIntensity_image_GFP_REPMEAN <- rowMeans(ReporterSumDataAll[, grep('CHOP', colnames(ReporterSumDataAll), value = TRUE)], na.rm = TRUE)
	ICAM1_Cytoplasm_Intensity_MeanIntensity_image_GFP_REPMEAN <- rowMeans(ReporterSumDataAll[, grep('ICAM1', colnames(ReporterSumDataAll), value = TRUE)], na.rm = TRUE)
	P21_Nuclei_Intensity_MeanIntensity_image_GFP_REPMEAN <- rowMeans(ReporterSumDataAll[, grep('P21', colnames(ReporterSumDataAll), value = TRUE)], na.rm = TRUE)
	SRXN1_Cytoplasm_Intensity_MeanIntensity_image_GFP_REPMEAN <- rowMeans(ReporterSumDataAll[, grep('SRXN1', colnames(ReporterSumDataAll), value = TRUE)], na.rm = TRUE)

	ReporterSumDataAllMeans <- data.frame(ReporterSumDataAll[, 1:5], 
										  CHOP_Nuclei_Intensity_MeanIntensity_image_GFP_REPMEAN,
										  ICAM1_Cytoplasm_Intensity_MeanIntensity_image_GFP_REPMEAN,
										  P21_Nuclei_Intensity_MeanIntensity_image_GFP_REPMEAN,
										  SRXN1_Cytoplasm_Intensity_MeanIntensity_image_GFP_REPMEAN)

	# Final set of objects
	DILI <- list(raw = dat, 
				 rawMergedReplicateSumData = ReporterSumData, 
				 parsedDat = ReporterSumDataAll,
				 parsedDatMeans = ReporterSumDataAllMeans)

	save(DILI, file = paste0(inputDir, 'DILI_08June2016.RData'))


## Fitting: gbm
	# Replicates merged
		# Use all data in single data.frame() for prediction

		# Use seperate reporters for prediction

	# Replicates seperate
		# Use all data in single data.frame() for prediction

		# Use seperate reporters for prediction
			
