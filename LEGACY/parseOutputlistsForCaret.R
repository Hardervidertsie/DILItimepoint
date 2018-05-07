## Parse reporter lines outputlist files for input into caret package
## 15 June 2016 - Wouter den Hollander

## Settings
	options(stringsAsFactors = FALSE)
	library(data.table)
	library(preprocessCore)

	inputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'
	inputFiles <- list.files(inputDir, pattern = '.Rdata')


## Functions
	parseOutputListSumData <- function(outputList, .debug = FALSE) {
		if(.debug) {
			load("/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/outputListHSPA1B.Rdata")
			.outputList <- outputList
		}

		# Internal parser function
		.parser <- function(dt) {
			# Treatment names
			dt$treatment <- gsub('_', '', toupper(gsub(' ', '', dt$treatment)))

			# Standardize plateID in sumData
			dt$plateID <- toupper(gsub('_wells', '', dt$plateID))

			# Create uniqueID in sumData
			dt$uniqueID <- paste(dt$treatment, dt$timeID, dt$dose_uM, dt$plateID, sep = '_')

			# Split uniqueID and create CMAX columns
			uniqueID_splitted <- do.call(rbind, strsplit(dt$uniqueID, '_'))
			dt$CMAX <- uniqueID_splitted[ ,grep('CMAX', toupper(uniqueID_splitted)[1, ])]
			dt$CMAX_numeric <- as.numeric(gsub('CMAX', '', dt$CMAX))

			# Create columns to match against
			dt$matchDose <- paste(dt$treatment, dt$dose_uM, dt$CMAX, sep = '_')
			dt$matchDoseTime <- paste(dt$treatment, dt$dose_uM, dt$CMAX, dt$timeID,sep = '_')

			# return
			return(dt)
		}

		# Assign names to replicates
		names(outputList) <- c('replicate_1', 'replicate_2', 'replicate_3')

		# Pipetting issues (replID > 3)
		pipIssues <- outputList$replicate_3$metaCSVData[which(outputList$replicate_3$metaCSVData$replID > 3), ]
		if(nrow(pipIssues) > 0) pipIssues <- .parser(pipIssues)

		# outputList parse
		parsedList <- lapply(outputList, function(x) .parser(x$sumData))

		# Flag pipetting issues in parsedList
		parsedList$replicate_1$pipIssue <- parsedList$replicate_1$uniqueID %in% pipIssues$uniqueID
		parsedList$replicate_2$pipIssue <- parsedList$replicate_2$uniqueID %in% pipIssues$uniqueID
		parsedList$replicate_3$pipIssue <- parsedList$replicate_3$uniqueID %in% pipIssues$uniqueID

		# Return
		return(parsedList)
	}

	combineReplicates <- function(parsedList, .debug = FALSE) {
		if(.debug) {
			parsedList <- dat$HSPA1B
			.parsedList <- parsedList
		}

		# Get all unique matchDoseTime 
		matchDoseTime <- sort(unique(unlist(lapply(parsedList, function(x) x$matchDoseTime))))
		
		# Create feature data.frame
		features <- data.frame(do.call(rbind, strsplit(matchDoseTime, '_')))
			colnames(features) <- c('treatment', 'dose_uM', 'CMAX', 'timeID')

			features$dose_uM <- as.numeric(features$dose_uM)
			features$timeID <- as.numeric(features$timeID)
			features$CMAX_numeric <- as.numeric(gsub('CMAX', '', features$CMAX))
			features$matchDoseTime <- matchDoseTime

		# Get GFP data 
		REP1_pipIssue <- parsedList$replicate_1[['pipIssue']][match(features$matchDoseTime, parsedList$replicate_1$matchDoseTime)]
		REP1_GFP <- parsedList$replicate_1[[grep('GFP', colnames(parsedList$replicate_1), value = TRUE)]][match(features$matchDoseTime, parsedList$replicate_1$matchDoseTime)]
			
		REP2_pipIssue <- parsedList$replicate_2[['pipIssue']][match(features$matchDoseTime, parsedList$replicate_2$matchDoseTime)]
		REP2_GFP <- parsedList$replicate_2[[grep('GFP', colnames(parsedList$replicate_2), value = TRUE)]][match(features$matchDoseTime, parsedList$replicate_2$matchDoseTime)]
		
		REP3_pipIssue <- parsedList$replicate_3[['pipIssue']][match(features$matchDoseTime, parsedList$replicate_3$matchDoseTime)]
		REP3_GFP <- parsedList$replicate_3[[grep('GFP', colnames(parsedList$replicate_3), value = TRUE)]][match(features$matchDoseTime, parsedList$replicate_3$matchDoseTime)]

		# Normalised GFP data
		GFP <- cbind(REP1_GFP, REP2_GFP, REP3_GFP)
		GFP <- apply(GFP, 2, function(x) {
			x[is.nan(x)] <- NA
			x
		})


		quantileGFP <- normalize.quantiles(GFP)
		zScoreGFP <- apply(GFP, 2, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
		minMaxGFP <- apply(GFP, 2, function(x) (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

		colnames(GFP) <- paste0('rawGFP_', c('REP1', 'REP2', 'REP3'))
		colnames(quantileGFP) <- paste0('quantileGFP_', c('REP1', 'REP2', 'REP3'))
		colnames(zScoreGFP) <- paste0('zScoreGFP_', c('REP1', 'REP2', 'REP3'))
		colnames(minMaxGFP) <- paste0('minMaxGFP_', c('REP1', 'REP2', 'REP3'))

		# Create average replicates column
		rawGFP_mean <- rowMeans(GFP, na.rm = TRUE)
		quantileGFP_mean <- rowMeans(quantileGFP, na.rm = TRUE)
		zScoreGFP_mean <- rowMeans(zScoreGFP, na.rm = TRUE)
		minMaxGFP_mean <- rowMeans(minMaxGFP, na.rm = TRUE)

		# Final object
		df <- data.frame(features, 
						 GFP, rawGFP_mean, 
						 quantileGFP, quantileGFP_mean,
						 zScoreGFP, zScoreGFP_mean,
						 minMaxGFP, minMaxGFP_mean,
						 REP1_pipIssue, REP2_pipIssue, REP3_pipIssue)

		return(df)
	}


## Load & parse outputLists
	for(i in 1:length(inputFiles)) {
		if(!exists('dat')) dat <- list()

		load(paste0(inputDir, inputFiles[[i]]))

		dat[[i]] <- parseOutputListSumData(outputList)
		names(dat)[[i]] <- gsub('outputList', '', gsub('.Rdata', '', inputFiles[[i]]))

		rm(outputList)
		cat(paste0('Loaded ', names(dat)[[i]], '.\n'))
	}	


## Check discrepancies of parsed data
	sapply(dat, function(x) sapply(x, nrow))
	  #		           BIP CHOP HSPA1B ICAM1  P21 SRXN1
	  #	replicate_1   3528 3510   3332  2348 3510  3290
	  #	replicate_2   3528 3510   3528  2348 3510  3120
	  #	replicate_3   3528 3510   3528  2348 3510  3504

	nonComplete <- lapply(dat, function(x) which(table(unlist(lapply(x, function(y) y$matchDoseTime))) != 3))
	sapply(nonComplete, length)
   	  #	BIP   CHOP HSPA1B  ICAM1    P21  SRXN1 
      #	  0     18    196     12     18    628 
		
	issuesCHOP <- lapply(dat$CHOP, function(x) x[x$matchDoseTime %in% names(nonComplete$CHOP)])		# DMEM multiples
	issuesICAM1 <- lapply(dat$ICAM1, function(x) x[x$matchDoseTime %in% names(nonComplete$ICAM1)])	# DMEM multiples
	issuesP21 <- lapply(dat$P21, function(x) x[x$matchDoseTime %in% names(nonComplete$P21)])		# DMEM multiples

	issuesHSPA1B <- lapply(dat$HSPA1B, function(x) x[x$matchDoseTime %in% names(nonComplete$HSPA1B)])	# $replicate_1 missings
	issuesSRXN1 <- lapply(dat$SRXN1, function(x) x[x$matchDoseTime %in% names(nonComplete$SRXN1)])		# All sorts of known issues


## Generate single data.frame for each reporter
	datCombined <- lapply(dat, combineReplicates)
	
	save(datCombined, file = paste0(inputDir, 'parsedOutputLists.RData'))

## Generate single data.frame containing all data
	# Get all unique matchDoseTime 
	matchDoseTime <- sort(unique(unlist(lapply(datCombined, function(x) x$matchDoseTime))))
		
	# Create feature data.frame
	features <- data.frame(do.call(rbind, strsplit(matchDoseTime, '_')))
		colnames(features) <- c('treatment', 'dose_uM', 'CMAX', 'timeID')

		features$dose_uM <- as.numeric(features$dose_uM)
		features$timeID <- as.numeric(features$timeID)
		features$CMAX_numeric <- as.numeric(gsub('CMAX', '', features$CMAX))
		features$matchDoseTime <- matchDoseTime

	quantileGFP <- sapply(datCombined, function(x) {
		x[[grep('quantileGFP_mean', colnames(x), value = TRUE)]][match(features$matchDoseTime, x$matchDoseTime)]
	})










