## Parse reporter lines outputlist files for input into caret package
## 06 July 2016 - Wouter den Hollander

## TO DO
##	Results Anke
##		Fraction > 3x DMSO, Integrated Intensity, Cell/Nuclear count, PI fraction, AnV fraction

## Settings
	options(stringsAsFactors = FALSE)
	
	library(data.table)
	library(preprocessCore)

	library(ggplot2)
	library(reshape2)

	library(gridGraphics)
	library(grid)
	library(gridExtra)

	library(parallel)

	inputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'
	inputFiles <- list.files(inputDir, pattern = '.Rdata')

## Functions
	

	.parseOutputlist <- function(outputList, onlyPipettingIssues = FALSE, .debug = FALSE, TimesAboveDMSO = 2, ...) {
		if(.debug) {
			.outputList <- 'outputListP21.Rdata'
			load(paste0(inputDir, .outputList))
			TimesAboveDMSO <- 2
		}

		# Internal functions
		.parser <- function(dt, catProgress = TRUE, ...) {
			# Treatment names
			if(catProgress) cat('Parseing treatments...\n')
			dt$treatment <- gsub('_', '', toupper(gsub(' ', '', dt$treatment)))
			if(catProgress) cat('1 of 5 steps done.\n')

			# Standardize plateID in sumData
			if(catProgress) cat('Parseing plateIDs...\n')
			dt$plateID <- toupper(gsub('_wells', '', dt$plateID))
			if(catProgress) cat('2 of 5 steps done.\n')

			# Create uniqueID in sumData
			if(catProgress) cat('Parseing uniqueIDs...\n')
			dt$uniqueID <- paste(dt$treatment, dt$timeID, dt$dose_uM, dt$plateID, sep = '_')
			if(catProgress) cat('3 of 5 steps done.\n')

			# Split uniqueID and create CMAX columns
			if(catProgress) cat('Parseing CMAX column...\n')
			uniqueID_splitted <- do.call(rbind, strsplit(dt$uniqueID, '_'))
			dt$CMAX <- uniqueID_splitted[ ,grep('CMAX', toupper(uniqueID_splitted)[1, ])]
			dt$CMAX_numeric <- as.numeric(gsub('CMAX', '', dt$CMAX))
			if(catProgress) cat('4 of 5 steps done.\n')

			# Create columns to match against
			if(catProgress) cat('Parseing matching columns...\n')
			dt$matchDose <- paste(dt$treatment, dt$dose_uM, dt$CMAX, sep = '_')
			dt$matchDoseTime <- paste(dt$treatment, dt$dose_uM, dt$CMAX, dt$timeID,sep = '_')
			if(catProgress) cat('5 of 5 steps done.\n\n\n')

			# return
			return(dt)
		}

		# Assign names to replicates
		names(outputList) <- c('replicate_1', 'replicate_2', 'replicate_3')

		# Pipetting issues (replID > 3)
		pipIssues <- outputList$replicate_3$metaCSVData[which(outputList$replicate_3$metaCSVData$replID > 3), ]
		if(nrow(pipIssues) > 0) pipIssues <- .parser(pipIssues)
		if(onlyPipettingIssues) return(pipIssues)

		# outputList parse
		cat('Parseing myDTs.\n')
		parsedDTs <- lapply(outputList, function(x) .parser(x$myDT))	# CHOP myDT too big; mclapply() error: long vectors not supported yet

		cat('Parseing sumDatas.\n')
		parsedSDs <- mclapply(outputList, function(x) .parser(x$sumData), mc.cores = 6)

		# Add fraction >TimesAboveDMSO & Cellcount per timepoint * treatment combination
		for(rep in names(parsedSDs)) {
			if(.debug) .parsedDTs <- parsedDTs
			
			cat('Calculating fraction of cells above background for', rep, '...\n')
			plateIDs <- unique(parsedSDs[[rep]]$plateID)

			if(any(grepl('Cytoplasm_Intensity_IntegratedIntensity_image_GFP', colnames(parsedDTs[[rep]])))) {
				featureToCheckAboveDMSO <- 'Cytoplasm_Intensity_IntegratedIntensity_image_GFP'
			}	

			if(any(grepl('Nuclei_Intensity_IntegratedIntensity_image_GFP', colnames(parsedDTs[[rep]])))) {
				featureToCheckAboveDMSO <- 'Nuclei_Intensity_IntegratedIntensity_image_GFP'
			}	

			# Get DMSO median & mean for featureToCheckAboveDMSO
			DMSOstats <- do.call(rbind, mclapply(plateIDs, function(.pID) {
				.pIDind  <- which(parsedDTs[[rep]]$plateID == .pID)
				.dmsoInd <- which(parsedDTs[[rep]]$treatment == 'DMSO')

				medianDMSO <- median(parsedDTs[[rep]][intersect(.pIDind, .dmsoInd), ][[featureToCheckAboveDMSO]], na.rm = TRUE)
				meanDMSO  <- mean(parsedDTs[[rep]][intersect(.pIDind, .dmsoInd), ][[featureToCheckAboveDMSO]], na.rm = TRUE)
				
				c(medianDMSO, meanDMSO)
			}, mc.cores = 6))

				colnames(DMSOstats) <- c('median', 'mean')
				rownames(DMSOstats) <- plateIDs

			parsedDTs[[rep]]$plateDMSOmedian <- DMSOstats[parsedDTs[[rep]]$plateID, 'median']
			parsedDTs[[rep]]$plateDMSOmean <- DMSOstats[parsedDTs[[rep]]$plateID, 'mean']

			parsedDTs[[rep]]$abovePlateDMSOmedian <- parsedDTs[[rep]][[featureToCheckAboveDMSO]] > (TimesAboveDMSO * parsedDTs[[rep]]$plateDMSOmedian)
			parsedDTs[[rep]]$abovePlateDMSOmean <- parsedDTs[[rep]][[featureToCheckAboveDMSO]] > (TimesAboveDMSO * parsedDTs[[rep]]$plateDMSOmean)

			# Add fraction above DMSO to parsedSDs
			treatments <- parsedSDs[[rep]]$uniqueID
			
			fractionAboveMedian <- unlist(mclapply(treatments, function(.treatment) {
				isAboveMedian <- parsedDTs[[rep]]$abovePlateDMSOmedian[parsedDTs[[rep]]$uniqueID == .treatment]
				length(which(isAboveMedian)) / length(isAboveMedian)
			}, mc.cores = 6))

			fractionAboveMean <- unlist(mclapply(treatments, function(.treatment) {
				isAboveMean <- parsedDTs[[rep]]$abovePlateDMSOmean[parsedDTs[[rep]]$uniqueID == .treatment]
				length(which(isAboveMean)) / length(isAboveMean)
			}, mc.cores = 6))

			parsedSDs[[rep]]$fractionAboveDMSOmedian <- fractionAboveMedian
			parsedSDs[[rep]]$fractionAboveDMSOmean   <- fractionAboveMean	

			cat('\n')	

			## Number of objects
			numberOfObjects <- sapply(parsedSDs[[rep]]$matchDoseTime, function(x) sum(parsedDTs[[rep]]$matchDoseTime == x))
			parsedSDs[[rep]]$ObjectsInImages <- numberOfObjects

			# PI
			parsedSDs[[rep]]$fractionPI <- sapply(parsedSDs[[rep]]$matchDoseTime, function(x) {
				PI <- parsedDTs[[rep]]$PI_masked_primaryID_AreaShape_Area[which(parsedDTs[[rep]]$matchDoseTime == x)]
				sum(!(is.na(PI)))/length(PI)
			} )

			# AnV
			parsedSDs[[rep]]$fractionAnV <- sapply(parsedSDs[[rep]]$matchDoseTime, function(x) {
				AnV <- parsedDTs[[rep]]$AnV_masked_primaryID_AreaShape_Area[which(parsedDTs[[rep]]$matchDoseTime == x)]
				sum(!(is.na(AnV)))/length(AnV)
			} )
		}

		# Return
		rslt <- list(myDTs = parsedDTs,
					 summaryDatas = parsedSDs,
					 pipettingIssues = pipIssues)
		
		return(rslt)
	}

	plotDR <- function(reporter, compound, variableToPlot = c('fractionAboveDMSOmean', 'fractionAnV', 'fractionPI'), yLab = '', yLim = c(0, 1)) {
		plotData <- melt(do.call(rbind, lapply(reporter$summaryDatas, function(x) as.data.frame(x[x$treatment == compound, ])))[, c('treatment', 'timeID', 'CMAX_numeric', variableToPlot)], 
						 id.vars = c('timeID', 'CMAX_numeric', 'treatment'))
		plotData$CMAX_numeric <- as.factor(plotData$CMAX_numeric)

		P <- ggplot(data = plotData, aes(x = CMAX_numeric, y = value, group = timeID, colour = timeID)) +
				geom_point() + geom_smooth(se = FALSE) +
				ylab(yLab) + xlab('CMAX') + labs(title = compound) 

			if(!all(is.na(yLim))) P <- P + ylim(yLim)

		return(P)
	}



## Extra plots & compounds (i.e. V2): above 2x DMSO
	outputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/ResultsAnkeV2/'
	compOI <- c('ALLYLALCOHOL', 'RIFAMPICIN', 'DICLOFENAC', 'VALPROICACID', 
				'ACETAMINOPHEN', 'PIOGLITAZONE', 'AZATHIOPRINE', 'NITROFURANTOIN', 
				'CIPROFLOXACIN', 'TOLBUTAMIDE', 'CLOFIBRATE', 'FENOFIBRATE', 
				'PHENYTOIN')

	newCompounds <- c(compOI, 'EDROPHONIUM', 'CLOTRIMAZOLE', 'MECLIZINE', 'PRIMIDONE')	

	# BIP
		load(paste0(inputDir, 'outputListBIP.Rdata'))
		BIP <- .parseOutputlist(outputList)
		
		pdf(paste0(outputDir, 'BIP_FractionGFP_above2xDMSO_.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(BIP, x, variableToPlot = 'fractionAboveDMSOmean'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'BIP_FractionPI_nonNA.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(BIP, x, variableToPlot = 'fractionPI'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'BIP_FractionAnV_nonNA.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(BIP, x, variableToPlot = 'fractionAnV'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'BIP_ImageObjects.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(BIP, x, variableToPlot = 'ObjectsInImages', yLim = NA) )
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'BIP_meanImageGFPintensity.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(BIP, x, variableToPlot = 'Cytoplasm_Intensity_MeanIntensity_image_GFP'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

	# CHOP
		load(paste0(inputDir, 'outputListCHOP.Rdata'))
		CHOP <- .parseOutputlist(outputList)
		
		pdf(paste0(outputDir, 'CHOP_FractionGFP_above2xDMSO_.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(CHOP, x, variableToPlot = 'fractionAboveDMSOmean'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'CHOP_FractionPI_nonNA.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(CHOP, x, variableToPlot = 'fractionPI'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'CHOP_FractionAnV_nonNA.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(CHOP, x, variableToPlot = 'fractionAnV'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'CHOP_ImageObjects.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(CHOP, x, variableToPlot = 'ObjectsInImages', yLim = NA) )
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'CHOP_meanImageGFPintensity.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(CHOP, x, variableToPlot = 'Nuclei_Intensity_MeanIntensity_image_GFP'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

	# ICAM1
		load(paste0(inputDir, 'outputListICAM1.Rdata'))
		ICAM1 <- .parseOutputlist(outputList)
		
		pdf(paste0(outputDir, 'ICAM1_FractionGFP_above2xDMSO_.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(ICAM1, x, variableToPlot = 'fractionAboveDMSOmean'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'ICAM1_FractionPI_nonNA.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(ICAM1, x, variableToPlot = 'fractionPI'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'ICAM1_FractionAnV_nonNA.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(ICAM1, x, variableToPlot = 'fractionAnV'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'ICAM1_ImageObjects.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(ICAM1, x, variableToPlot = 'ObjectsInImages', yLim = NA) )
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'ICAM1_meanImageGFPintensity.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(ICAM1, x, variableToPlot = 'Cytoplasm_Intensity_MeanIntensity_image_GFP'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

	# P21
		load(paste0(inputDir, 'outputListP21.Rdata'))
		P21 <- .parseOutputlist(outputList)
		
		pdf(paste0(outputDir, 'P21_FractionGFP_above2xDMSO_.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(P21, x, variableToPlot = 'fractionAboveDMSOmean'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'P21_FractionPI_nonNA.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(P21, x, variableToPlot = 'fractionPI'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'P21_FractionAnV_nonNA.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(P21, x, variableToPlot = 'fractionAnV'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'P21_ImageObjects.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(P21, x, variableToPlot = 'ObjectsInImages', yLim = NA) )
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()

		pdf(paste0(outputDir, 'P21_meanImageGFPintensity.pdf'), width = 6, height = 16)
			plots.tmp <- lapply(newCompounds, function(x) plotDR(P21, x, variableToPlot = 'Nuclei_Intensity_MeanIntensity_image_GFP'))
			grid.arrange(grobs = plots.tmp, ncol = 2, clip = TRUE)
		dev.off()



	

######## LEGACY #########
	## Create single summaryData data.frame() with all reporters
		# BIP
		load(paste0(inputDir, 'outputListBIP.Rdata'))
			if(FALSE) .tmp <- outputList
		BIP <- .parseOutputlist(outputList)

		# CHOP
		load(paste0(inputDir, 'outputListCHOP.Rdata'))
		CHOP <- .parseOutputlist(outputList)

		# ICAM1
		load(paste0(inputDir, 'outputListICAM1.Rdata'))
		ICAM1 <- .parseOutputlist(outputList)

		# P21
		load(paste0(inputDir, 'outputListP21.Rdata'))
		P21 <- .parseOutputlist(outputList)

		# HSPA1B
		load(paste0(inputDir, 'outputListHSPA1B.Rdata'))
		HSPA1B <- .parseOutputlist(outputList)

	## Data & plotting for Anke
		outputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/ResultsAnke/'
		compOI <- c('ALLYLALCOHOL', 'RIFAMPICIN', 'DICLOFENAC', 'VALPROICACID', 
					'ACETAMINOPHEN', 'PIOGLITAZONE', 'AZATHIOPRINE', 'NITROFURANTOIN', 
					'CIPROFLOXACIN', 'TOLBUTAMIDE', 'CLOFIBRATE', 'FENOFIBRATE', 
					'PHENYTOIN')

		newCompounds <- c(compOI, 'EDROPHONIUM', 'CLOTRIMAZOLE', 'MECLIZINE', 'PRIMIDONE')
					

		# BIP
			load(paste0(inputDir, inputFiles[1]))
			BIP <- parseOutputlist(outputList)
			compOI.BIP <- lapply(compOI, function(x) plotDR(BIP, x))

				pdf(paste0(outputDir, 'BIP.pdf'), width = 6, height = 16)
					grid.arrange(grobs = compOI.BIP, ncol = 2, clip = TRUE)
				dev.off()

		# CHOP
			load(paste0(inputDir, inputFiles[2]))
			CHOP <- parseOutputlist(outputList)
			compOI.CHOP <- lapply(compOI, function(x) plotDR(CHOP, x))

				pdf(paste0(outputDir, 'CHOP.pdf'), width = 6, height = 16)
					grid.arrange(grobs = compOI.CHOP, ncol = 2, clip = TRUE)
				dev.off()

		# ICAM1
			load(paste0(inputDir, inputFiles[4]))
			ICAM1 <- parseOutputlist(outputList)
			compOI.ICAM1 <- lapply(compOI, function(x) plotDR(ICAM1, x))

				pdf(paste0(outputDir, 'ICAM1.pdf'), width = 6, height = 16)
					grid.arrange(grobs = compOI.ICAM1, ncol = 2, clip = TRUE)
				dev.off()

		# P21
			load(paste0(inputDir, inputFiles[5]))
			P21 <- parseOutputlist(outputList)
			compOI.P21 <- lapply(compOI, function(x) plotDR(P21, x))

				pdf(paste0(outputDir, 'P21.pdf'), width = 6, height = 16)
					grid.arrange(grobs = compOI.P21, ncol = 2, clip = TRUE)
				dev.off()

		# HSPA1B
			load(paste0(inputDir, inputFiles[3]))
			HSPA1B <- parseOutputlist(outputList)
			compOI.HSPA1B <- lapply(compOI, function(x) plotDR(HSPA1B, x))

				pdf(paste0(outputDir, 'HSPA1B.pdf'), width = 6, height = 16)
					grid.arrange(grobs = compOI.HSPA1B, ncol = 2, clip = TRUE)
				dev.off()

		# SRXN1
			load(paste0(inputDir, inputFiles[6]))
			SRXN1 <- parseOutputlist(outputList)
			compOI.SRXN1 <- lapply(compOI, function(x) plotDR(SRXN1, x))

				pdf(paste0(outputDir, 'SRXN1.pdf'), width = 6, height = 16)
					grid.arrange(grobs = compOI.SRXN1, ncol = 2, clip = TRUE)
				dev.off()
















