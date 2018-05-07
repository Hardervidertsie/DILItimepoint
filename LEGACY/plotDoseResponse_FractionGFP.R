## Settings
	options(stringsAsFactors = FALSE)
	
	library(ggplot2)

	library(data.table)
	library(preprocessCore)
	library(reshape2)

	inputDir <- '/Users/Wouter/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'


## Data
	load(paste0(inputDir, 'parsedBIP.RData'))
	load(paste0(inputDir, 'parsedSRXN1.RData'))
	load(paste0(inputDir, 'parsedP21.RData'))
	load(paste0(inputDir, 'parsedCHOP.RData'))
	load(paste0(inputDir, 'parsedICAM1.RData'))
	load(paste0(inputDir, 'parsedHSPA1B.RData'))
	load(paste0(inputDir, 'parsedBTG2.RData'))
	load(paste0(inputDir, 'parsedHMOX1.RData'))

  # Parseing
	HSPA1B <- HSPA1B[1:2]	# Replicate three failed


## Functions
	plotDoseResponse <- function(cmpnd, dmsoCutoff = c(1, 2, 3, 5), dmsoFeature = c('median', 'mean'), reporters = c('BIP', 'SRXN1', 'P21', 'CHOP', 'ICAM1', 'HSPA1B', 'BTG2', 'HMOX1'), .debug = FALSE) {
		if(.debug) {
			cmpnd <- 'VALPROICACID'
			dmsoCutoff <- 2
			dmsoFeature <- 'median'
		}

		tmp <- do.call(rbind, lapply(reporters, function(.reporter) do.call(rbind, lapply(get(.reporter), function(x) { 
			.tmp <- data.frame(x$sumData[grep(cmpnd, x$sumData$treatment), ]) 
			 # rownames(.tmp) <- .tmp$matchDoseTime
			  .tmp$reporter <- .reporter

			.tmp <- .tmp[, c('reporter', 'treatment', 'timeID', 'repID', 'CMAX_numeric', grep(paste0('Fraction_Above_', dmsoCutoff, '_', dmsoFeature), colnames(.tmp), value = TRUE))]
			  colnames(.tmp)[6] <- 'GFP'

			.tmp
		}))))

		tmp$CMAX_numeric <- factor(tmp$CMAX_numeric, ordered = TRUE)


		P <- ggplot(tmp, aes(CMAX_numeric, GFP, colour = reporter)) +
			geom_boxplot() +
			# geom_smooth(method = 'loess') +
			facet_grid(reporter ~ timeID) +
			ylim(0, 1) +
			theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8), legend.position = 'none') +
			xlab('CMAX') + ylab('Fraction GFP+ cells') + ggtitle(cmpnd)

		return(P)	
	}


## Nanette Valproaic Acid - Mon 13 Feb 2017
  # Heatmap GFP+ fraction cells	
	pdf('/Users/Wouter/Desktop/VPA_timepointScreen.pdf', height = 8, width = 6)
		plotDoseResponse(cmpnd = 'VALPROICACID', dmsoCutoff = 2, dmsoFeature = 'median')
	dev.off()


## SOT 2017 poster
	pdf('/Users/Wouter/Desktop/SOT2017_doseResponse_timepointScreen.pdf', height = 6, width = 5)
		plotDoseResponse(cmpnd = 'ACETAMINOPHEN', dmsoCutoff = 2, dmsoFeature = 'median')
		plotDoseResponse(cmpnd = 'BENZBROMARONE', dmsoCutoff = 2, dmsoFeature = 'median')
		plotDoseResponse(cmpnd = 'CISPLATIN', dmsoCutoff = 2, dmsoFeature = 'median')
		plotDoseResponse(cmpnd = 'CLOZAPINE', dmsoCutoff = 2, dmsoFeature = 'median')
		plotDoseResponse(cmpnd = 'FLUOXETINE', dmsoCutoff = 2, dmsoFeature = 'median')
		plotDoseResponse(cmpnd = 'ACETAMINOPHEN', dmsoCutoff = 2, dmsoFeature = 'median')
		plotDoseResponse(cmpnd = 'METHOTREXATE', dmsoCutoff = 2, dmsoFeature = 'median')
		plotDoseResponse(cmpnd = 'XIMELAGATRAN', dmsoCutoff = 2, dmsoFeature = 'median')
		plotDoseResponse(cmpnd = 'OMEPRAZOLE', dmsoCutoff = 2, dmsoFeature = 'median')
	dev.off()
  	

























