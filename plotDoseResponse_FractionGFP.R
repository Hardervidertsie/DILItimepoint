## Settings
	options(stringsAsFactors = FALSE)
	
	library(ggplot2)

	library(data.table)
	library(preprocessCore)
	library(reshape2)

	inputDir <- '/Users/Wouter/stack/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'


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
	plotDoseResponse <- function(cmpnd, feature = c('FractionOfCells','Intensity'), dmsoCutoff = c(1, 2, 3, 5), dmsoFeature = c('median', 'mean'), reporters = c('BIP', 'SRXN1', 'P21', 'CHOP', 'ICAM1', 'HSPA1B', 'BTG2', 'HMOX1'), .debug = FALSE) {
		if(.debug) {
			cmpnd <- 'ROTENONE'
			dmsoCutoff <- 2
			dmsoFeature <- 'median'
		}

		tmp <- do.call(rbind, lapply(reporters, function(.reporter) do.call(rbind, lapply(get(.reporter), function(x) { 
			.tmp <- data.frame(x$sumData[which(x$sumData$treatment == cmpnd), ]) 
			 # rownames(.tmp) <- .tmp$matchDoseTime
			  .tmp$reporter <- .reporter

			if(feature == 'FractionOfCells') .feature <- grep(paste0('Fraction_Above_', dmsoCutoff, '_', dmsoFeature), colnames(.tmp), value = TRUE)
			if(feature == 'Intensity')       .feature <- grep('image_GFP', colnames(.tmp), value = TRUE)

			.tmp <- .tmp[, c('reporter', 'treatment', 'timeID', 'repID', 'CMAX_numeric', .feature)]
			  colnames(.tmp)[6] <- 'GFP'

			.tmp
		}))))

		tmp$CMAX_numeric <- factor(tmp$CMAX_numeric, ordered = TRUE)


		P <- ggplot(tmp, aes(CMAX_numeric, GFP, colour = reporter)) +
			geom_boxplot() +
			# ylim(0, 1) +
			theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8), legend.position = 'none') +
			xlab('CMAX') + 
			# ylab('Fraction GFP+ cells') + 
			ggtitle(cmpnd)

		if(feature == 'FractionOfCells') P <- P + ylim(0, 1) + ylab('Fraction GFP+ cells') + facet_grid(reporter ~ timeID) 
		if(feature == 'Intensity')       P <- P + ylab('GFP Intensity') + facet_grid(reporter ~ timeID, scales = 'free_y')
			

		return(P)	
	}

## Wanda - 03 Apr 2017
  # Compounds
    wanda <- c('ROTENONE', 'ANTIMYCINA', 'OLIGOMYCINA', 'FCCP', 'CCCP', 'CYCLOSPORINA',		# Compounds of interest
    		   'TUNICAMYCIN', 'THAPSIGARGIN', 'ETOPOSIDE', 'CISPLATIN', 'DEM', 'CDDO-ME')	# Controls of interest

  # Plots
    pdf('/Users/Wouter/Desktop/_Wanda07Apr2017_FractionGFP.pdf', height = 8, width = 6)
      for(i in wanda) print(plotDoseResponse(cmpnd = i, feature = 'FractionOfCells', dmsoCutoff = 2, dmsoFeature = 'median'))
    dev.off()

    pdf('/Users/Wouter/Desktop/_Wanda07Apr2017_IntensityGFP.pdf', height = 8, width = 6)
      for(i in wanda) print(plotDoseResponse(cmpnd = i, feature = 'Intensity', dmsoCutoff = 2, dmsoFeature = 'median'))
    dev.off()

## Anke - 17 May 2017
	anke <- c('VALPROICACID', 'DICLOFENAC', 'ACETAMINOPHEN', 'CIPROFLOXACIN', 'NITROFURANTOIN', 'TOLCAPONE', 'AZATHIOPRINE', 'TROGLITAZONE',
			  'NEFAZODONE', 'KETOCONAZOLE', 'OMEPRAZOLE', 'PHENYTOIN', 'AMIODARONE', 'VERAPAMIL', 'BUSPIRONE', 'FAMOTIDINE', 'DMSO', 'TUNICAMYCIN')

	pdf('/Users/Wouter/Desktop/Anke_17May2017_FractionGFP.pdf', height = 8, width = 6)
      for(i in anke) print(plotDoseResponse(cmpnd = i, feature = 'FractionOfCells', dmsoCutoff = 2, dmsoFeature = 'median'))
    dev.off()

## Nanette Valproaic Acid - Mon 13 Feb 2017
  # Heatmap GFP+ fraction cells	
	pdf('/Users/Wouter/Desktop/VPA_timepointScreen.pdf', height = 8, width = 6)
		plotDoseResponse(cmpnd = 'VALPROICACID', dmsoCutoff = 2, dmsoFeature = 'median')
	dev.off()

## Suus - Mon 24 Apr 2017
	suus <- c('VALPROICACID', 'COLCHICINE', 'GRISEOFULVIN', 'DICLOFENAC', 'CARBAMAZAPINE')

	pdf('/Users/Wouter/Desktop/TMP_suus_mean.pdf', height = 8, width = 6)
      for(i in suus) print(plotDoseResponse(cmpnd = i, feature = 'FractionOfCells', dmsoCutoff = 2, dmsoFeature = 'mean'))
    dev.off()

    getData <- function(cmpnd, reporters = c('BIP', 'SRXN1', 'P21', 'CHOP', 'ICAM1', 'HSPA1B', 'BTG2', 'HMOX1')) {
    	if(FALSE) {
    		cmpnd <- 'VALPROICACID'
    		reporters <- c('BIP', 'SRXN1', 'P21', 'CHOP', 'ICAM1', 'HSPA1B', 'BTG2', 'HMOX1')
    	}

    	tmp <- lapply(reporters, function(.reporter) do.call(rbind, lapply(get(.reporter), function(x) { 
			.tmp <- data.frame(x$sumData[which(x$sumData$treatment == cmpnd), ]) 
			.tmp$reporter <- .reporter

			.tmp
		})))

		names(tmp) <- sapply(tmp, function(x) unique(x$reporter))

		return(tmp)
    }

    suusSumData <- lapply(suus, getData)

    save(suusSumData, file = '/Users/Wouter/stack/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/sumData_Suus_24Apr2017.RData')

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
  	
## Check PI staining correlation across reporters
	reporters <- c('BIP', 'SRXN1', 'P21', 'CHOP', 'ICAM1', 'HSPA1B', 'BTG2', 'HMOX1')
	sumData_PI <- lapply(reporters, function(.reporter) {
		lapply(get(.reporter), function(x) {
			tmp <- as.data.frame(x$sumData)[, c('treatment', 'CMAX', 'timeID', 'PI_fractionNonNA')]
			tmp$ID <- apply(tmp[, 1:3], 1, function(z) paste(z, collapse = '_'))

			tmp
		} )
	})

	toCheck <- names(which(table(unlist(lapply(sumData_PI, function(x) lapply(x, function(y) y$ID)))) == 23))	# 23 is the highest - has most measurements accross all reporters

	dat <- lapply(sumData_PI, function(x) do.call(cbind, lapply(x, function(y) {
		tmp <- y[which(y$ID %in% toCheck), ] 
		tmp[order(tmp$ID), ]
	}) ))

	dat <- lapply(dat, function(x) if(ncol(x) == 15) x[, c(1, 2, 3, 5, 4, 9, 14)] else x[, c(1, 2, 3, 5, 4, 9)] )
	dat <- lapply(dat, function(x) {
		x$meanReplicateFractionPI <- apply(x[, 5:ncol(x)], 1, function(y) mean(as.numeric(y)))
		x
	})
	
	names(dat) <- reporters

	dat <- do.call(cbind, lapply(dat, function(x) x[, c('ID', 'meanReplicateFractionPI')]))
	  rownames(dat) <- dat[, 1]

	datMeanPI <- dat[, grep('meanRepl', colnames(dat))]
	  colnames(datMeanPI) <- do.call(rbind, strsplit(colnames(datMeanPI), '\\.'))[, 1]

	library(pheatmap)
	library(RColorBrewer)

	pheatmap(datMeanPI[grep('_72', rownames(datMeanPI)), ], 
			 color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100),
			 clustering_distance_cols = 'correlation',
			 clustering_rows = FALSE)

	pheatmap(cor(datMeanPI[grep('_72', rownames(datMeanPI)), ]), 
			 colorRampPalette((brewer.pal(n = 7, name = "RdYlBu")))(100),
			 clustering_rows = FALSE,
			 clustering_rows = FALSE,
			 breaks = c(seq(0, 0.8, 0.016), seq(0.804, 1, 0.004)))

## PI/AnV dose-response like curves for Isoude - 28 Jun 2017
	reporters <- c('BIP', 'SRXN1', 'P21', 'CHOP', 'ICAM1', 'HSPA1B', 'BTG2', 'HMOX1')
	sumData_PI <- lapply(reporters, function(.reporter) {
		lapply(get(.reporter), function(x) {
			tmp <- as.data.frame(x$sumData)[, c('treatment', 'CMAX', 'timeID', 'PI_fractionNonNA')]
			tmp$ID <- apply(tmp[, 1:3], 1, function(z) paste(z, collapse = '_'))

			tmp
		} )
	})

	plotDoseResponse_PI <- function(cmpnd, reporters = c('BIP', 'SRXN1', 'P21', 'CHOP', 'ICAM1', 'HSPA1B', 'BTG2', 'HMOX1'), .debug = FALSE) {
		if(.debug) {
			cmpnd <- 'BUTHIONINESULFOXIMINE'
			reporters <- c('BIP', 'SRXN1', 'P21')
		}

		tmp <- do.call(rbind, lapply(reporters, function(.reporter) do.call(rbind, lapply(get(.reporter), function(x) { 
			.tmp <- data.frame(x$sumData[which(x$sumData$treatment == cmpnd), ]) 
			 # rownames(.tmp) <- .tmp$matchDoseTime
			  .tmp$reporter <- .reporter

			.feature <- 'PI_fractionNonNA'
			

			.tmp <- .tmp[, c('reporter', 'treatment', 'timeID', 'repID', 'CMAX_numeric', .feature)]
			  colnames(.tmp)[6] <- 'PI'

			.tmp
		}))))

		tmp$CMAX_numeric <- factor(tmp$CMAX_numeric, ordered = TRUE)
		tmp$PI <- as.numeric(tmp$PI)


		P <- ggplot(tmp, aes(CMAX_numeric, PI, colour = reporter)) +
			geom_boxplot() +
			# ylim(0, 1) +
			theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8), legend.position = 'none') +
			xlab('CMAX') + 
			# ylab('Fraction GFP+ cells') + 
			ggtitle(cmpnd)

		 P <- P + ylim(0, 1) + ylab('Fraction PI+ cells') + facet_grid(reporter ~ timeID) 	

		return(P)	
	}

	compounds <- sort(unique(BIP[[1]]$sumData$treatment))
	compounds <- compounds[-which(compounds == 'STAUSPORIN' | compounds == 'SODIUM-ARSENITE' | compounds == 'BUTHIONINESULFOXIMINE' | compounds == 'COCL2C' | compounds == 'FLURBIPROFEN')]

	pdf('/Users/Wouter/Desktop/doseResponse_PI.pdf', width = 9, height = 9)
		mclapply(compounds, function(.cmpnd) {
			cat(.cmpnd, '\n')
			plotDoseResponse_PI(cmpnd = .cmpnd,)
		}, mc.cores = 6)
	dev.off()

	plotDoseResponse_AnV <- function(cmpnd, reporters = c('BIP', 'SRXN1', 'P21', 'CHOP', 'ICAM1', 'HSPA1B', 'BTG2', 'HMOX1'), .debug = FALSE) {
		if(.debug) {
			cmpnd <- 'BUTHIONINESULFOXIMINE'
			reporters <- c('BIP', 'SRXN1', 'P21')
		}

		tmp <- do.call(rbind, lapply(reporters, function(.reporter) do.call(rbind, lapply(get(.reporter), function(x) { 
			.tmp <- data.frame(x$sumData[which(x$sumData$treatment == cmpnd), ]) 
			 # rownames(.tmp) <- .tmp$matchDoseTime
			  .tmp$reporter <- .reporter

			.feature <- 'AnV_fractionNonNA'
			

			.tmp <- .tmp[, c('reporter', 'treatment', 'timeID', 'repID', 'CMAX_numeric', .feature)]
			  colnames(.tmp)[6] <- 'AnV'

			.tmp
		}))))

		tmp$CMAX_numeric <- factor(tmp$CMAX_numeric, ordered = TRUE)
		tmp$AnV <- as.numeric(tmp$AnV)


		P <- ggplot(tmp, aes(CMAX_numeric, AnV, colour = reporter)) +
			geom_boxplot() +
			# ylim(0, 1) +
			theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8), legend.position = 'none') +
			xlab('CMAX') + 
			# ylab('Fraction GFP+ cells') + 
			ggtitle(cmpnd)

		 P <- P + ylim(0, 1) + ylab('Fraction AnV+ cells') + facet_grid(reporter ~ timeID) 	

		return(P)	
	}

	compounds <- sort(unique(BIP[[1]]$sumData$treatment))
	compounds <- compounds[-which(compounds == 'STAUSPORIN' | compounds == 'SODIUM-ARSENITE' | compounds == 'BUTHIONINESULFOXIMINE' | compounds == 'COCL2C' | compounds == 'FLURBIPROFEN')]

	pdf('/Users/Wouter/Desktop/doseResponse_AnV.pdf', width = 9, height = 9)
		lapply(compounds, function(.cmpnd) {
			cat(.cmpnd, '\n')
			plotDoseResponse_AnV(cmpnd = .cmpnd,)
		})
	dev.off()





writeOut_sumData <- function(reporter, basename) {
	lapply(1:length(reporter), function(replicate) {
	  	write.table(reporter[[replicate]]$sumData, file = paste0(basename, '_replicate_', replicate, '.txt'), sep = '\t', quote = FALSE)
	})

}

writeOut_sumData(BIP, 'BIP')
writeOut_sumData(BTG2, 'BTG2')
writeOut_sumData(CHOP, 'CHOP')
writeOut_sumData(HMOX1, 'HMOX1')
writeOut_sumData(HSPA1B, 'HSPA1B')
writeOut_sumData(ICAM1, 'ICAM1')
writeOut_sumData(P21, 'P21')
writeOut_sumData(SRXN1, 'SRXN1')


















