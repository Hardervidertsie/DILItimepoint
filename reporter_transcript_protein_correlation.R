## Settings
	options(stringsAsFactors = FALSE)
	
	library(ggplot2)

	library(data.table)
	library(preprocessCore)
	library(reshape2)

	library(stringr)
	library(parallel)
	library(reshape2)

	library(DESeq2)

	
## GFP Data
	inputDir_GFP <- '/Users/Wouter/stack/ownCloud/TOX/ExperimentData/DILIscreen/RDataFiles_EndPointDILI/'

	load(paste0(inputDir_GFP, 'parsedBIP.RData'))
	load(paste0(inputDir_GFP, 'parsedSRXN1.RData'))
	load(paste0(inputDir_GFP, 'parsedP21.RData'))
	load(paste0(inputDir_GFP, 'parsedCHOP.RData'))
	load(paste0(inputDir_GFP, 'parsedICAM1.RData'))
	load(paste0(inputDir_GFP, 'parsedHSPA1B.RData'))
	load(paste0(inputDir_GFP, 'parsedBTG2.RData'))
	load(paste0(inputDir_GFP, 'parsedHMOX1.RData'))

  # Parseing
	HSPA1B <- HSPA1B[1:2]	# Replicate three failed


## BioSpyder data
	load('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/PHH50/RESULTS/HepG2rp.RData')

	probesOI <- c(BIP = 'HSPA5_14023', SRXN1 = 'SRXN1_20094', P21 = 'CDKN1A_1219', CHOP = 'DDIT3_16736', ICAM1 = 'ICAM1_3187', HSPA1B = 'HSPA1B_3136', BTG2 = 'BTG2_13191', HMOX1 = 'HMOX1_3041')
	probesOI <- probesOI[-which(names(probesOI) == 'HSPA1B' | names(probesOI) == 'ICAM1')]	# No BioSpyder data

## Functions
	getData <- function(reporter, timepoint = 24, compound = c('TUN', 'CPT', 'TNF', 'DEM'), .debug = FALSE) {
		if(.debug) {
			reporter <- 'SRXN1'
			timepoint <- 24
			compound <- 'DEM'
		}

		compound_HepG2 <- c(TUN = 'TUNICAMYCIN', CPT = 'CISPLATIN', TNF = 'TNF-A', DEM = 'DEM')[compound]

		if(!exists(reporter)) stop(paste0('Can not find ', reporter, ' data.\n' ))
		if(!exists('meta_HepG2rp')) stop('Can not find BioSpyder meta data.\n')


		data_BS <- meta_HepG2rp[which(meta_HepG2rp$CELL_ID == reporter & 
									  meta_HepG2rp$TIMEPOINT == timepoint & 
									  meta_HepG2rp$TREATMENT == compound), ]
		data_BS <- data_BS[order(data_BS$CONCENTRATION), ]
		data_BS$cpm <- cpm_HepG2rp[probesOI[reporter], data_BS$SAMPLE_NAME]
		
		data_HepG2 <- do.call(rbind, lapply(get(reporter), function(.replicate) {
			.sumData <- .replicate$sumData
			.sumData <- .sumData[which(.sumData$treatment == compound_HepG2[compound] & .sumData$timeID == timepoint), ]

			.sumData
		} ))

		data_HepG2 <- data_HepG2[order(data_HepG2$dose_uM), ]
		if(!is.numeric(data_HepG2$dose_uM)) data_HepG2$dose_uM <- as.numeric(as.character(data_HepG2$dose_uM))

		if(cor(data_BS$CONCENTRATION, data_HepG2$dose_uM, method = 'spearman')) {  
			return(list(BioSpyder = data_BS, HepG2 = as.data.frame(data_HepG2)))
		} else {
			return(NA)
		}
	}

	checkCor <- function(BS_HepG2_list, returnData = TRUE, HepG2_fractionGFP = c(1, 2, 3, 5), .debug = FALSE) {
		if(.debug) {
			BS_HepG2_list <- getData(reporter = 'BTG2', compound = 'CPT')
			returnData <- TRUE
			HepG2_fractionGFP <- c(1, 2, 3, 5)
		}

		pearson <- sapply(HepG2_fractionGFP, function(.cutoff) {
			cor(BS_HepG2_list$BioSpyder$cpm, BS_HepG2_list$HepG2[, paste0('Fraction_Above_', .cutoff, '_meanIntregated_plateDMSO')], method = 'pearson')
		})
		names(pearson) <- c('cutoff_1', 'cutoff_2', 'cutoff_3', 'cutoff_5')

		spearman <- sapply(HepG2_fractionGFP, function(.cutoff) {
			cor(BS_HepG2_list$BioSpyder$cpm, BS_HepG2_list$HepG2[, paste0('Fraction_Above_', .cutoff, '_meanIntregated_plateDMSO')], method = 'spearman')
		})
		names(pearson) <- c('cutoff_1', 'cutoff_2', 'cutoff_3', 'cutoff_5')

		if(returnData) {
			return(list(pearson = pearson, spearman = spearman, data = BS_HepG2_list))
		} else {
			return(list(pearson = pearson, spearman = spearman))
		}
	}


## Magic
	cors2check <- data.frame(compound = sapply(names(probesOI), function(x) names(which.max(table(meta_HepG2rp[which(meta_HepG2rp$CELL_ID == x), 'TREATMENT'])))),
							 probe    = probesOI,
							 reporter = names(probesOI))

	rslt <- lapply(1:nrow(cors2check), function(i) {
		cat(i, '\n')
		.d <- getData(reporter = cors2check$reporter[i],
					  compound = cors2check$compound[i])

		checkCor(.d)
	})





















