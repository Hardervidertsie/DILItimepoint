## Acquire SMILES IDs

## Options
  options(stringsAsFactors = FALSE)


## Data
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
	
	DILI[DILI$DILIConcern == 'N/A', c('doseDILI', 'binaryDILI')] <- NA
	
	rownames(DILI) <- DILI$treatment

  smiles <- read.delim('/Users/Wouter/Desktop/smiles.txt')
  smiles$treatment <- toupper(smiles$LabelCompoundName)

## Magic
  compound <- data.frame(compound = DILI$treatment)
    compound$SMILES <- sapply(compound$compound, function(x) smiles$SMILES[which(smiles$treatment == x)] )
    compound$SMILES[which(sapply(compound$SMILES, length) == 0)] <- NA
    compound$SMILES <- unlist(compound$SMILES)

## Save
  write.table(compound, file = '/Users/Wouter/stack/ownCloud/TOX/ExperimentData/DILIscreen/SMILES.txt', quote = FALSE, sep = '\t', row.names = FALSE)
