

# function input: output from `_0_parse_outputlists.R'`

require(dplyr)
require(data.table)


##==## create_summaries
# function with parsed lists containing 3 reps as input
# select relevant columns
# calculate means and fractions
# create 1 table from the 3 replicates
# store table


# load data, check number of unique dose_uM and CMAX value combinations, save if not 6 combinations per compound.
# select relevant columns
# output result
clean_load_sel_count_dose_cmax <- function( parsedlistpath, name, cyto ) {
  load(paste0(parsedlistpath, name, '/parsed', name, '.Rdata'))
  
  myDT_list <- lapply( get(name), '[[', 'myDT' )
  
  count_concentrations <- lapply( myDT_list, function(x) {
    unique(x[, list(treatment, dose_uM, CMAX, replID)])
  }
  )
  
  count_concentrations_over6 <- lapply(count_concentrations,
                                       function(x) { x[ , .N, by = treatment][N!=6]})
  for(i in seq_along(count_concentrations_over6)){
    count_concentrations_over6[[i]]$repID <- i
  }
  count_concentrations <-  do.call( 'rbind', count_concentrations_over6)
  
  write.table(count_concentrations[ treatment %in% count_concentrations_over6$treatment, ],
              file = paste0('../generated/validation/more6conc_', name, '.txt'), sep = '\t', row.names = FALSE)

  # select columns, if a column not found an error will be returned by data.table
  
  myDT_list_sel <- lapply(myDT_list, function(x) {
    x[ , c('treatment', 'dose_uM', 'CMAX', 'cell_line', 'plateID', 'timeID', 'replID', 'imageCountParentObj', 'locationID', 
              'PI_masked_primaryID_AreaShape_Area', 'AnV_masked_primaryID_AreaShape_Area', 
              ifelse(TRUE, 'Cytoplasm_Intensity_MeanIntensity_image_GFP', 'Nuclei_Intensity_MeanIntensity_image_GFP' )
              ,
              ifelse(TRUE, 'Cytoplasm_Intensity_IntegratedIntensity_image_GFP', 'Nuclei_Intensity_IntegratedIntensity_image_GFP')
                )]
  } 
)
  
  return(myDT_list_sel)
    
}

##==## calculate GFP count values, base GFP intensity counts on log transformed single cell integrated intensity values
calc_GFP_Counts <- function(inputdata, ICAM1 = FALSE) {
  
 # get variable colnames 
  meanIntColname <- lapply(inputdata, function(x) {
    colnames(x)[ which(grepl( 'Intensity_MeanIntensity_image_GFP' , colnames(x) ))]
  }
  )
  
  integrIntColname <- lapply(inputdata, function(x) {
    colnames(x)[ which(grepl( 'Intensity_IntegratedIntensity_image_GFP' , colnames(x) ))]
  }
  )
  if(length(unique(meanIntColname)) != 1| length(unique(integrIntColname)) !=1  ){
    stop("replid colname error")
  }
  
  meanIntColname <- unlist(unique(meanIntColname))
  integrIntColname <- unlist(unique(integrIntColname))
  
  # calculate count threshold values based on DMSO values (or TNF for ICAM1)
  
  
  if(length(meanIntColname) != 0 & length(integrIntColname) == 0 ){
      dmso_platemeans <-lapply(inputdata, function(x) {
        x[ treatment %in% "DMSO",  mean(log2(get(meanIntColname) + 1), na.rm = TRUE) , by = plateID]
      }
      )
      dmso_platemeans <- lapply(dmso_platemeans, function(x) {setnames(x, "V1", "DMSO_meanint_plate_lognorm")})
 
  } else if(length(integrIntColname) != 0 ){
    dmso_platemeans<- dmso_platemeans <- lapply(inputdata, function(x) {
      x[ treatment %in% "DMSO" ,  mean(log2(get(integrIntColname) + 1), na.rm = TRUE) , by = plateID]
    }
    )
    dmso_platemeans <- lapply(dmso_platemeans, function(x) {setnames(x, "V1", "DMSO_intint_plate_lognorm")})
  } else {
    stop("No mean or integr column names match")
  }
   
  
  # join tables
  dmso_platemeans<-lapply(dmso_platemeans, function(x) setkey(x, "plateID")) 
  inputdata <- lapply(inputdata, function(x) setkey(x, "plateID"))
  
  inputdata <- mapply(function(x, y) {y[x]}, x = dmso_platemeans, y = inputdata, SIMPLIFY = FALSE)
  
  # count GFP cells above 1X 2X, 3X and 4X DMSO platelognorm
  lognorm <- function(x) log2(x + 1)
  if(!ICAM1){
  if( "DMSO_meanint_plate_lognorm" %in% colnames(inputdata[[1]])  & !"DMSO_intint_plate_lognorm" %in% colnames(inputdata[[1]])){
    inputdata <- lapply(inputdata, function(x) x[ , c("GFP_pos1m", "GFP_pos2m", "GFP_pos3m", "GFP_pos4m") :=
                                                           list(lognorm(get(meanIntColname)) > 1*DMSO_meanint_plate_lognorm,
                                                                lognorm(get(meanIntColname)) > 2*DMSO_meanint_plate_lognorm,
                                                                lognorm(get(meanIntColname)) > 3*DMSO_meanint_plate_lognorm,
                                                                lognorm(get(meanIntColname)) > 4*DMSO_meanint_plate_lognorm)]
    )
    
                                                  
  
  } else if("DMSO_intint_plate_lognorm" %in% colnames(inputdata[[1]])){
    inputdata <- lapply(inputdata, function(x) x[ , c("GFP_pos1i", "GFP_pos2i", "GFP_pos3i", "GFP_pos4i") := 
                                                    list(lognorm(get(integrIntColname)) > 1*DMSO_intint_plate_lognorm, 
                                                         lognorm(get(integrIntColname)) > 2*DMSO_intint_plate_lognorm,
                                                         lognorm(get(integrIntColname)) > 3*DMSO_intint_plate_lognorm,
                                                         lognorm(get(integrIntColname)) > 4*DMSO_intint_plate_lognorm)]
                        )
  } else{
    stop("No lognorm column exists")
  }
    
  } else{ #ICAM1
    
    if( "DMSO_meanint_plate_lognorm" %in% colnames(inputdata[[1]])  & !"DMSO_intint_plate_lognorm" %in% colnames(inputdata[[1]])){
      inputdata <- lapply(inputdata, function(x) x[ , c("GFP_pos1m", "GFP_pos2m", "GFP_pos3m", "GFP_pos4m",
                                                        "GFP_neg1m", "GFP_neg2m", "GFP_neg3m", "GFP_neg4m") :=
                                                      list(lognorm(get(meanIntColname)) > 1*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(meanIntColname)) > 2*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(meanIntColname)) > 3*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(meanIntColname)) > 4*DMSO_meanint_plate_lognorm,
                                                           
                                                           lognorm(get(meanIntColname)) < 1*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(meanIntColname)) < 1/2*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(meanIntColname)) < 1/3*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(meanIntColname)) < 1/4*DMSO_meanint_plate_lognorm)]
      )
      
       
      
    } else if("DMSO_intint_plate_lognorm" %in% colnames(inputdata[[1]])){
      inputdata <- lapply(inputdata, function(x) x[ , c("GFP_pos1i", "GFP_pos2i", "GFP_pos3i", "GFP_pos4i",
                                                        "GFP_neg1i", "GFP_neg2i", "GFP_neg3i", "GFP_neg4i") := 
                                                      list(lognorm(get(integrIntColname)) > 1*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(integrIntColname)) > 2*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(integrIntColname)) > 3*DMSO_meanint_plate_lognorm,
                                                           lognorm(get(integrIntColname)) > 4*DMSO_meanint_plate_lognorm,
                                                           
                                                           lognorm(get(integrIntColname)) < 1*DMSO_intint_plate_lognorm,
                                                           lognorm(get(integrIntColname)) < 1/2*DMSO_intint_plate_lognorm,
                                                           lognorm(get(integrIntColname)) < 1/3*DMSO_intint_plate_lognorm,
                                                           lognorm(get(integrIntColname)) < 1/4*DMSO_intint_plate_lognorm)]
      )
    } else{
      stop("No lognorm column exists")
    }
    
    
  } #//end ICAM1
  return(inputdata)
  }

# calculate cell relative cell counts compared to plate DMSO
norm_cell_counts <- function(inputdata) {
  
  dmso_platemeans <-lapply(inputdata, function(x) {
    x[ treatment %in% "DMSO",  mean(imageCountParentObj, na.rm = TRUE),  , by = plateID]
  }
  )
  
  dmso_platemeans <- lapply(dmso_platemeans, function(x) {setnames(x, "V1", "DMSO_mean_cell_counts")})

# join tables
dmso_platemeans<-lapply(dmso_platemeans, function(x) setkey(x, "plateID")) 
inputdata <- lapply(inputdata, function(x) setkey(x, "plateID"))

inputdata <- mapply(function(x, y) {y[x]}, x = dmso_platemeans, y = inputdata, SIMPLIFY = FALSE)
  
inputdata <- lapply(inputdata, function(x) {
  x[, norm_cell_counts := imageCountParentObj / DMSO_mean_cell_counts]
})
  
return(inputdata)
}

# calculate dead cells (more than 2 pixels evidence)
calc_cell_death <- function(inputdata){
 # change cell death column names

inputdata <-  lapply(inputdata, function(x) {
  setnames(x, "PI_masked_primaryID_AreaShape_Area", "PI_pos")
  setnames(x, "AnV_masked_primaryID_AreaShape_Area", "Annexin_pos")
})

   # replace NA by 0
inputdata <- lapply( inputdata, function(x) {
                       x[ , c("PI_pos", "Annexin_pos"):=
                         list(ifelse(is.na(PI_pos), 0, PI_pos ),
                         ifelse(is.na(Annexin_pos), 0, Annexin_pos))]
    }
    )
                         

inputdata <- lapply( inputdata, function(x) {
  x[ , c("PI_pos", "Annexin_pos"):=
       list(ifelse(is.na(PI_pos), 0, PI_pos ),
            ifelse(is.na(Annexin_pos), 0, Annexin_pos))]
}
)

# Assign dead cells with more than 2 pixels to TRUE, rest to FALSE   
inputdata <- lapply( inputdata, function(x) {
  x[ , c("PI_pos", "Annexin_pos"):=
       list(ifelse(PI_pos > 2, TRUE, FALSE ),
            ifelse(Annexin_pos > 2, TRUE, FALSE))]
}
)

return(inputdata)
}

# calculate summaries 
calc_summaries <- function(inputdata) {
  # remove locationID
  inputdata <- lapply(inputdata, function(x) { x[, locationID := NULL]})
    
  lapply(inputdata, function(x) {
    x[ , lapply(.SD, mean, na.rm = TRUE), 
       by = c("treatment", "dose_uM", "CMAX", "cell_line", "plateID", "timeID", "replID") ]
  })

    }

calc_ICAM1_diff <- function(inputdata, ICAM1 = TRUE) {
  if(ICAM1){
    
    
    inputdata <- lapply(inputdata, function(x) {
      x[ , c("GFP_diff1i",  "GFP_diff2i", "GFP_diff3i", "GFP_diff4i")  := 
           list(GFP_pos1i - GFP_neg1i, 
                GFP_pos2i - GFP_neg2i,
                GFP_pos3i - GFP_neg3i,
                GFP_pos4i - GFP_neg4i
                ) ]
    })
    
    
  } else{
    print("Only applicable for ICAM1 reporter")
  }
return(inputdata)
  }

##==## run functions

parsedlistpath <- 'D:/analysis/DILI timepoint/data/'
dir(parsedlistpath)

name = 'BTG2'
cyto = TRUE

loaded_data <- clean_load_sel_count_dose_cmax(parsedlistpath = 'D:/analysis/DILI timepoint/data/', name = 'BIP', cyto = TRUE)

outputdata <- calc_GFP_Counts(inputdata = loaded_data, ICAM = FALSE)
outputdata <- norm_cell_counts(inputdata = outputdata)
outputdata <- calc_cell_death(inputdata = outputdata)
summ_data <- calc_summaries(inputdata = outputdata)

# if ICAM1:
summ_data <- calc_ICAM1_diff(inputdata = summ_data, ICAM = TRUE)

summ_data <- do.call('rbind', summ_data)

write.table(summ_data, file = paste0('../generated/results/processed summaries/', name, '.txt'), sep = '\t', row.names = FALSE)





##==## load and validate summaries: treatment dose_uM CMAX comparison over cell lines/ replicates


##==## basic quality control plotting of all data
