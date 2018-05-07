# dose response script for DILI data


# first figure to generate: per compound dose response fingerprints. normalize all features per plate relative to DMSO/DMEM t0
# then summarize the time dimension to meaningfull features

options(stringsAsFactors = FALSE)

library(data.table)
library(ggplot2)
library(grid)
library(plyr)
require(stringr)
require(reshape2)
require(NMF)

rm(list=ls())

theme_sharp<-function(base_size = 12, base_family = "serif") {
  theme(line = element_line(colour = "black", size = 0.5, linetype = 1,
                            lineend = "butt"),
        rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
        text = element_text(family = base_family, face = "plain",
                            colour = "black", size = base_size,
                            hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9),
        axis.text = element_text(size = rel(0.8), colour = "black"),
        strip.text = element_text(size = rel(0.8)),
        axis.line = element_line(),
        axis.text.x = element_text(vjust = 1),
        axis.text.y = element_text(hjust = 1),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(vjust=0),
        axis.title.y = element_text(angle = 90,vjust=0.3),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks.margin = unit(0.1, "cm"),
        legend.background = element_rect(colour = "black",size=0.25),
        legend.margin = unit(0.2, "cm"),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.key.size = unit(1.2, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(size = rel(0.8)),
        legend.text.align = NULL,
        legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0),
        legend.title.align = 0.5,
        legend.position = "right",
        legend.direction = NULL,
        legend.justification = "center",
        legend.box = NULL,
        
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey",size=0.5),
        panel.grid.minor = element_line(colour = "grey", size = 0.25),
        panel.margin = unit(0.25, "lines"),
        
        strip.background = element_rect(fill = "white", colour = NA),
        strip.text.x = element_text(),
        strip.text.y = element_text(angle = -90),
        
        plot.background = element_rect(colour = "white"),
        plot.title = element_text(size = rel(1.8),vjust=1),
        plot.margin = unit(c(1, 1, 0.5, 0.5), "lines"),
        
        complete = TRUE
  )
}

setwd("J:/Workgroups/FWN/LACDR/TOX/Data Annaloes/Data Johnny/Heatmap/")
if(FALSE) setwd('/Users/Wouter/mnt/win_shares/Public/Workgroups/FWN/LACDR/TOX/Data Annaloes/Data Johnny/Heatmap/')   # Wouter

raw.data =list()
dir.files <- dir("summary files")
dir.files <- paste("summary files/",dir.files[grepl("(.txt)$", dir.files)], sep = "")


for (i in seq_along(dir.files)){
  raw.data[[i]] <- read.table(file = dir.files[i], header = TRUE, sep ="\t")
  raw.data[[i]]$variable  <- colnames(raw.data[[i]])[ ncol(raw.data[[i]])]
  colnames(raw.data[[i]])[ ncol(raw.data[[i]])-1] <- "value"
  }
raw.data <- do.call('rbind', raw.data)
head(raw.data)
tail(raw.data)



unique(raw.data$variable)
raw.data$variable <- gsub("count_AnV_masked_primaryID_AreaShape_Area_larger_0_", "AnVFraction", raw.data$variable)
raw.data$variable <- gsub("count_PI_masked_primaryID_AreaShape_Area_larger_0_", "PIFraction", raw.data$variable)
raw.data$variable <- gsub("Nuclei_Intensity_IntegratedIntensity_image_GFP", "IntegratedIntensity", raw.data$variable)
raw.data$variable <- gsub("Cytoplasm_Intensity_IntegratedIntensity_image_GFP", "IntegratedIntensity", raw.data$variable)
raw.data$variable <- gsub("Nuclei_Intensity_MeanIntensity_image_GFP", "MeanIntensity", raw.data$variable)
raw.data$variable <- gsub("Cytoplasm_Intensity_MeanIntensity_image_GFP", "MeanIntensity", raw.data$variable)
raw.data$variable <- gsub("imageCountParentObj", "CellCount", raw.data$variable)



# fix cell line names
unique(raw.data$cell_line)

#zorg dat zelfde concentratie echt zelfde concentatie is (40.0!= 40.00)
unique(raw.data$dose_uM)



raw.data<- as.data.table(raw.data)

## eerst Int Int FC DMSO en normalizeren 0-1 per plateID

raw.data.Int<-subset(raw.data, raw.data$variable=="IntegratedIntensity")  

unique(raw.data.Int$plateID)

dataIn = raw.data.Int
normVar = "IntegratedIntensity"
buffer <- dataIn[ variable%in% normVar ] # select 1 variable


Checkcontrol <- buffer[which(buffer$treatment == 'DMSO' & buffer$dose_uM == 0.2)] # select the DMSO control
if(any(colnames(dataIn) %in% "timeID")){ # if time course data
  controlV<-Checkcontrol[ timeID==24 ,mean(value, na.rm=TRUE), by = plateID]   #plateID is echt een unique plaat
} else { # single time point
  controlV<-Checkcontrol[  ,mean(value, na.rm=TRUE), by = plateID] 
}
unique(controlV$plateID)

setkey(controlV, "plateID")
setkey(buffer, "plateID")
buffer<- buffer[controlV]
buffer[, FC_value:= value/V1 ]  #V1 is controlV value
buffer[, value:=NULL]
buffer[, variable:= paste("mmnFC", normVar, sep ="_")]


buffer[ ,plateID:= gsub("_rep[1-30]{1,2}", "", plateID) ] # plateID zonder rep erin
buffer<-buffer[  , mean(FC_value,na.rm=TRUE), by = c("treatment", "variable", 
                                                        "dose_uM", "cell_line", "timeID" , "plateID")] # gemiddelde van de reps voor FC_value


buffer[, plateID:=NULL] #gemiddelde van reps is genomen en timeID is er nog dus plateID kan weg

setnames(buffer, old="V1", new="mFC_value") #V1 was nu gemiddelde FC --> hernoemd

minValue <- buffer[,min(mFC_value, na.rm=TRUE), by="cell_line"] # mininum waarde per cellijn
maxValue <- buffer[,max(mFC_value, na.rm=TRUE), by="cell_line"] # maximum waarde per cellijn
setkey(minValue, "cell_line")
setkey(maxValue, "cell_line")
setkey(buffer, "cell_line")
buffer<- buffer[minValue]
setnames(buffer, "V1", "minValue")
buffer<- buffer[maxValue]
setnames(buffer, "V1", "maxValue")
buffer[, mmnvalue:= (mFC_value - minValue)/ 
               (maxValue-minValue)]           # min max normalizatie

setnames(buffer, "mmnvalue", "value")
buffer[, mFC_value:=NULL]
buffer[, minValue:=NULL]
buffer[, maxValue:=NULL]      #alles weg behalve de DMSO FC en dan min max normalizatie waarde --> is nu value gemaakt al 
buffer[, timeID:=NULL] ## alleen 72h dus geen time ID nodig


raw.data.bInt <- buffer # buffer weer gewoon raw.data noemen. bInt = bewerkte Integrated Intensity

unique(raw.data.bInt$cell_line)


## zelfde voor Mean Intensity eerst Meant Int FC DMSO en normalizeren 0-1 per plateID
raw.data.MeanInt<-subset(raw.data, raw.data$variable=="MeanIntensity")  


dataIn = raw.data.MeanInt
normVar = "MeanIntensity"
buffer <- dataIn[ variable%in% normVar ] # select 1 variable

Checkcontrol <- buffer[which(buffer$treatment == 'DMSO' & buffer$dose_uM == 0.2)] # select the DMSO control
if(any(colnames(dataIn) %in% "timeID")){ # if time course data
  controlV<-Checkcontrol[ timeID==24 ,mean(value, na.rm=TRUE), by = plateID]   #plateID is echt een unique plaat
} else { # single time point
  controlV<-Checkcontrol[  ,mean(value, na.rm=TRUE), by = plateID] 
}

setkey(controlV, "plateID")
setkey(buffer, "plateID")
buffer<- buffer[controlV]
buffer[, FC_value:= value/V1 ]  #V1 is controlV value
buffer[, value:=NULL]
buffer[, variable:= paste("mmnFC", normVar, sep ="_")]


buffer[ ,plateID:= gsub("_rep[1-30]{1,2}", "", plateID) ] # plateID zonder rep erin
buffer<-buffer[  , mean(FC_value,na.rm=TRUE), by = c("treatment", "variable", 
                                                     "dose_uM", "cell_line", "timeID" , "plateID")] # gemiddelde van de reps voor FC_value
buffer[, plateID:=NULL] #gemiddelde van reps is genomen en timeID is er nog dus plateID kan weg

setnames(buffer, old="V1", new="mFC_value") #V1 was nu gemiddelde FC --> hernoemd


minValue <- buffer[,min(mFC_value, na.rm=TRUE), by="cell_line"] # mininum waarde per cellijn
maxValue <- buffer[,max(mFC_value, na.rm=TRUE), by="cell_line"] # maximum waarde per cellijn
setkey(minValue, "cell_line")
setkey(maxValue, "cell_line")
setkey(buffer, "cell_line")
buffer<- buffer[minValue]
setnames(buffer, "V1", "minValue")
buffer<- buffer[maxValue]
setnames(buffer, "V1", "maxValue")
buffer[, mmnvalue:= (mFC_value - minValue)/ 
         (maxValue-minValue)]           # min max normalizatie

setnames(buffer, "mmnvalue", "value")
buffer[, mFC_value:=NULL]
buffer[, minValue:=NULL]
buffer[, maxValue:=NULL]      #alles weg behalve de DMSO FC en dan min max normalizatie waarde --> is nu value gemaakt al
buffer[, timeID:=NULL]

raw.data.bMeanInt <- buffer # buffer weer gewoon raw.data noemen. bInt = bewerkte Integrated Intensity


## mean van count maken, per cell line en tijdspunt
raw.data.Count<-subset(raw.data, raw.data$variable=="countPositives")  

raw.data.mCount<-raw.data.Count[  , mean(value,na.rm=TRUE), by = c("treatment", "variable", 
                                                                   "dose_uM", "cell_line" )]

setnames(raw.data.mCount, "V1", "value")

## mean van count en bewerkte Integrated Intensity en bewerkte Mean Intensity samenvoegen
raw.data.allInt <- rbind(raw.data.bInt, raw.data.bMeanInt,raw.data.mCount)


unique(raw.data.allInt$variable)


##check min/max 
min(raw.data.allInt[ variable == "mmnFC_IntegratedIntensity" & cell_line=="HepG2 SRXN1" , value], na.rm=TRUE)


## mean van PI en AnV over alle cellijnen berekenen
raw.data.AnVPI<-subset(raw.data, raw.data$variable %in% c("AnVFraction", "PIFraction"))  
raw.data.mAnVPI<-raw.data.AnVPI[  , mean(value,na.rm=TRUE), by = c("treatment", "variable", 
                                                     "dose_uM")] 

#test#  raw.data.test<-subset(raw.data.total, raw.data.total$treatment=="DEM")


setnames(raw.data.mAnVPI, "V1", "value")



## variable en cell_line aan elkaar --> anders een kolom te veel vergeleken met PI/AnV
raw.data.allInt$variable[ raw.data.allInt$variable %in% "countPositives" | raw.data.allInt$variable %in% "mmnFC_IntegratedIntensity" | raw.data.allInt$variable %in% "mmnFC_MeanIntensity"]<-
  paste(raw.data.allInt$variable[ raw.data.allInt$variable %in% "countPositives"| raw.data.allInt$variable %in% "mmnFC_IntegratedIntensity" | raw.data.allInt$variable %in% "mmnFC_MeanIntensity"]
        , raw.data.allInt$cell_line[raw.data.allInt$variable %in% "countPositives"| raw.data.allInt$variable %in% "mmnFC_IntegratedIntensity" | raw.data.allInt$variable %in% "mmnFC_MeanIntensity"], sep="_")

raw.data.allInt[, cell_line:=NULL] ## cell_line kolom weghalen

## gemiddelde mm cell count
#raw.data.cellcount<-subset(raw.data, raw.data$variable=="mmCellCount")  
#raw.data.mcellcount<-raw.data.cellcount[  , mean(value,na.rm=TRUE), by = c("treatment", "variable", 
                                                                 # "dose_uM", "cell_line" )]
#setnames(raw.data.mcellcount, "V1", "value")

#raw.data.mcellcount$variable[ raw.data.mcellcount$variable %in% "mmCellCount"]<-
 # paste(raw.data.mcellcount$variable[ raw.data.mcellcount$variable %in% "mmCellCount"]
  #      , raw.data.mcellcount$cell_line[raw.data.mcellcount$variable %in% "mmCellCount"], sep="_")
#raw.data.mcellcount[, cell_line:=NULL]


## min max van cell count
raw.data.CC<-subset(raw.data, raw.data$variable=="CellCount")  

unique(raw.data.Int$plateID)

raw.data.CC[ ,plateID:= gsub("_rep[1-30]{1,2}", "", plateID) ] # plateID zonder rep erin
raw.data.CC<-raw.data.CC[  , mean(value,na.rm=TRUE), by = c("treatment", "variable", 
                                                     "dose_uM", "cell_line", "timeID" )] # gemiddelde van de reps 

setnames(raw.data.CC, old="V1", new="m_value") 

minValue <- raw.data.CC[,min(m_value, na.rm=TRUE), by="cell_line"] # mininum waarde per cellijn
maxValue <- raw.data.CC[,max(m_value, na.rm=TRUE), by="cell_line"] # maximum waarde per cellijn
setkey(minValue, "cell_line")
setkey(maxValue, "cell_line")
setkey(raw.data.CC, "cell_line")
raw.data.CC<- raw.data.CC[minValue]
setnames(raw.data.CC, "V1", "minValue")
raw.data.CC<- raw.data.CC[maxValue]
setnames(raw.data.CC, "V1", "maxValue")
raw.data.CC[, mmnvalue:= (m_value - minValue)/ 
         (maxValue-minValue)]           # min max normalizatie

setnames(raw.data.CC, "mmnvalue", "value")
raw.data.CC[, m_value:=NULL]
raw.data.CC[, minValue:=NULL]
raw.data.CC[, maxValue:=NULL]      #alles weg behalve de DMSO FC en dan min max normalizatie waarde --> is nu value gemaakt al 
raw.data.CC[, timeID:=NULL] ## alleen 72h dus geen time ID nodig


raw.data.mmCC <- raw.data.CC # buffer weer gewoon raw.data noemen. bInt = bewerkte Integrated Intensity
unique(raw.data.mmCC$cell_line)


raw.data.mmCC$variable[ raw.data.mmCC$variable %in% "CellCount" ]<-
  paste(raw.data.mmCC$variable[ raw.data.mmCC$variable %in% "CellCount"]
        , raw.data.mmCC$cell_line[raw.data.mmCC$variable %in% "CellCount"], sep="_")

raw.data.mmCC[, cell_line:=NULL]

raw.data.mmCC$variable <- gsub("CellCount", "mmCellCount", raw.data.mmCC$variable)

## gemiddelde mm cell count alle cell lijnen samen
#raw.data.mtcellcount<-raw.data.tcellcount[  , mean(value,na.rm=TRUE), by = c("treatment", "variable", 
#                                                                   "dose_uM")] 
#setnames(raw.data.mtcellcount, "V1", "value")


##  IntegratedIntensity, Countpositives, MeanIntensity, relative cell count in 1 tabel
raw.data.semitotal <- rbind(raw.data.allInt,raw.data.mmCC)

## AnV en PI  1 tabel erbij
raw.data.total <- rbind(raw.data.semitotal,raw.data.mAnVPI )

unique(raw.data.total$variable)

unique(raw.data.total$treatment)

### controles weghalen want die wil BOB toch niet in de heatmap..
raw.data.total<-raw.data.total[!treatment== "DEM"]
raw.data.total<-raw.data.total[!treatment== "Azathioprine"]
raw.data.total<-raw.data.total[!treatment== "Thapsigargin"]
#raw.data.total<-raw.data.total[!treatment== "Tunicamycin"]
raw.data.total<-raw.data.total[!treatment== "Omeprazol"]
raw.data.total<-raw.data.total[!treatment== "Cisplatine"]
raw.data.total<-raw.data.total[!treatment== "Etoposide"]
raw.data.total<-raw.data.total[!treatment== "DMSO"]
raw.data.total<-raw.data.total[!treatment== "mock"]

unique(raw.data.total$treatment)


## dose levels maken voor Heatmap (omdat je niet overal precies dezelfde doses hebt)
raw.data.df <- as.data.frame(raw.data.total)


raw.data.df <- raw.data.df[  order( raw.data.df[ , "variable"],
                                    raw.data.df[ , "treatment"],
                                    raw.data.df[ , "dose_uM"]  
),]


unique(raw.data.df$variable)
counts.d <- ddply(raw.data.df, .(variable, treatment ), summarize,count.d.l = length(dose_uM))
head(counts.d)    #--> sommige compounds meer of minder dan 10 concentraties...
raw.data <- as.data.table(raw.data.df)

raw.data$dose.f <- 'kaas'
setkeyv(raw.data, c("treatment", "variable"))

for (i in 1 : nrow(counts.d))
{
  
  raw.data[ list(counts.d$treatment[i],  counts.d$variable[i]),
            dose.f:=gl(counts.d$count.d.l[i], 1)]
}


all.data.w <- dcast(data= raw.data,  treatment + variable ~as.numeric(dose.f), value.var = "value")
head(all.data.w)
tail(all.data.w)
#raw.data.test<-subset(all.data.w, all.data.w$treatment=="Alpidem")


## cluster the compounds based on similarity of fingerprint dose response vectors
all.data.w<- as.data.table(all.data.w)
colnames(all.data.w)[3:12] <- paste("dose", 1:10, sep = "_")

all.data.w[, treatment:=factor(treatment)]
longVecsPerComp <- split(all.data.w, all.data.w$treatment)
(longVecsPerComp[[1]])



gatherCompVecs = list()
for( j in seq_along(longVecsPerComp)){
  cat('j=', j , '\n')
  gathervecs = list() ### deze hier want niet elke compound heeft evenveel variabelen
  for(i in 1:nrow(longVecsPerComp[[j]])){
    cat('i=', i , '\n')
    variable_bla <- longVecsPerComp[[j]][i, variable]
    gathervecs[[i]] <- longVecsPerComp[[j]][i,  list(dose_1, dose_2, dose_3, dose_4, dose_5,  # selecteer kolom namen annotatie in data.table
                                                     dose_6, dose_7, dose_8, dose_9, dose_10)]
    colnames(gathervecs[[i]]) <- paste(variable_bla, colnames(gathervecs[[i]]))
    
  }
  gatherCompVecs[[j]] <- as.data.frame(t(unlist(gathervecs)))
  row.names(gatherCompVecs[[j]]) <- c(as.character(unique(longVecsPerComp[[j]]$treatment)))
  names(gatherCompVecs)[[j]] <- as.character(unique(longVecsPerComp[[j]]$treatment))
  }

length(gatherCompVecs) # elke entry in deze lijst bevat de vectoren van de verschillende metingen

######!!!!###2016-07-08#####
gatherCompVecsDF <- rbind.fill(gatherCompVecs)
#gatherCompVecsDF <- do.call(rbind, gatherCompVecs)   # by wouter

# loses the row.names, set again
rownames(gatherCompVecsDF) <- unlist(lapply(gatherCompVecs, function(i) {rownames(i)}))


#########################################################################################################
#########################################################################################################
#### The following function is written/adapted by Wouter den Hollander - 18 July 2106

## clustering function with slick feature selection
performClustering <- function(data, returnTree = FALSE, reporter = c('CHOP', 'SRXN1', 'p21'), feature = c('AnVFraction', 'PIFraction', 'countPositives', 'mmCellCount', 'mmnFC_IntegratedIntensity', 'mmnFC_MeanIntensity'), clusterMethod = c('euclidean', 'correlation'), .debug = FALSE) {
  if (.debug) {
    data <- gatherCompVecsDF
    reporter <- c('CHOP', 'SRXN1', 'p21')
    feature <- c('AnVFraction', 'PIFraction', 'countPositives', 'mmCellCount', 'mmnFC_IntegratedIntensity', 'mmnFC_MeanIntensity')
    clusterMethod <- c('euclidean')
  }

  ### Remarks
  ## 1) This function is very script specific, given the format of the data, reporters and features.

  ### Arguments
  ## returnTree:  logical that indicates whether the dendrogram OR the heatmap row order should be returned
  
  ### Values
  ## Returns the cluster dendrogram if returnTrue is set to TRUE
  ## Returns the order of the compounds by which is clustered. Clustering is done only on the selected reporter * feature.

  ### Function
  # Checks / parsing
    # Check clusterMethod
      if(length(clusterMethod) != 1) stop('Select a single method to cluster by.\n')
      if(!(grepl('euclidean', clusterMethod) | grepl('correlation', clusterMethod))) stop('Select either euclidean or correlation as clusterMethod.\n')

    # Parse 
      reporterSelection <- lapply(reporter, function(x) grep(x, colnames(data)))
        names(reporterSelection) <- reporter
      featureSelection  <- lapply(feature, function(x) grep(x, colnames(data)))
        names(featureSelection) <- feature

      fractionFeatureInd <- grep('Fraction', names(featureSelection))
      reporterSpecificFeatureInd <- which(!grepl('Fraction', names(featureSelection)))
      columnSelection <- c(unlist(featureSelection[fractionFeatureInd]),
                           intersect(unlist(featureSelection[reporterSpecificFeatureInd]), unlist(reporterSelection)))

      if(length(columnSelection) == 0) stop('Something went wrong with selecting the data to cluster on!\n')

    # data matrix (superstitious)
      data <- as.matrix(data) 

    # Subsetted data to cluster on
      dataSubset <- data[, columnSelection]

      rowsWithOnlyNA <- which(rowSums(is.na(dataSubset)) == ncol(dataSubset))
        if(length(rowsWithOnlyNA) > 0) {
          cat(paste0('The following treatments only contain NA among the selected features, and will be discarded:\n', names(rowsWithOnlyNA), '\n'))
          data <- data[-rowsWithOnlyNA, ]
        }
     
      dataSubset <- data[, columnSelection]

  # Do the magic
    if(clusterMethod == 'euclidean') {
      matrixToClusterOn <- dist(dataSubset)
    }

    if(clusterMethod == 'correlation') {
       matrixToClusterOn <- dist(cor(t(dataSubset), use = 'pairwise.complete.obs'))
    }

    clustering <- hclust(matrixToClusterOn, method = 'ward.D')
    rowOrdering <- clustering$labels[clustering$order]

  # Return
    if(returnTree) {
        return(clustering)
      } else {
        return(rowOrdering)
      }
}


## Wrapper functions: less prone to c/p errors & writes plot to PDF
CreateHeatmap <- function(reporter, feature, clusterMethod) {
  ## Forced column order & heatmap annotations/aesthetics
  colOrder <- match(read.table(file= 'column names.txt', sep = '\t')[, 1], 
              colnames(gatherCompVecsDF))

  require(RColorBrewer)

  columnAnnotation <- do.call(rbind, strsplit(colnames(gatherCompVecsDF)[colOrder], 'dose_'))
  columnAnnotation <- list(dose = as.numeric(columnAnnotation[, 2]), 
                           feature = columnAnnotation[, 1])

  columnAnnotationColors <- list(dose = c(brewer.pal(9, "BuPu"), '#000000'),
                                 feature = rainbow(length(unique(columnAnnotation$feature))))
  
  colBreaks <- c(seq(-0.01, 0.3,length.out=15), seq(0.301,1.01, length.out = 14))
  heatmapColors <- brewer.pal(9, "YlOrRd")

  clustering <- performClustering(data = gatherCompVecsDF, 
                                         reporter = reporter,
                                         feature = feature,
                                         clusterMethod = clusterMethod)

  rowColorAnnotation <- read.table(file = 'rowAnot.txt', sep ="\t", header = FALSE, row.names = 1) 
  rowColorAnnotation <- rowColorAnnotation[clustering, , drop = FALSE]
    colnames(rowColorAnnotation) <- 'DILI'

  pdfFile <- paste0(clusterMethod, 
                    '__', 
                    paste0(reporter, collapse = '_'), 
                    '__', 
                    paste0(feature, collapse = '_'), '.pdf') 

  pdf(pdfFile, width = 18, height = 22)
    aheatmap(gatherCompVecsDF[clustering, colOrder], Rowv = NA, Colv = NA, 
             color = heatmapColors, breaks = colBreaks,
             fontsize = 8, width = 10, height = 10, legend = TRUE,
             annRow = rowColorAnnotation , annCol = columnAnnotation, annColors = columnAnnotationColors,
             cexRow= 2, cexCol = 2)
  dev.off()

  cat('Saved heatmap to\n', paste0(getwd(), '/', pdfFile, '.\n'))
}

CreateDendrogram <- function(reporter, feature, clusterMethod) {
  dendrogram <- performClustering(data = gatherCompVecsDF, 
                                         reporter = reporter,
                                         feature = feature,
                                         clusterMethod = clusterMethod,
                                         returnTree = TRUE)

  pdfFile <- paste0('dendrogram__',
                    clusterMethod, 
                    '__', 
                    paste0(reporter, collapse = '_'), 
                    '__', 
                    paste0(feature, collapse = '_'), 
                    '.pdf') 

  pdf(pdfFile, width = 22, height = 6)
    plot(dendrogram)
  dev.off()

  cat('Saved dendrogram to\n', paste0(getwd(), pdfFile, '.\n'))
}


## Examples
# Cluster by euclidean distance, only on IntegratedIntensity from CHOP.
CreateHeatmap(reporter = 'CHOP',
              feature = 'mmnFC_IntegratedIntensity',
              clusterMethod = 'euclidean')

CreateDendrogram(reporter = 'CHOP',
              feature = 'mmnFC_IntegratedIntensity',
              clusterMethod = 'euclidean')


# Cluster by correlation, on IntegratedIntensity & MeanIntensity from CHOP & SRXN1.
CreateHeatmap(reporter = c('CHOP', 'SRXN1'),
              feature = c('mmnFC_IntegratedIntensity', 'mmnFC_MeanIntensity'),
              clusterMethod = 'correlation')

CreateDendrogram(reporter = c('CHOP', 'SRXN1'),
              feature = c('mmnFC_IntegratedIntensity', 'mmnFC_MeanIntensity'),
              clusterMethod = 'correlation')


# Cluster by correlation, on all data (select all features & reporters)
CreateHeatmap(reporter = c('CHOP', 'SRXN1', 'p21'),
              feature = c('AnVFraction', 'PIFraction', 'countPositives', 'mmCellCount', 'mmnFC_IntegratedIntensity', 'mmnFC_MeanIntensity'),
              clusterMethod = 'correlation')

CreateDendrogram(reporter = c('CHOP', 'SRXN1', 'p21'),
              feature = c('AnVFraction', 'PIFraction', 'countPositives', 'mmCellCount', 'mmnFC_IntegratedIntensity', 'mmnFC_MeanIntensity'),
              clusterMethod = 'correlation')





