#!/usr/bin/env Rscript
#
# Plot histogram of copy number
#
# plotChrCNAhist.R -l chrLength.txt ASCAT.count -t eps -o cohortA
# R CMD BATCH .R
# source(".R")
#
#library(MASS)
library(optparse)
library(tools)
#
# Copy number plot
#
lineDiffHist <- function( countD, chrStart, k, seq, colorR ) {
  for( i in 1:22 ) {
     chrData <- subset( countD, chr==i )
     if ( length(chrData[,1]) == 0 ) {
       next
     }
     for( j in 1:length(chrData[,1]) ) {
        if ( chrData[j,k] > 0 ) {
          yV = sum( chrData[j,seq], na.rm=TRUE )
          lines( c(chrStart[i]+chrData[j,2], chrStart[i]+chrData[j,3]),
                 c(yV, yV), col=colorR, lwd=2 )
       }
     }
  }
}
#
plotCNAhistogram <- function( chrLength, countD, cT, cR, title, outputPrefix, outType ) {
  cS <- 0.5
  xLabel <- "Chromosome position"
  yLabel <- "Proportion"
  #
  chrStart <- rep( 0, nrow(chrLength)+1 )
  xtickPos <- rep( 0, nrow(chrLength) )
  for( i in 2:length(chrStart) ) {
     j <- i - 1
     chrStart[i] <- chrStart[j] + as.numeric(chrLength[j,2])
     xtickPos[j] <- (chrStart[i] + chrStart[j])/2
  }
  # For x axe
  xmin <- 1
  xmax <- max( chrStart )
  xticks <- chrStart
  # For y axe
  ymin <- 0
  ymax <- 1.0
  #
  outName <- paste0( outputPrefix, ".hist.", outType )
  if ( outType == "svg" ) {
    svg( outName, height=4.5, width=18 )
  } else if ( outType == "png" ) {
    png( outName, height=600, width=1800 )
  } else if ( outType == "pdf" ) {
    pdf( outName, height=4.5, width=18 )
  } else {
    postscript( outName, paper='special', height=4.5, width=18, horizontal=F )
  }
  plot( NA, NA, xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=cS, xaxt="n",
        xlab="", ylab="", main="" )
  mtext( xLabel, side=1, cex=1.0, line=2 )
  mtext( yLabel, side=2, cex=1.0, line=2.4 )
  mtext( title, side=3, cex=1.5, line=1 )
  abline( v=chrStart, col="yellow", cex=cS )
  # xlabel
  par(xpd=TRUE)
  for( i in 1:22 ) {
     text( xtickPos[i], -ymax/14, paste(i) )
  }
  #
  for( k in 4:length(countD[1,]) ) {
    print(c("Ploting ", outputPrefix, k) )
    countN <- countD[countD[,k]>0.0000001,]
    if ( length(countN[,1]) > 0 ) {
      lineDiffHist( countN, chrStart, k, c(4:k), cR[k-3] )
    }
  }
  legend( x=(xmax-0.005*(xmax-xmin)), y=(ymax+0.0*(ymax-ymin)),
          rev(cT), col=rev(cR), lty=1, cex=0.8, bty="n" )
  dev.off()
}
#
# Setting options
#
option_list <- list(
    make_option(c("-l", "--length-of-genome"), type="character", default=FALSE,
                dest="chrLenFile", help="Length of chromosomes"),
    make_option(c("-t", "--type-of-image"), type="character", default="png",
                dest="outType", help="eps|png|svg|pdf"),
    make_option(c("-o", "--output-prefiex"), type="character", default=FALSE,
                dest="prefix", help="header of output file")
)
parser <- OptionParser(usage = "%prog [options] copy-number-count-file", option_list=option_list)
opt <- parse_args(parser, positional_arguments=1)
# Check
if ( length(opt$args) < 1 ) {
   print("Copy number count file must be specified.")
   q(status=1)
}
# Setting
if ( opt$options$prefix == FALSE ) {
   opt$options$prefix <- sub(".count.txt","",file_path_sans_ext(basename(opt$args[1])))
}
#print(opt)
#print(paste(opt$args[1],opt$options$prefix,opt$options$outType))
# 
# Reading data
#
#countList <- read.table( "LM_Tumor.count", header=TRUE )
countList <- read.table( opt$args[1], header=TRUE )
#> colnames(countList)
# [1] "chr"   "start" "end"   "m2"    "m1"    "p1"    "p2"    "p3"    "o2LOH"
#[10] "loss"  "gain"
chrs = paste(c(1:22))
countList <- countList[countList[,1]%in%chrs,]
#countList$chr <- as.numeric(countList$chr)
if ( opt$options$chrLenFile == FALSE ) {
  chrLength <- cbind(chrs,NA)
  for( i in 1:length(chrs) ) {
    chrLength[i,2] <- max(countList[countList$chr==chrs[i],3], na.rm=T)
  }
} else {
  #chrLength <- read.table( "../../ascat-2.5.2/gcProcessing/hg19.chrom.sizes", header=F )
  chrLength <- read.table( opt$options$chrLenFile, header=F )
  chrLength[,1] <-  gsub("chr","",chrLength[,1])
  chrLength <- chrLength[chrLength[,1]%in%chrs,]
}
sL <- sort(as.numeric(chrLength[,1]), index.return = TRUE)
chrLength <- chrLength[sL$ix,]
#
# Loss
#
cT <- c( "0", "1", ">=2 LOH" )
cR <- c( "blue", "cyan", "green" )
title <- paste( "Histogram of copy number (loss):", opt$options$prefix )
outputPrefix <- paste0( opt$options$prefix, "Loss" )
plotCNAhistogram( chrLength, countList[,c(1,2,3,4,5,9)], cT, cR, title, outputPrefix, opt$options$outType )
#
# Gain
#
cT <- c( ">=5", "4", "3" )
cR <- c( "brown", "red", "orange" )
title <- paste( "Histogram of copy number (gain):", opt$options$prefix )
outputPrefix <- paste0( opt$options$prefix, "Gain" )
plotCNAhistogram( chrLength, countList[,c(1,2,3,8,7,6)], cT, cR, title, outputPrefix, opt$options$outType )
