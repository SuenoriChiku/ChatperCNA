#!/usr/bin/env Rscript
#
# Plot copy numbers for segment file by order
#
# colorOrderCNA.R  -l chrLength.txt cohortADendrogramEuclWard.order ASCAT.segments -t eps -o cohortA
# % R CMD BATCH .R
# #source(".R")
#library(MASS)
library(optparse)
library(tools)
#
#
# Copy number plot
#
colorCopy <- function( seg, chrStart, ybase, cR ) {
  for( i in 1:22 ) {
     chrData <- subset( seg, chr==i )
     if ( length(chrData[,1]) > 0 ) {
       for( j in 1:length(chrData[,1]) ) {
         cI <- chrData[j,10]+1
         lines( c(chrStart[i]+chrData[j,3], chrStart[i]+chrData[j,4]),
                c(ybase, ybase), col=cR[cI], lty=1, lwd=3 )
	 if ( chrData[j,11]==1 ) {
            lines( c(chrStart[i]+chrData[j,3], chrStart[i]+chrData[j,4]),
                   c(ybase, ybase), col=cR[7], lty=1, lwd=1 )
         }
       }
     }
  }
}
#
plotIndividualCopyNumber <- function( chrLength, ID, ASCATSeg, cT, cR, outputPrefix, outType ){
  titleName <- paste( "Copy numbers:", outputPrefix )
  xLabel <- "Chromosome position"
  yLabel <- "Sample"
  cS <- 0.5
  #
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
  ytickLabels <- ID
  ymin <- 0
  ymax <- length(ytickLabels)
  yticks <- 1:ymax
  #
  outName <- paste0( outputPrefix, ".OrderCNA.", outType )
  h <- 8
  if ( length(ID) > 110 ) {
    h <- 10
  }
  if ( outType == "svg" ) {
    svg( outName, height=h, width=18 )
  } else if ( outType == "png" ) {
    png( outName, height=(h*100), width=1800 )
  } else if ( outType == "pdf" ) {
    pdf( outName, height=h, width=18 )
  } else {
    postscript( outName, paper='special', height=h, width=18, horizontal=F )
  }
  plot( NA, NA, xlim=c(xmin,xmax), ylim=c(ymin,ymax), cex=cS, axes=FALSE,
        xlab="", ylab="", main="" )
  mtext( xLabel, side=1, cex=1.0, line=1 )
  mtext( yLabel, side=2, cex=1.0, line=1 )
  mtext( titleName, side=3, cex=1.5, line=1 )
  par(xpd=FALSE)
  abline( v=chrStart, col="yellow", cex=cS )
  abline( h=yticks,   col="lightgray", cex=cS )
  par(xpd=TRUE)
  # xlabel
  yPos <- ymin - (ymax-ymin)*0.05
  print( yPos )
  for( i in 1:22 ) {
     text( xtickPos[i], yPos, paste(i) )
  }
  # ylabel
  #axis( side=2, at=yticks, labels=ytickLabels )
  xV <- - chrStart[23]/40 # 80
  cS <- 0.5
  if ( ymax <= 30 ) {
    xV <- - chrStart[23]/40 # 40
    cS <- 0.8
  } else if ( ymax <= 50 ) {
    xV <- - chrStart[23]/40 # 40
    cS <- 0.75
  }
  for( i in yticks ) {
     text( xV, yticks[i], ytickLabels[i], cex=cS )
  }
  #
  for( i in yticks ) {
    aSeg <- ASCATSeg[paste(ASCATSeg[,1])==ID[i],]
    if ( length(aSeg[,1]) == 0 ) {
       next
    }
    print( c(i,ID[i]) )
    colorCopy( aSeg, chrStart, yticks[i], cR )  
  }
  legend( x=(xmax+0.00*(xmax-xmin)), y=(ymax+0.14*(ymax-ymin)),
  	  cT, col=cR, lty=1, lwd=2, cex=0.8, bty="n")
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
parser <- OptionParser(usage = "%prog [options] clustering-order-file CNA-segment-file", option_list=option_list)
opt <- parse_args(parser, positional_arguments=2)
# Check
if ( length(opt$args) < 1 ) {
   print("Copy number count file must be specified.")
   q(status=1)
}
# Setting
if ( opt$options$prefix == FALSE ) {
   opt$options$prefix <- sub("DendrogramEuclWard.order","",file_path_sans_ext(basename(opt$args[1])))
}
#print(opt)
#print(paste(opt$args[1],opt$options$prefix,opt$options$outType))
# 
# Reading data
#
#orderDen <- read.table( "LM_TumorDendrogramEuclWard.order", header=TRUE )
#ASCATSeg <- read.table( "LM_Tumor.segments", header=TRUE )
orderDen <- read.table( opt$args[1], header=TRUE )
ASCATSeg <- read.table( opt$args[2], header=TRUE )
#> colnames(ASCATSeg)
#[1] "sample"   "chr"      "startpos" "endpos"   "nMajor"   "nMinor"
#colnames(ASCATSeg)
# [1] "ID"         "chr"        "start"      "end"        "startSNP"  
# [6] "endSNP"     "numSNPs"    "nA"         "nB"         "copyNumber"
#[11] "o2LOH"
#
chrs = paste(c(1:22))
ASCATSeg <- ASCATSeg[ASCATSeg[,2]%in%chrs,]
if ( ncol(ASCATSeg) < 7  ) {
   ASCATSeg[,7:11] <- NA
   ASCATSeg[,8:9] <- ASCATSeg[,5:6]
   ASCATSeg[,10] <- ASCATSeg[,8] + ASCATSeg[,9]
   ASCATSeg[,11] <- 0
   ASCATSeg[ASCATSeg[,8] > 1 & ASCATSeg[,9] == 0,11] <- 1
}
#
if ( opt$options$chrLenFile == FALSE ) {
  chrLength <- cbind(chrs,NA)
  for( i in 1:length(chrs) ) {
    chrLength[i,2] <- max(ASCATSeg[ASCATSeg$chr==chrs[i],4], na.rm=T)
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
#ID <- sub("\\.", "-", orderDen$x)
ID <- orderDen$x
cmax <- 5
ASCATSeg[ASCATSeg[,10] > cmax,10] <- cmax
cT <- c("0", "1", "2", "3", "4", ">5", ">=2 LOH")
cR <- c("darkblue", "cyan", "lightgray", "orange", "red", "brown", "green")
plotIndividualCopyNumber( chrLength, ID, ASCATSeg, cT, cR, opt$options$prefix, opt$options$outType )
