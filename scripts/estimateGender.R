#!/usr/bin/env Rscript
#
# Estimate gender
#
# estimateGender.R -l Germline_LRR.txt -t eps -o cohortA Germline_BAF.txt 
# R --vanilla --quiet --args  < estimateGender.R
# R CMD BATCH .R
# source(".R")
#
library(MASS)
library(optparse)
library(tools)
#
# PropX vs PropY
#
plotPropXY <- function( X, Y, gender, yLab, outputPrefix, outType ) {
  cS <- 0.8
  xLab <- "Proportion of hetero on X chromosome"
  title <- paste("Hetero X vs", yLab)
  lePos <- "topright"
  cT <- c("Male", "Female")
  cR <- c("blue", "red")
  cP <- c(1, 2)
  #
  xmin <- min( X, na.rm=TRUE )
  xmax <- max( X, na.rm=TRUE )
  ymin <- min( Y, na.rm=TRUE )
  ymax <- max( Y, na.rm=TRUE )
  #
  outName <- paste0( outputPrefix, ".gender.", outType )
  if ( outType == "svg" ) {
    svg( outName, height=6, width=6)
  } else if ( outType == "png" ) {
    png( outName, height=600, width=600 )
  } else if ( outType == "pdf" ) {
    pdf( outName, height=6, width=6)
  } else {
    postscript( outName, paper='special', height=6, width=6, horizontal=F)
  }
  #
  plot(NA,NA,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xLab,ylab=yLab,main=title)
  mList <- gender == "XY"
  points( X[mList], Y[mList], col=cR[1], pch=cP[1], cex=cS )
  points( X[!mList], Y[!mList], col=cR[2], pch=cP[2], cex=cS )
  legend( lePos, cT, col=cR, pch=cP, cex=1.0, bty="n" )
  dev.off()
}
#
histPropX <- function( X, gender, outputPrefix, outType ) {
  xLab <- "Proportion of hetero on X chromosome"
  yLab1 <- "Count of male"
  yLab2 <- "Count of female"
  title <- "Histogram of proportion of hetero call on X"
  lePos <- "top"
  cT <- c("Male", "Female")
  cR <- c("blue", "red")
  #
  xmin <- min( X, na.rm=TRUE )
  xmax <- max( X, na.rm=TRUE )
  #
  outName <- paste0( outputPrefix, ".gender.", outType )
  if ( outType == "svg" ) {
    svg( outName, height=6, width=6)
  } else if ( outType == "png" ) {
    png( outName, height=600, width=600 )
  } else if ( outType == "pdf" ) {
    pdf( outName, height=6, width=6)
  } else {
    postscript( outName, paper='special', height=6, width=6, horizontal=F)
  }
  #
  par( ps=14, mar=c(par("mar")[1:3],3.5) )
  mList <- gender == "XY"
  truehist( X[mList], prob=FALSE, col=cR[1], nbins="FD",
            xlim=c(xmin,xmax), xlab=xLab, ylab=yLab1, main=title )
  par(new=T)
  truehist( X[!mList], prob=FALSE, col=cR[2], nbins="FD", axes=F,
            xlim=c(xmin,xmax), xlab="", ylab="", main="" )
  mtext( side=4, line=2.5, text=yLab2 )
  axis(4)
  legend( lePos, cT, col=cR, lty=1, cex=1.0, bty="n" )
  dev.off()
}
#
# Reading data
#
option_list <- list(
    make_option(c("-l", "--lrr"), type="character", default=FALSE,
                dest="lrr", help="Y chromosome log R ratio file"),
    make_option(c("-t", "--type-of-image"), type="character", default="png",
                dest="outType", help="eps|png|svg|pdf"),
    make_option(c("-o", "--output-prefiex"), type="character", default=FALSE,
                dest="prefix", help="header of output file")
)
parser <- OptionParser(usage = "%prog [options] Germline-BAF-file", option_list=option_list)
opt <- parse_args(parser, positional_arguments=1)
# Check
if ( length(opt$args) < 1 ) {
   print("X chromosome BAF file is required.")
   q(status=1)
}
# Setting
if ( opt$options$prefix == FALSE ) {
   opt$options$prefix <- sub("_BAF","",file_path_sans_ext(basename(opt$args[1])))
}
#print(opt)
#print(paste(opt$args[1],opt$options$lrr,opt$options$prefix,opt$options$outType))
# Reading data
baf <- read.table( opt$args[1], header=T, row.names=1, comment.char="", sep = "\t", check.names=F )
# baf <- read.table( "LM_BAF.txt", header=T, row.names=1, comment.char="", sep = "\t", check.names=F )
# baf <- read.table( "UM_BAF.txt", header=T, row.names=1, comment.char="", sep = "\t", check.names=F )
# head(colnames(baf))
# [1] "chr"        "pos"        "GSM2693270" "GSM2693271" "GSM2693272"
# [6] "GSM2693273"
bafX <- baf[baf[,1]=="X",]
propHeteroX <- apply(bafX[,3:ncol(bafX)], 2, function(x) sum(0.15<x & x<0.85,na.rm=T)/sum(!is.na(x),na.rm=T))
# Estimate gender
if ( opt$options$lrr == FALSE ) {
   #genderClust <- kmeans(propHeteroX,2)
   gender <- propHeteroX
   gender[propHeteroX < 0.05] <- "XY"
   gender[propHeteroX >= 0.05] <- "XX"
} else {
   # lrr <- read.table( "LM_LogR.txt", header=T, row.names=1, comment.char="", sep = "\t", check.names=F )
   # lrr <- read.table( "UM_LogR.txt", header=T, row.names=1, comment.char="", sep = "\t", check.names=F )
   lrr <- read.table( opt$options$lrr, header=T, row.names=1, comment.char="", sep = "\t", check.names=F )
   lrrY <- lrr[lrr[,1]=="Y",]
   #meanLrrY <- apply(lrrY[,3:ncol(lrr)], 2, mean, na.rm=T )
   medianLrrY <- apply(lrrY[,3:ncol(lrr)], 2, median, na.rm=T )
   # algorithm=c("Hartigan-Wong","Lloyd","MacQueen")
   clustData <- cbind(scale(propHeteroX,center=F),scale(medianLrrY,center=F))
   genderClust <- kmeans(clustData,2,nstart=25)
   #plot(cbind(propHeteroX,medianLrrY),col=genderClust$cluster)
   gender <- genderClust$cluster
   if ( genderClust$centers[1,1] > genderClust$centers[2,1] ) {
     gender[gender==1] <- "XX"
     gender[gender==2] <- "XY"
   } else {
     gender[gender==1] <- "XY"
     gender[gender==2] <- "XX"
   }
}
outFile <- paste0(opt$options$prefix, ".gender.txt")
#y <- cbind(gender,propHeteroX)
write.table( cbind(gender,propHeteroX), file=outFile, quote=FALSE, sep="\t", col.names=NA, row.names=T )
#
if ( opt$options$lrr == FALSE ) {
   histPropX( propHeteroX, gender, opt$options$prefix, opt$options$outType )
} else {
  plotPropXY( propHeteroX, medianLrrY, gender, "Median of LogR on Y", opt$options$prefix, opt$options$outType )
}
