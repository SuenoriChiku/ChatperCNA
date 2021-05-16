#!/usr/bin/env Rscript
#
# Plot PCA result
# plotPCAsatter.R sample-information.txt SVD.u|v
# % R CMD BATCH .R
# #source(".R")
#
#library(MASS)
library(optparse)
library(tools)
#
# Sample
#
plotSampleID <- function( xData, yData, sampleID, titlePrefix, xLab, yLab, outPrefix, outType ) {
  title <- paste(titlePrefix, "SampleID")
  xmin = min( xData, na.rm=TRUE )
  xmax = max( xData, na.rm=TRUE )
  ymin = min( yData, na.rm=TRUE )
  ymax = max( yData, na.rm=TRUE )
  #
  outName <- paste0( outPrefix, ".SampleID.", outType )
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
  text( xData, yData, sampleID, cex=0.6 )
  dev.off()
}
#
plotScatter <- function( xData, yData, pheno, phenoName, titlePrefix, xLab, yLab, outPrefix, outType ) {
  xmin = min( xData, na.rm=TRUE )
  xmax = max( xData, na.rm=TRUE )
  ymin = min( yData, na.rm=TRUE )
  ymax = max( yData, na.rm=TRUE )
  fp <- factor(pheno)
  cT <- levels(fp)
  nL <- nlevels(fp)
  cR <- rainbow(nL+5)[1:nL]
  nP <- 20
  fN <- floor(nL/nP)
  cP <- c(rep(1:nP,fN),1:(nL-fN*nP))
  cS <- 0.5
  yM <- 0.25
  title <- paste0( titlePrefix, ": phenoName" )
  outName <- paste( outPrefix, phenoName, outType, sep="." )
  if ( outType == "eps" ) {
    postscript( outName, paper='special', height=6, width=6, horizontal=F)
  } else if ( outType == "png" ) {
    png( outName, height=600, width=600 )
    cS <- 1.0
    yM <- 0.1
  } else if ( outType == "svg" ) {
    svg( outName, height=6, width=6 )
  } else {
    pdf( outName, height=6, width=6 )
  }
  plot(NA,NA,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xLab,ylab=yLab,main=title)
  for( i in 1:nL ) {
    points( xData[fp==cT[i]], yData[fp==cT[i]],col=cR[i], pch=cP[i], cex=cS )
  }
  par(xpd=TRUE)
  legend( xmax+(xmax-xmin)*0.03, ymax+(ymax-ymin)*yM, cT, col=cR, pch=cP, cex=0.5, bty="n")
  dev.off()
}
#
# Reading data
#
option_list <- list(
    make_option(c("-m", "--main-tile-prefix"), type="character", default="PCA scatter plot",
                dest="tilePrefix", help="Tile prefix of th figure"),
    make_option(c("-t", "--type-of-image"), type="character", default="png",
                dest="outType", help="eps|png|svg|pdf"),
    make_option(c("-o", "--output-prefiex"), type="character", default=FALSE,
                dest="prefix", help="header of output file")
)
parser <- OptionParser(usage = "%prog [options] experiments-information-file SVD.u|v-file", option_list=option_list)
opt <- parse_args(parser, positional_arguments=2)
# Check
if ( length(opt$args) < 2 ) {
   print("Information fila and PCA result file are required.")
   q(status=1)
}
# Setting
if ( opt$options$prefix == FALSE ) {
  opt$options$prefix <- sub(".[u|v]$","",file_path_sans_ext(basename(opt$args[2])))
}
#print(opt)
#print(paste(opt$args[1],opt$options$prefix,opt$options$outType))
info <- read.table( opt$args[1], header=T, sep = "\t" )
results <- read.table( opt$args[2], header=T, sep = "\t" )
#info <- read.table( "../misc/GSE100787_sampleInfo.txt", header=T, sep = "\t" )
#results <- read.table( "LM_Tumor_OrgSample.u", header=T, sep = "\t" )
#> colnames(info)
#[1] "Accession" "SampleID"  "tissue"    "chip"
#> colnames(results)
#[1] "V1" "V2" "V3" "V4"
info <- info[info[,2] %in% rownames(results),]
results <- results[rownames(results) %in% info[,2], ]
#
cLabels <- colnames(info)
xLabH <- "Score"
yLabH <- "Score"
xLab <- paste( xLabH, 1 )
yLab <- paste( yLabH, 2 )
outPrefix12 <- paste0( opt$options$prefix, 12 )
plotSampleID( results[,1], results[,2], info[,2], opt$options$tilePrefix, xLab, yLab, outPrefix12, opt$options$outType )
for( i in 3:ncol(info)) {
   plotScatter( results[,1], results[,2], info[,i], cLabels[i], opt$options$tilePrefix, xLab, yLab, outPrefix12, opt$options$outType )
}
#
xLab <- paste( xLabH, 3 )
yLab <- paste( yLabH, 4 )
outPrefix34 <- paste0( opt$options$prefix, 34 )
plotSampleID( results[,3], results[,4], info[,2], opt$options$tilePrefix, xLab, yLab, outPrefix34, opt$options$outType )
for( i in 3:ncol(info)) {
   plotScatter( results[,3], results[,4], info[,i], cLabels[i], opt$options$tilePrefix, xLab, yLab, outPrefix34, opt$options$outType )
}
