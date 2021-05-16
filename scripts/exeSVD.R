#!/usr/bin/env Rscript
#
# Principal component analysis
#
# % exeSVD.R columeSample-rowProbe.txt
# % R CMD BATCH .R
# #source(".R")
#
library(optparse)
library(tools)
library(svd)
#library(data.table)
options(warn=2)
#library(MASS)
#
# Function
#
writePCAresult <- function( data, sFlag, noe, prefix, outType ) {
  if ( is.na(noe) ) {
    noe <- min(dim(data))
  }
  if ( sFlag == "TRUE" ) {
    data <- scale(data, center=T, scale=T)
    #data <- t(scale(t(data), center=T, scale=T))
  } else {
    data <- scale(data, center=T, scale=F)
    #data <- t(scale(t(data), center=T, scale=F))
  }
  pcaRes <- propack.svd( data, neig=noe, opts=list(kmax=500) )
  outName = paste0( prefix, ".d" )
  write.table( pcaRes$d, file=outName, quote=FALSE, sep="\t" )
  outName = paste0( prefix, ".u" )
  rownames(pcaRes$u) <- rownames(data)
  write.table( pcaRes$u, file=outName, quote=FALSE, sep="\t" )
  outName = paste0( prefix, ".v" )
  rownames(pcaRes$v) <- colnames(data)
  write.table( pcaRes$v, file=outName, quote=FALSE, sep="\t" )
  # screeplot
  outName <- paste0( prefix, ".d.", outType )
  if ( outType == "svg" ) {
    svg( outName, height=6, width=6 )
  } else if ( outType == "png" ) {
    png( outName, height=600, width=600 )
  } else if ( outType == "pdf" ) {
    pdf( outName, height=6, width=6 )
  } else {
    postscript( outName, paper='special', height=6, width=6, horizontal=F )
  }
  barplot( pcaRes$d, main="Eigenvalues", xlab="Principal components")
  axis( 1, at=1:noe, labels=1:noe, tick=F )
  dev.off()
}
# Old
writePrcompResult <- function( data, sFlag, prefix ) {
  pcaRes <- prcomp( data, scale=sFlag )
  outName = paste0( prefix, "Eigen.txt" )
  write.table( pcaRes$sdev[1:8], file=outName, quote=FALSE, sep="\t" )
  outName = paste0( prefix, "Score.txt" )
  write.table( pcaRes$x[,1:8], file=outName, quote=FALSE, sep="\t" )
  outName = paste0( prefix, "Rotation.txt" )
  write.table( pcaRes$rotation[,1:8], file=outName, quote=FALSE, sep="\t" )
  # screeplot
  outName = paste0( prefix, "Screeplot.eps" )
  postscript(outName, paper='special', height=6, width=6, horizontal=F)
  screeplot( pcaRes, main="Eigenvalues", xlab="Principal components")
  dev.off()
}
#
# Reading data
#
option_list <- list(
    make_option(c("-n", "--number-of-eigenvalues"), type="numeric", default=NA,
                dest="noe", help="number of desired eigentriples"),
    make_option(c("-t", "--type-of-image"), type="character", default="png",
                dest="outType", help="eps|png|svg|pdf"),
    make_option(c("-o", "--output-prefiex"), type="character", default=FALSE,
                dest="prefix", help="header of output file")
)
parser <- OptionParser(usage = "%prog [options] data-file", option_list=option_list)
opt <- parse_args(parser, positional_arguments=1)
# Check
if ( length(opt$args) < 1 ) {
   print("Input data file is required.")
   q(status=1)
}
# Setting
if ( opt$options$prefix == FALSE ) {
  opt$options$prefix <- sub("_LogR","",file_path_sans_ext(basename(opt$args[1])))
}
#print(opt)
#print(paste(opt$args[1],opt$options$prefix,opt$options$outType))
data <- read.table( opt$args[1], header=T, row.names=1, comment.char="", sep = "\t", check.names=F )
#data <- read.table( gzfile("../LRR-BAF/LM_Tumor_LogR.txt.gz"), header=TRUE )
#data <- read.table("../LRR-BAF/LM_Tumor_LogR.txt",header=T, row.names=1, comment.char="", sep = "\t", check.names=F )
#
#> head(colnames(data))
#[1] "chr"   "pos"   "LM16M" "LM16T" "LM7M"  "LM7T" 
chrs = paste(c(1:22))
data <- data[data[,1]%in%chrs,]
data <- data[complete.cases(data),3:ncol(data)]
writePCAresult( t(data), FALSE, opt$options$noe, paste(opt$options$prefix,"_OrgSample",sep=""), opt$options$outType )
writePCAresult( t(data), TRUE, opt$options$noe, paste(opt$options$prefix,"_ScaSample",sep=""), opt$options$outType )
writePCAresult( data, FALSE, opt$options$noe, paste(opt$options$prefix,"_OrgProbe",sep=""), opt$options$outType )
writePCAresult( data, TRUE, opt$options$noe, paste(opt$options$prefix,"_ScaProbe",sep=""), opt$options$outType )
#
#writePrcompResult( t(data), FALSE, paste(opt$options$prefix,"_OrgSamp",sep="") )
#writePrcompResult( t(data), TRUE,  paste(opt$options$prefix,"_ScaSamp",sep="") )
#writePrcompResult( data, FALSE, paste(opt$options$prefix,"_OrgProb",sep="") )
#writePrcompResult( data, TRUE,  paste(opt$options$prefix,"_ScaProb",sep="") )
