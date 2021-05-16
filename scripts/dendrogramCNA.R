#!/usr/bin/env Rscript
#
# Make hierarchical clusterings
#
# dendrogramCNA.R ASCAT.matrix lungOrg -t eps -o cohortA
# % R CMD BATCH .R
# #source(".R")
#
#library(MASS)
#library(stats)
library(optparse)
library(tools)
library(mclust)
library(amap)
#library(hopach)
#
# Sample dendrogram
#
plotSampleDendrogram <- function( sample.dist, distName, titleHead, outHead, outType ) {
  #
  # Ward
  #
  outName <- paste0( outHead, substr(distName,1,4), "Ward.", outType )
  if ( outType == "svg" ) {
    svg( outName, height=12, width=18 )
  } else if ( outType == "png" ) {
    png( outName, height=1200, width=1800 )
  } else if ( outType == "pdf" ) {
    pdf( outName, height=12, width=18 )
  } else {
    postscript( outName, paper='special', height=12, width=18, horizontal=F )
  }
  par(mar=c(10.5,4,2,0))
  title <- paste( titleHead, ": ", distName, ", Ward", sep="" )
  sample.hist <- hclust(sample.dist,"ward")
  #
  outTable <-paste0( outHead, substr(distName,1,4), "Ward.order")
  ids <- sample.hist$labels[sample.hist$order]
  #ids <- sub( "\\.", "-", sample.hist$labels[sample.hist$order] )
  write.table(ids, file=outTable, quote=FALSE)
  print(sample.hist$labels[sample.hist$order])
  print(sample.hist$merge)
  print( cor(sample.dist, cophenetic(sample.hist)) )
  #
  cS <- 0.5
  if ( attr(sample.dist,"Size") <= 30 ) {
    cS <- 0.8
  } else if ( attr(sample.dist,"Size") <= 50 ) {
    cS <- 0.75
  }
  plot( sample.hist, xlab="", main=title, sub="", cex=cS, hang=-1 )
  #cutree( sample.hist, k=3)
  dev.off()
}
#
# Sample dendrogram
#
plotSampleDendroManager <- function( cnaData, titleHead, outHead, outType ) {
  #
  # Euclidian distance
  #
  sample.dist <- dist( t(cnaData) )
  plotSampleDendrogram( sample.dist, "Euclid", titleHead, outHead, outType )
#  #
#  # Cosine distance
#  #
#  sample.dist <- Dist( t(cnaData), method="pearson" )
#  plotSampleDendrogram( sample.dist, "Cosine", titleHead, outHead, outType )
#  #
#  # Manhattan distance
#  #
#  sample.dist <- dist( t(cnaData), "manhattan" )
#  plotSampleDendrogram( sample.dist, "Manhattan", titleHead, outHead, outType )
}
#
# Setting options
#
option_list <- list(
    make_option(c("-t", "--type-of-image"), type="character", default="png",
                dest="outType", help="eps|png|svg|pdf"),
    make_option(c("-o", "--output-prefiex"), type="character", default=FALSE,
                dest="prefix", help="header of output file")
)
parser <- OptionParser(usage = "%prog [options] CNA-matrix-file", option_list=option_list)
opt <- parse_args(parser, positional_arguments=1)
# Check
if ( length(opt$args) < 1 ) {
   print("Copy number matrix file must be specified.")
   q(status=1)
}
# Setting
if ( opt$options$prefix == FALSE ) {
   opt$options$prefix <- sub(".matrix","",file_path_sans_ext(basename(opt$args[1])))
}
#print(opt)
#print(paste(opt$args[1],opt$options$prefix,opt$options$outType))
#
# Reading data
#
#cna <- read.table( "LM_Tumor.matrix", header=TRUE )
cna <- read.table( opt$args[1], header=TRUE )
#
# Make data
#
numS <- length( cna[1,] )
rownames(cna) <- paste( "chr", cna$chr, ":", cna$start, "-", cna$end, sep="" )
cnaData <- cna[,c(4:numS)]
#
# Calculation
#
titleHead <- paste( opt$options$prefix, "dendrogram" )
outHead <- paste0( opt$options$prefix, "Dendrogram" )
outHead <- gsub( " ", "", outHead )
plotSampleDendroManager( cnaData, titleHead, outHead, opt$options$outType )
