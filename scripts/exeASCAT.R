#!/usr/bin/env Rscript
#
# ASCAT tumor only mode
#
# ASCATtumorOnly.R GC.txt Illumina2.5M Tumor_LogR.txt Tumor_BAF.txt
#
# AffySNP6 Custom10k
# Affy10k Affy100k Affy250k_sty Affy250k_nsp AffyOncoScan AffyCytoScanHD
# HumanCNV370quad HumanCore12 HumanCoreExome24 HumanOmniExpress12
# Illumina109k IlluminaCytoSNP Illumina610k Illumina660k Illumina700k
# Illumina1M Illumina2.5M IlluminaOmni5 IlluminaOmniExpressExome
#
#library(MASS)
library(ASCAT)
library(optparse)
library(tools)
#
# Reading data
#
option_list <- list(
    make_option(c("-l", "--germline-logR-file"), type="character",default=FALSE,
                dest="Germline_LogR_file", help="Germline Log R Ratio file"),
    make_option(c("-b", "--germline-BAF-file"), type="character", default=FALSE,
                dest="Germline_BAF_file", help="Germline B Allele Freq file"),
    make_option(c("-c", "--gc-file"), type="character", default=FALSE,
                dest="GC_file", help="ASCAT GC content file"),
    make_option(c("-n", "--gender-file"), type="character", default=FALSE,
                dest="Gender_file", help="ASCAT gender file"),
    make_option(c("-p", "--platform "), type="character", default=FALSE,
                dest="platform", help="SNPs array platform: AffySNP6|Custom10k|Illumina109k|IlluminaCytoSNP|Illumina610k|Illumina660k|Illumina700k|Illumina1M|Illumina2.5M|IlluminaOmni5|Affy10k|Affy100k|Affy250k_sty|Affy250k_nsp|AffyOncoScan|AffyCytoScanHD|HumanCNV370quad|HumanCore12|HumanCoreExome24|HumanOmniExpress12|IlluminaOmniExpressExome"),
    make_option(c("-m", "--multiple-samples"), action="store_true",default=FALSE,
                dest="multi", help="Use asmultipcf segmentation"),
    make_option(c("-t", "--type-of-image"), type="character", default="png",
                dest="outType", help="eps|png|svg|pdf"),
    make_option(c("-o", "--output-prefiex"), type="character", default=FALSE,
                dest="prefix", help="header of output file")
)
parser <- OptionParser(usage = "%prog [options] Tumor-LogR-file Tumor-BAF-file", option_list=option_list)
opt <- parse_args(parser, positional_arguments=2)
# Check
if ( length(opt$args) < 2 ) {
   print("Input Tumor_LogR_file and Tumor_BAF_file.")
   q(status=1)
}
if ( (opt$options$Germline_LogR_file == FALSE | opt$options$Germline_BAF_file == FALSE)
     & opt$options$platform == FALSE ) {
   print("Please specify Germline_LogR_file and Germline_BAF_file, or platform. Genotypes are required.")
   q(status=1)
}
# Setting
Tumor_LogR_file <- opt$args[1]
Tumor_BAF_file <- opt$args[2]
outPrefix <- opt$options$prefix
if ( opt$options$prefix == FALSE ) {
   outPrefix <- sub("_BAF","",file_path_sans_ext(basename(opt$args[2])))
}
# Reading data
gender <- NULL
if ( opt$options$Gender_file != FALSE ) {
   print("Reading gender file.")
   genderData <- read.table(opt$options$Gender_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
   gender <- genderData$gender
}
if ( opt$options$Germline_LogR_file != FALSE & opt$options$Germline_BAF_file != FALSE ) {
   ascat.bc = ascat.loadData(Tumor_LogR_file,Tumor_BAF_file,opt$options$Germline_LogR_file,opt$options$Germline_BAF_file, gender=gender)
} else {
   ascat.bc = ascat.loadData(Tumor_LogR_file,Tumor_BAF_file, gender=gender)
}
# GC correction
if ( opt$options$GC_file != FALSE ) {
   ascat.bc <- ascat.GCcorrect(ascat.bc, opt$options$GC_file )
} else {
   print("Warning. GC correction is skipped.")
}
# Plot raw data
ascat.plotRawData(ascat.bc)
# Germline genotype
gg <- NULL
if ( (opt$options$Germline_LogR_file == FALSE | opt$options$Germline_BAF_file == FALSE)
     & opt$options$platform != FALSE ) {
   print("Estimating germline genotypes")
   gg <- ascat.predictGermlineGenotypes(ascat.bc, platform=opt$options$platform)
}
# Segmentation
if ( opt$options$multi ) {
  print("Mutiple sample segmentation method is selected")
  ascat.bc = ascat.asmultipcf(ascat.bc, ascat.gg=gg)
} else {
  print("Segmentation")
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)
}
# Plot segmented data
ascat.plotSegmentedData(ascat.bc)
# Run ASCAT
print("Run ASCAT")
ascat.output = ascat.runAscat(ascat.bc)
# Outputs
print("Outputing")
outName <- paste( outPrefix, ".RData", sep="" )
save( ascat.output, file=outName )
#ascat.output$nA                    ascat.output$failedarrays
#ascat.output$nB                    ascat.output$nonaberrantarrays
#ascat.output$aberrantcellfraction  ascat.output$segments
#ascat.output$ploidy                ascat.output$segments_raw
#ascat.output$psi                   ascat.output$distance_matrix
#ascat.output$goodnessOfFit
#
if ( length(ascat.output$failedarrays) > 0 ) {
   outName <- paste( outPrefix, ".failure", sep="" )
   write.table( ascat.output$failedarrays, file=outName, quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE )
}
#
if ( length(ascat.output$nonaberrantarrays) > 0 ) {
   outName <- paste( outPrefix, ".noaberrant", sep="" )
   write.table( ascat.output$nonaberrantarrays, file=outName, quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE )
}
#
outName <- paste( outPrefix, ".segments", sep="" )
write.table( ascat.output$segments, file=outName, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE )
#
# aberrantcellfraction and psi
#
gP <- cbind( ascat.output$aberrantcellfraction, ascat.output$ploidy, ascat.output$psi )
rownames(gP) <- colnames(ascat.output$nA)
colnames(gP) <- c("aberr","ploidy","psi")
outName <- paste( outPrefix, ".aberr", sep="" )
write.table( gP, file=outName, quote=FALSE, sep="\t", col.names=NA, row.names=TRUE )
q(status=0)
#
# LiverMetastasis
#
genderData <- read.table("LM.gender.txt", header=T, row.names=1)
gender <- genderData$gender
ascat.bc = ascat.loadData("LM_Tumor_LogR.txt","LM_Tumor_BAF.txt", gender=gender)
ascat.bc <- ascat.GCcorrect(ascat.bc, "GC_Omni2.5pos.txt")
ascat.plotRawData(ascat.bc)
gg <- ascat.predictGermlineGenotypes(ascat.bc, platform="Illumina2.5M")
ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)
#
# UvealMelanoma
#
genderData <- read.table("UM.gender.txt", header=T, row.names=1)
gender <- genderData$gender
ascat.bc = ascat.loadData("UM_Tumor_LogR.txt","UM_Tumor_BAF.txt","UM_Germline_LogR.txt","UM_Germline_BAF.txt")
ascat.bc <- ascat.GCcorrect(ascat.bc, "GC_Illumina660k.txt" )
ascat.plotRawData(ascat.bc)
#gg <- ascat.predictGermlineGenotypes(ascat.bc, platform="Illumina660k")
ascat.bc = ascat.aspcf(ascat.bc) #, ascat.gg=gg
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)
