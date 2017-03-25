library( "GenomicFeatures" )
library( "Rsamtools" )
library( "GenomicAlignments" )

setwd("/share/kevin2/jsanjak/RNAseq")

#hse <- makeTxDbFromGFF( "ref/dmel-all-r6.13.gtf", format="gtf" )
#exonsByGene <- exonsBy( hse, by="gene" )
Bodyfls <- list.files( "data/RNAseq/Body", pattern="*bowtie.sort.bam$", full=TRUE, recursive=T )
#Embryofls <- list.files( "data/RNAseq/Embryo", pattern="*bowtie.sort.bam$", full=TRUE, recursive=T )
BodybamLst <- BamFileList( Bodyfls)#, yieldSize=100000 )
Bodyse <- summarizeOverlaps( exonsByGene, BodybamLst, mode="Union", singleEnd=FALSE,ignore.strand=TRUE,fragments=TRUE )

### this is a big "counts table" of raw reads counts, once you have it you don't need the bam files again
### save se so you don't have to run the above again
save(Bodyse,file="data/RNAseq/Bodyse.rda")
