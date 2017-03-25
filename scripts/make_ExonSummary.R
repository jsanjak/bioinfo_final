args=(commandArgs(TRUE))
library( "GenomicFeatures" )
library( "Rsamtools" )
library( "GenomicAlignments" )

tiss<-args[1]

hse <- makeTxDbFromGFF( "ref/dmel-all-r6.13.gtf", format="gtf" )
exonsByGene <- exonsBy( hse, by="gene" )
tiss.fls <- list.files( paste("data/RNAseq/",tiss,sep=""), pattern="*bowtie.sort.bam$", full=TRUE, recursive=T )
tiss.bamLst <- BamFileList( tiss.fls)#, yieldSize=100000 )
tiss.se <- summarizeOverlaps( exonsByGene, tiss.bamLst, mode="Union", singleEnd=FALSE,ignore.strand=TRUE,fragments=TRUE )

### this is a big "counts table" of raw reads counts, once you have it you don't need the bam files again
### save se so you don't have to run the above again
save(tiss.se,file=paste("data/RNAseq/",tiss,"se.rda",sep=""))
