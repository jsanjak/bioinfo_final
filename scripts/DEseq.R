args=(commandArgs(TRUE))
library( "GenomicFeatures" )
library( "Rsamtools" )
library( "GenomicAlignments" )
library( "DESeq2" )
library(dplyr)
library(data.table)
library(MASS)
library(lme4)

#setwd("/share/kevin2/jsanjak/RNAseq")

gene_j<-args[1]
chrom<-args[2]
tiss<-args[3]
load(paste("data/RNAseq/",tiss,"se.rda",sep=""))

genes <- rownames(assays(tiss.se)$counts)
RILgenos <- fread("ref/HMM_ALL.txt")
colnames(RILgenos) <- c("chr","pos","RILcode","B1","B2","B3","B4","B5","B6","B7","AB8")
load('ref/poslist.rda')

ofile<-paste("data/eqtl/",tiss,"/",gene_j,".",chrom,".txt",sep="")
write("gene tissue chr pos LRstat logP",ofile,ncol=6)
for (i in seq(nrow(pos_chrom))){
    print(c(j,i))
    tiss.se_pos <- tiss.se
    #get the genotypes at this position
    geno_pos <- filter(RILgenos,chr==chrom,pos==pos_chrom[i,2])
    #Tack the genotypes onto the sample info
    sampleInfo_pos <- left_join(sampleInfo,dplyr::select(geno_pos,-chr,-pos),by="RILcode")

    sampleInfo_pos <- do.call(rbind,lapply(c("Flowcell01","Flowcell02","Flowcell03","Flowcell04","Flowcell05"),function(x)mutate(sampleInfo_pos,Flowcell=x)))

    #get the runs that actually happened
    sampleInfo_pos <- sampleInfo_pos %>% mutate(run= paste(Flowcell,".Sample",SampleNumber,".RIL.",RILcode,".",i7index,".bowtie.sort.bam",sep=""))
    seIdx <- match(colnames(tiss.se_pos), sampleInfo_pos$run)
    colData(tiss.se_pos) <- cbind( colData(tiss.se_pos), sampleInfo_pos[ na.omit(seIdx), ] )

    #Note that we do not include B1 genotype as this will be our "baseline"
    fulldata.gene <-data.frame(ge= as.numeric(assays(tiss.se_pos)$counts[gene_j,]),dplyr::select(as.data.table(colData(tiss.se_pos)),Flowcell,RILcode,B2,B3,B4,B5,B6,B7,AB8))

    mnull <- glm.nb("ge ~ 1", data = fulldata.gene)
    m1 <- glm.nb("ge ~ B2 + B3 + B4 + B5 + B6 + B7 + AB8" , data = fulldata.gene)

    #lmmnull <- glmer.nb("ge ~ (1 | Flowcell)",data = fulldata.gene)
    #lmm1 <- glmer.nb("ge ~  B2 + B3 + B4 + B5 + B6 + B7 + AB8 + (1 | Flowcell)",data = fulldata.gene)
    anova(mnull, m1, test="LRT")

    LRstat <- anova(mnull,m1)[["LR stat."]][2]
    logP <- -pchisq(LRstat,7,lower.tail=F,log.p=T)/(log(10))
    output <-c(gene_j,tiss,chrom,pos_chrom[i,2],LRstat,logP)
    write(output,ofile,ncol=length(output),append=T)
}
