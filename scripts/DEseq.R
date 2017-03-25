args=(commandArgs(TRUE))
library( "GenomicFeatures" )
library( "Rsamtools" )
library( "GenomicAlignments" )
library( "DESeq2" )
library(dplyr)
library(data.table)
library(MASS)
library(lme4)

setwd("/share/kevin2/jsanjak/RNAseq")

#gene_j<-args[1]

chrom<-"3R"#args[2]
load("data/RNAseq/Bodyse.rda")
genes <- rownames(assays(Bodyse)$counts)
RILgenos <- fread("ref/HMM_ALL.txt")
colnames(RILgenos) <- c("chr","pos","RILcode","B1","B2","B3","B4","B5","B6","B7","AB8")
load('ref/poslist.rda')

for ( j in seq(100) ){
  gene_j <- genes[j]
  ofile<-paste("data/eqtl/Body.",gene_j,".",chrom,".txt",sep="")#)#args[3]
  write("gene tissue chr pos LRstat logP",ofile,ncol=6)
  for (i in seq(nrow(pos_chrom))){
    print(c(j,i))
    Bodyse_pos <- Bodyse
    #get the genotypes at this position
    geno_pos <- filter(RILgenos,chr==chrom,pos==pos_chrom[i,2])
    #Tack the genotypes onto the sample info
    sampleInfo_pos <- left_join(sampleInfo,dplyr::select(geno_pos,-chr,-pos),by="RILcode")
    
    sampleInfo_pos <- do.call(rbind,lapply(c("Flowcell01","Flowcell02","Flowcell03","Flowcell04","Flowcell05"),function(x)mutate(sampleInfo_pos,Flowcell=x)))
    
    #get the runs that actually happened
    sampleInfo_pos <- sampleInfo_pos %>% mutate(run= paste(Flowcell,".Sample",SampleNumber,".RIL.",RILcode,".",i7index,".bowtie.sort.bam",sep=""))
    seIdx <- match(colnames(Bodyse_pos), sampleInfo_pos$run)
    colData(Bodyse_pos) <- cbind( colData(Bodyse_pos), sampleInfo_pos[ na.omit(seIdx), ] )
    
    #Note that we do not include B1 genotype as this will be our "baseline"
    fulldata.gene <-data.frame(ge= as.numeric(assays(Bodyse_pos)$counts[gene_j,]),dplyr::select(as.data.table(colData(Bodyse_pos)),Flowcell,RILcode,B2,B3,B4,B5,B6,B7,AB8))
    
    mnull <- glm.nb("ge ~ 1", data = fulldata.gene)
    m1 <- glm.nb("ge ~ B2 + B3 + B4 + B5 + B6 + B7 + AB8" , data = fulldata.gene)
    
    #lmmnull <- glmer.nb("ge ~ (1 | Flowcell)",data = fulldata.gene)
    #lmm1 <- glmer.nb("ge ~  B2 + B3 + B4 + B5 + B6 + B7 + AB8 + (1 | Flowcell)",data = fulldata.gene)
    anova(mnull, m1, test="LRT")
    
    LRstat <- anova(mnull,m1)[["LR stat."]][2]
    logP <- -pchisq(LRstat,7,lower.tail=F,log.p=T)/(log(10))
    #Logb x = Loga x/Loga b
    #LOD_m <-  (-2*(logLik(mnull) - logLik(m1)))/(2*log(10))
    #LOD_lmm <-  (-2*(logLik(lmmnull) - logLik(lmm1)))/(2*log(10))
    output <-c(gene_j,"Body",chrom,pos_chrom[i,2],LRstat,logP)
    write(output,ofile,ncol=length(output),append=T)
  }
}
#0 + B1 + B2 + B3 + B4 + B5 + B6 + B7 + AB8
####
#ddsFull <- DESeqDataSet( Bodyse, design = ~ Flowcell )    # names must match sampleInfo colnames

#dds <- DESeq(ddsFull)

#res_PB <- results(dds, contrast=c("TissueCode","P","B"))
#res_PE <- results(dds, contrast=c("TissueCode","P","E"))
#res_PH <- results(dds, contrast=c("TissueCode","P","H"))
#res_HB <- results(dds, contrast=c("TissueCode","H","B"))
#res_HE <- results(dds, contrast=c("TissueCode","H","E"))
#res_BE <- results(dds, contrast=c("TissueCode","B","E"))



