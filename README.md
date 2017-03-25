## bioinfo_final
# Tissue specific eQTL

For the final project, I decided to try to develop a tissue specific eQTL mapping pipeline using the RIL RNAseq data. I pulled the [RIL data](http://wfitch.bio.uci.edu/~dspr/Data/index.html) 
from Tony Long's website. See A and B RIL HMM genotypes can be found in the ref folder.

# Pipline
* Organize data see: scripts/ADL_RNASEQ_structure.py
* Build reference genome indices see scripts/build_indices.sh
* Align data with bowtie/tophat see: scripts/ run_bowtie_tophat and submit_bowtie_tophat
* Make exon count tables using GFF file see: scripts/make_ExonSummary.R and submit_ExonSummary.sh
* Extract specific genotypes of RILS used in the experiment: scripts/extract_RILS.sh
* Loop over genes, genomics position and tissues and fit a log-negative binomial generalized linear regression
model to the data using the founder haplotypes as the explanatory variables. Then analyze significance via likelihood ratio test
against the null model of no genetic effect. see scripts/DEseq.R

# Todo
* Finish analyzing all the tissues
* Account for genomic inflation of test statistic
* Consider fitting the flowcell as a randome effect covariate
* Compare eQTL peaks across tissues and see which ones are tissue specific.
* Consider calculating the genome-wide genetic correlation of gene expression between tissues by gene
    * This would involve calculating the genetic relatedness matrix and fitting the multivariate LMM to estimate variance components


