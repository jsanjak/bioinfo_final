#!/bin/sh
#$ -q krt,krti,sf,bio,pub64
#$ -js 100 
#$ -N HMM_geno_subset

#This script just subsets the HMM genotype data for ease of use later.
#It subsets the data to just the RILs used in this experiment and only the additive genotype probabilities.
#RILS.txt is just a list of ril codes from the RNAseq sample coding file
#the HMM files are from the DSPR website
cd /share/kevin2/jsanjak/RNAseq/ref
awk 'NR==FNR {a[$1]; next} {for(i in a){ if( $3 == i ) {print $1 " " $2 " " $3 " " $40 " " $41 " " $42 " " $43 " " $44 " " $45 " " $46 " " $47 " " $48} }}' RILS.txt HMMregA_R2.txt  >> HMM_ALL.txt
awk 'NR==FNR {a[$1]; next} {for(i in a){ if( $3 == i ) {print $1 " " $2 " " $3 " " $40 " " $41 " " $42 " " $43 " " $44 " " $45 " " $46 " " $47 " " $48} }}' RILS.txt HMMregB_R2.txt  >> HMM_ALL.txt
