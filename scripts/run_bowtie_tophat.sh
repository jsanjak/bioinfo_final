#!/bin/sh

#$ -q krt,krti,bio,pub64
#$ -pe openmp 8

module load bwa/0.7.12
module load samtools/1.3
module load bowtie2/2.2.7
module load tophat/2.1.0

cd ref 

folder=$1
read_base=$2
ref="dmel-all-chromosome-r6.13.fasta"
ref_base=$(echo $ref | cut -d. -f1)


mkdir -p $folder/alignments/$read_base

multimapping=4    # number of reads reported for multimapping
bowtie2 -k $multimapping -X2000 --mm --threads 8 -x $ref_base -1 $folder/${read_base}.F.fq.gz  -2 $folder/${read_base}.R.fq.gz | samtools view -bS - > $folder/alignments/${read_base}/${read_base}.bowtie.bam

samtools sort $folder/alignments/${read_base}/${read_base}.bowtie.bam -o $folder/alignments/${read_base}/${read_base}.bowtie.sort.bam
samtools index $folder/alignments/${read_base}/${read_base}.bowtie.sort.bam

tophat -p 8 -G dmel-all-r6.13.gtf -o $folder/alignments/${read_base} $ref_base  $folder/${read_base}.F.fq.gz  $folder/${read_base}.F.fq.gz
samtools sort  $folder/alignments/${read_base}/accepted_hits.bam -o  $folder/alignments/${read_base}/accepted_hits.sort.bam

