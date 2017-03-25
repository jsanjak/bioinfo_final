#!/bin/sh
#$ -q krt,krti,bio,pub64
#$ -N build_indexes

module load bwa/0.7.12
module load samtools/1.3
module load bowtie2/2.2.7

cd ref

ref="dmel-all-chromosome-r6.13.fasta"
prefix=$(echo $ref | cut -d. -f1)
bwa index -p $prefix $ref
samtools faidx  $ref
java -d64 -Xmx128g -jar /data/apps/picard-tools/1.87/CreateSequenceDictionary.jar  R=$ref O=${prefix}.dict
bowtie2-build  $ref $prefix

