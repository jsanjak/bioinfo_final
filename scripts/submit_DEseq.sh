#!/bin/sh

while read gene
do 
    for chrom in "3R" # "2L" "2R" "3R" "3L" "X" "4"
    do
        for tiss in "Body" #"Embryo" "Head" "Pupa"
        do
            qsub -N DEseq.${gene}.${chrom}.${tiss} run_DEseq.sh $gene $chrom $tiss 
        done
    done    
done < ref/FBgenes_test.txt #FBgenes.txt
