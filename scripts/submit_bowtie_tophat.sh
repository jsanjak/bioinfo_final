#!/bin/sh
cd data
cwd=`pwd`
for exp in RNAseq
do
    cd $exp
    for tiss in Body Embryo Head Pupa
    do
        cd $tiss
        ls *fq.gz | rev | cut -d. -f4- | rev | uniq | while read read_base
        do
            folder=$cwd/$exp/$tiss
            echo $folder 
            echo $read_base
            qsub -N bt2.align.${exp}.${tiss}.${read_base} ../../../run_bowtie_tophat.sh $folder $read_base 
        done
        cd ..
    done
    cd ..
done
cd ..

