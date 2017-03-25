#!/bin/sh
#$ -q krt,krti,bio,pub64

module load R

for tiss in "Body" #"Embryo" "Head" "Pupa"
do
R --vanilla --slave --args < scripts/make_ExonSummary.R $tiss
done
