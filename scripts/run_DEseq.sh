#!/bin/sh
#$ -q krt,krti,bio,pub64

module load R
R --vanilla --slave --args < DEseq.R $1 $2 $3
