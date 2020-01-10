#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
# Job name
#$ -N sc3
# Number of cpu cores required
#$ -pe smp 1
# RAM requirement per cpu core
#$ -l h_vmem=16G
# require computer to handle
# -l h=EG01
# Email
#$ -M spsc83@gmail.com
#$ -m bes

export LC_ALL=C
export MALLOC_ARENA_MAX=4


echo "==START=="; date

echo "python sc3.py ${1} ${2} ${3} ${4} ${5} $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18}"
python sc3.py ${1} ${2} ${3} ${4} ${5} $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18}
date; echo "==END=="
