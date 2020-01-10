#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
# Job name
#$ -N nb
# Number of cpu cores required
#$ -pe smp 10
# RAM requirement per cpu core
#$ -l h_vmem=16G
# require computer to handle
# -l h=EG01
# Email
#$ -M spsc83@gmail.com
#$ -m bes

export LC_ALL=C
export MALLOC_ARENA_MAX=4
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 data_root normal tumor" >&2
    exit 1
fi

nb=/data/eg01/xdong2/apps/novoBreak/nb_distribution
ref=~/projects/sccaller/ref/human_g1k_v37_decoy.fasta
np=8

dat_root=$1
tumor=$3
normal=$2
my_dir=$(pwd)
echo $normal
echo $tumor

echo "==START=="; date
echo "my_dir="$my_dir
./run_novoBreak.sh ${nb} ${ref} ${dat_root}/${tumor}.gatk.bam ${dat_root}/${normal}.gatk.bam ${np} ${normal}_${tumor}

date; echo "==END=="
