#!/bin/bash

bwa="/data/eg01/xdong2/apps/novoBreak/nb_distribution/bin/bwa"
ref="/data/eg03/sccaller/file_in_147/ref/human_g1k_v37_decoy.fasta"
r1=1
r2=1
if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <novoBreak_result> <working_dir> <output> <split_line_num> <ifq_dir> <sv2alignments> <sample_bam> <sample_fq1> <sample_fq2> <ghdbsnp>"
    echo "ghdbsnp: germline heterozygous dbsnp"
    exit 1
fi
novoBreak_result=$1
working_dir=$2
# sv2alignments="/data/eg03/sccaller/apps/novoBreak/hb_h24/ssake/sv2alignments"
sv2alignments=$6
output=$3
# sample_bam="/data/eg03/trafic/projects/2018-TE/bam_fibroblast/h24.gatk.bam"
sample_bam=$7
# sample_fq1="/data/eg03/sccaller/apps/novoBreak/data/hs-j258-h24-r1.fq"
sample_fq1=$8
# sample_fq2="/data/eg03/sccaller/apps/novoBreak/data/hs-j258-h24-r2.fq"
sample_fq2=$9
# ghdbsnp="/data/eg03/sccaller/apps/novoBreak/hb.ht.vcf"
ghdbsnp=${10}
split_line_num=$4
ifq_dir=$5

if [ ${working_dir: -1} == "/" ]; then
	splited_dir=${working_dir}"splited/"
else
	splited_dir=${working_dir}"/splited/"
	working_dir=${working_dir}"/"
fi	

# echo $splited_dir
# echo $working_dir
mkdir -p $working_dir
mkdir -p $splited_dir

split -l $split_line_num $novoBreak_result $splited_dir

for sv_file in $(ls $splited_dir)
do
   echo "submitting: splited/"$sv_file
   qsub qsub_sc3.sh sc3procedure $bwa $sv2alignments $ref "splited/"$sv_file $sample_bam $sample_fq1 $sample_fq2 $ghdbsnp $r1 $r2 $output $working_dir $ifq_dir
done


