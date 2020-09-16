PEA: Identify genome structural variations SVs from single-cell whole-genome sequencing data
=
Version 3.1.0

Updated date: 2020.08.01

Cite us:

Dong X et al. Identifying genome structural variations in single cells. In submission.

#####
Author and License
-
Authors: Xiao Dong, Yujue Wang

Email: biosinodx@gmail.com (X.D.), xiao.dong@einsteinmed.org (X.D.), spsc83@gmail.com (Y.W.), yujue.wang@einsteinmed.org(Y.W.)

Licensed under the GNU Affero General Public License version 3 or later


#####
Dependence
-
Linux OS (CentOS7 tested)

Python3 (version 3.6.9 tested)
Python packages: argparse, os, re, argpars, sys, string, scipy

R (version 3.6.0 tested)
R packages: VariantAnnotation, StructuralVariantAnnotation, stringr, Biostrings, BSgenome with BSgenome.Hsapiens.UCSC.hg38 (for hg38 genome)

Manta (version 1.6.0 tested)
https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

Samtools (version 1.9.0 tested)

BWA (version 0.7.17 tested)

#####
Installation guide
-
The scripts can be used directly without installing.


# Input data requirement

1. bam file(s): for germline SV calling, require bam file of the sample; for somatic SV calling, require bam files of both a single cell and bulk DNA of the cell's origin

2. hSNP in vcf format

3. reference genome in fasta format


# USAGE STEP 1. Calling candidate breakpoint ends (BPEs) and assembly of their contigs

This step is done with using Manta, but it can be replaced by other software tools. 

1.1 For germline SVs, run the following shell scripts in Linux,

mkdir -p ${sn}

${manta_dir}/bin/configManta.py \
  --normalBam=${bam} \
  --referenceFasta=${ref} \
  --runDir=${sn}

$(pwd)/${sn}/runWorkflow.py

In the above,

${manta_dir} is the path to the directory of Manta software.

${bam} is the path to sample bam file.

${ref} is the path to reference genome fasta file.

${sn} is the sample name.

1.2 For somatic SVs, i.e., SVs in a single cell but not in bulk DNA of its origin, run the following shell scripts in Linux,

mkdir -p ${bulk}.${cell}

${manta_dir}/bin/configManta.py \
  --normalBam=${bulk_bam} \
  --tumorBam=${cell_bam} \
  --referenceFasta=${ref} \
  --runDir=${bulk}.${cell}

$(pwd)/${bulk}.${cell}/runWorkflow.py

In the above,

${manta_dir} is the path to the directory of Manta software.

${bulk_bam} is the path to bulk DNA bam file.

${cell_bam} is the path to single cell bam file.

${ref} is the path to reference genome fasta file.

${bulk} is the bulk DNA sample name.

${cell} is the single cell sample name.


# USAGE STEP 2. SV-type annotation

2.1 for germline SVs, run the following shell scripts in Linux,

zcat ${sn}/results/variants/diploidSV.vcf.gz > ${sn}/results/variants/diploidSV.vcf

Rscript ./manta-annotate-concensus.R ${sn}

python3 ./manta_keepmate1.py -i ${sn}.manta.concensus.vcf -o ${sn}.manta.sc3input.vcf

In the above,

${sn} is the sample name.

2.2 for somatic SVs, run the following shell scripts in Linux, 

zcat ${bulk}.${cell}/results/variants/somaticSV.vcf.gz | grep "#" > ${bulk}.${cell}/results/variants/somaticSV.pass.vcf

zcat ${bulk}.${cell}/results/variants/somaticSV.vcf.gz | grep -v "#" | awk 'length($1)<=5 && $7=="PASS"' >> ${bulk}.${cell}/results/variants/somaticSV.pass.vcf

Rscript ./manta-annotate-concensus-somatic.R ${bulk}.${cell}

python3 ./manta_keepmate1.py -i ${bulk}.${cell}.manta.concensus.vcf -o ${bulk}.${cell}.manta.sc3input.vcf

In the above,

${bulk} is the bulk DNA sample name.

${cell} is the single cell sample name.


# USAGE STEP 3. enhanced reference genome construction, realignment and SV calling

3.1 for germline SVs, run the following shell scripts in Linux,

${vcf}=$(pwd)/${sn}.manta.sc3input.vcf

python3 pea_s3p1.py \
  --cellid ${sn} \
  --dellyvcf ${vcf} \
  --wkdir ${sn} \
  --cbam ${bam} \
  -g ${ref} \
  -s ${s} -e ${e}

python3 pea_s3p2.py \
  --cellid ${sn} \
  --dellyvcf ${vcf} \
  --wkdir ${sn} \
  --cbam ${bam} \
  -g ${ref} \
  --hsnp ${hsnp} \
  -s ${s} -e ${e}

In the above,

${sn} is the sample name.

${bam} is the path to sample bam file.

${ref} is the path to reference genome fasta file.

${hsnp} is the path to a list of heterozygous SNPs in the sample in vcf format.

${s} is the line number of the first SV in the ${bulk}.${cell}.manta.sc3input.vcf that a user want to analyze; ${e} is the line number of the last SV in the ${bulk}.${cell}.manta.sc3input.vcf that a user want to analyze.

3.2 for somatic SVs, run the following shell scripts in Linux, 

${vcf}=$(pwd)/${bulk}.${cell}.manta.sc3input.vcf

python3 pea_s3p1.py \
  --cellid ${cell} \
  --dellyvcf ${vcf} \
  --wkdir ${bulk}.${cell} \
  --cbam ${cell_bam} \
  -g ${ref} \
  -s ${s} -e ${e}

python3 pea_s3p2.py \
  --cellid ${cell} \
  --dellyvcf ${vcf} \
  --wkdir ${bulk}.${cell} \
  --cbam ${cell_bam} \
  -g ${ref} \
  --hsnp ${hsnp} \
  -s ${s} -e ${e}

In the above,

${bulk} is the bulk DNA sample name.

${cell} is the single cell sample name.

${cell_bam} is the path to single cell bam file.

${ref} is the path to reference genome fasta file.

${hsnp} is the path to a list of heterozygous SNPs in the sample in vcf format.

${s} is the line number of the first SV in the ${bulk}.${cell}.manta.sc3input.vcf that a user want to analyze; ${e} is the line number of the last SV in the ${bulk}.${cell}.manta.sc3input.vcf that a user want to analyze.


# Expected output.

The above output files of SVs 

Germline SVs: ${sn}/sv_temp_${s}_${e}/pool.vcf

Somatic SVs: ${bulk}.${cell}/sv_temp_${s}_${e}/pool.vcf

with the following information,

GT: genotype

thetaA: theta for the enhance reference genome C1-CSV pair

thetaB: theta for the enhance reference genome C2-CSV pair

refAAcount: no. reads supporting ref genotype in the enhance reference genome C1-CSV pair

refASVcount: no. reads supporting SV genotype in the enhance reference genome C1-CSV pair

refBBcount: no. reads supporting ref genotype in the enhance reference genome C2-CSV pair

refBSVcount: no. reads supporting SV genotype in the enhance reference genome C2-CSV pair

LA0:LA1:LB0:LB1: likelihoods for models h0, h1 based on the C1-CSV pair, and models h0, h1 based on the C2-CSV pair, respectfully.


# Release Notes

v3.1.0, 2020.08.01, revised significantly the PEA method to enable a full genome coverage on SV calling; and kept only PEA method in this release (without previous method for SNV and INDEL analysis)

v3.0.0, 2020.01.01, allowed genome structure variation calling using the PEA method.

v2.0.0, 2019.04.01, allowed parallel processing, optimized I/O, optimized pipeline, output in vcf format, and fixed bugs

v1.21, 2018.08.18, fixed bugs

v1.2, 2017.05.01, allowed INDEL calling, release version

v1.1.3, 2017.01.09, users can change the min mapQ, default to 40

v1.1.2, 2016.12.30, fixed bugs

v1.1.1, 2016.12.29, update read_mpileup to consider indels

V1.1, 2016.07.25, fixed bugs, release version

v1.0, 2016.04.26, release version

v0.0.4, 2016.04.26, fixed bugs in readling mpileup file

v0.0.3, 2016.04.22, read_mpilup function returns mindepth fails before returning reference genotype

v0.0.3, 2016.04.22, default mapQ change from 20 to 40

v0.0.2, 2016.04.19, fixed bugs - jump mpileup file column not fit problem.

v0.0.1, 2016.03, added likelihood ratio test based on null distribution from the data resampling.


