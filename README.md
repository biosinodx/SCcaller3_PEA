SCcaller with the PEA method
=
Single Cell Caller (SCcaller) - Identify single nucleotide variations (SNVs), short insertions, short deletions (INDELs) and genome structural variations (GSVs) from single cell sequencing data

Version 3.0.0

Updated date: 2020.01.01

Cite us:

For SNV and INDEL calling (the SCcaller method): Dong X et al. Accurate identification of single-nucleotide variants in whole-genome-amplified single cells. Nat Methods. 2017 May;14(5):491-493. doi: 10.1038/nmeth.4227.

For GSV calling (the PEA method): Dong X et al. Identifying genome structural variations in single cells. In submission.

#####
Author and License
-
Authors: Xiao Dong, Yujue Wang

Email: biosinodx@gmail.com (X.D.), xiao.dong@einsteinmed.org (X.D.), spsc83@gmail.com (Y.W.), yujue.wang@einsteinmed.org(Y.W.)

Licensed under the GNU Affero General Public License version 3 or later

#####
Dependence
-
    sv: python modules sys, logging, os, re, subprocess, time, copy, struct, itertools, operator, multiprocessing, socket, Queue, collections, numpy, math, datetime
        samtools v.1.9+ (Other versions not tested)
        bwa v.0.7.10-r806-dirty
        SSAKE v3.8.3
    si: python modules os, argparse, sys, subprocess, re, collections, itertools, logging, time, functiontools, random, string, math, numpy, multiprocessing, pysam(0.15.1)
        samtools v.1.9+ (Other versions not tested)

#####
Usage
-
    python sccaller.py <command> <parameters>
* Commands:
  * Structure variation calling (the PEA method)<br>
  
        sv  <file_ref>: reference file
            <cpu_num>: number of processes
            <job_size>: the candidate SVs will be divided in to pieces of job_size lines
            <single_cell_bam>: bam file of a single cell
            <single_cell_fq1>: the sequence file
            <single_cell_fq2>: the sequence file
            <bulk_bam>: bam file of bulk
            <hSNP>: hSNP file
            <r1>: the min number of reads support related SNP
            <r2>: the min number of reads do not support related SNP
            <output>: output file name
            <ifq_dir>: the directory contain the indexs of sequence files
  * SNVs and INDELs calling (the SCcaller method)<br>
    
          si  -b, --bam:      bam file of a single cell
              -f, --fasta:    fasta file of reference genome
              -o, --output:   output file name
              -s, --snp_in:   Candidate snp input file, either from dbsnp data or heterozygous snp (hsnp) data of the bulk, for known heterogous call. file type: bed (1-based) or vcf.
              -t, --snp_type: SNP type for --snp_in. It could be either "dbsnp" or "hsnp". When choosing dbsnp, --bulk bulk_bamfile is required.

            Optional arguments:
              --RD:           min. read depth of known heterogous SNP called from bulk when choosing -t dbsnp. Default: 20. Recommand: 10,15,20, depending on average read depth
              --bias:         default theta (bias) for SNVs whose theta cannot be estimated. Default=0.75
              --bulk:         bamfile of bulk DNA sequencing
              --bulk_min_depth:min. reads for bulk. Default: 20
              --bulk_min_mapq:min. mapQ for bulk. Default: 20
              --bulk_min_var: min. num. variant supporting reads for bulk. Default: 1
              --format:       output file format. bed or vcf. Default: vcf
              --head:         first chromosome as sorted as in fasta file to analyze (1-based). Default: the first chr. in the fasta
              --mapq:         min. mapQ. Default: 40
              --min_depth:    min. reads. Default: 10
              --minvar:       min. num. variant supporting reads. Default: 4
              --null:         min. allelic fraction considered. Default=0.03
              --tail:         last chromosome as sorted as in fasta file to analyze (1-based). Default: the last chr. in the fasta
              -d, --wkdir:    work dir. Default: ./
              -e, --engine:   pileup engine. samtools or pysam. Default: pysam
              -h, --help:     Help
              -l, --lamb:     lambda for bias estimation. Default=10000
              -n, --cpu_num:  num. processes. Default: 1
              -w, --work_num: num. splits per chromosome for multi-process computing. Default: 100


    * I. Basic usage: calling SNVs and INDELs from a cell

      * I.a When you have heterozygous SNPs pre-called from bulk DNA of the same subject,<br>
      
            python sccaller.py si\
              --bam cell.bam \ # bam file of a single cell
              --fasta ref.fa \ # reference genome in fasta format
              --output cell.vcf \ # output vcf file
              --snp_type hsnp \ # using heterozygous SNPs pre-called from bulk DNA
              --snp_in hsnp.vcf (or bed) \ # vcf or bed file of heterozygous SNPs pre-called from bulk DNA
              --cpu_num 8 \ # using 8 cpu threads
              --engine samtools # using samtools engine

      * I.b When you do not have heterozygous SNPs pre-called from bulk DNA of the same subject, obtain SNPs from dbSNP or other databases,<br>
      

            python sccaller.py si\
              --bam cell.bam \ # bam file of a single cell
              --fasta ref.fa \ # reference genome in fasta format
              --output cell.vcf \ # output vcf file
              --snp_type dbsnp \ # using SNPs from dbSNP database (or other database)
              --snp_in dbsnp.vcf (or bed) \ # vcf or bed file containing all SNPs in dbSNP (or other) database
              --cpu_num 8 \ # using 8 cpu threads
              --engine samtools # using samtools engine

    * II. Calling `somatic` SNVs and INDELs not present in bulk DNA

      * II.a Step 1. Calling SNVs and INDELs from a cell together with bulk DNA in input,<br>
      
            python sccaller.py si\
              --bam cell.bam \ # bam file of a single cell
              --bulk bulk.bam \ # bam file of bulk DNA
              --fasta ref.fa \ # reference genome in fasta format
              --output cell.vcf \ # output vcf file
              --snp_type hsnp \ # using heterozygous SNPs pre-called from bulk DNA
              --snp_in hsnp.vcf (or bed) \ # vcf or bed file of heterozygous SNPs pre-called from bulk DNA
              --cpu_num 8 \ # using 8 cpu threads
              --engine samtools # using samtools engine

      * II.b Step 2. Filtering out SNVs and INDELs observed in bulk or sequencing depth <= 20x in the single cell<br>
      
            grep '0/1' cell.vcf | grep 'True' | awk '$7=="." && length($5)==1' | awk -F "[:,]" '$8+$9>=20' > cell.somatic.snv.vcf
            grep '0/1' cell.vcf | grep 'True' | awk '$7=="." && length($5)>1' | awk -F "[:,]" '$8+$9>=20' > cell.somatic.indel.vcf

    * III. Notes on X/Y chromosome in males and ploidy<br>
      Please note, sccaller was designed assuming diploidy genome and two copies of chromosomes. It cannot be used for calling mutations from X/Y chromosome of a male subject.


#####
##RELEASE NOTES

v3.0.0, 2020.01.01, allow genome structure variation calling using the PEA method.

v2.0.0, 2019.04.01, allowing parallel processing, optimizing I/O, optimizing pipeline, output in vcf format, and fixing bugs

v1.21, 2018.08.18, fixing bugs

v1.2, 2017.05.01, allow INDEL calling, release version

v1.1.3, 2017.01.09, users can change the min mapQ, default to 40

v1.1.2, 2016.12.30, fixing bugs

v1.1.1, 2016.12.29, update read_mpileup to consider indels

V1.1, 2016.07.25, fixing bugs, release version

v1.0, 2016.04.26, release version

v0.0.4, 2016.04.26, fixed bugs in readling mpileup file

v0.0.3, 2016.04.22, read_mpilup function returns mindepth fails before returning reference genotype

v0.0.3, 2016.04.22, default mapQ change from 20 to 40

v0.0.2, 2016.04.19, fix bugs - jump mpileup file column not fit problem.

v0.0.1, 2016.03, add likelihood ratio test based on null distribution from the data resampling.    
