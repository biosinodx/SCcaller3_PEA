#!/bin/bash

if [ $# != 5 -a $# != 6 ]; then
	echo $0 \<main_path\> \<ref\> \<tumor_bam\> \<normal_bam\> \<n_CPUs:INT\> \[outputdir:-PWD\]
	exit 1
fi

nbbin=`readlink -f $1`
ref=`readlink -f $2`
tumor_bam=`readlink -f $3`
normal_bam=`readlink -f $4`
n_cpus=$5

if [ $# == 6 ]; then
	output=`readlink -f $6`
fi
lastdir=`pwd`
novobreak=$lastdir"/bin/novoBreak"
bwa=$nbbin/bin/bwa
samtools=$nbbin/bin/samtools
handle_fq=$lastdir/scripts/sc3.0.0.py
ssakepath=$nbbin/bin/ssake_v3-8-3/
echo "inside_dir="$lastdir
echo "ssakepath="$ssakepath
if [ $# == 6 ]; then
	mkdir -p $output
	cd $output
fi

# bam---> kmer so.bam ger.bam
$novobreak -i $tumor_bam -c $normal_bam -r $ref  -o kmer.stat 
echo after novobreak

mkdir -p group_reads
cd group_reads

if [ ! -f "somaticreads.srtnm.bam" ]; then 
	echo step 2: sorting bam file ...
    $samtools sort -n -@ $n_cpus -o somaticreads.srtnm.bam ../somaticreads.bam
else
	echo step 2: somaticreads.srtnm.bam already exist. Go to step 3.
fi

if [[ ! -f read1org.fq && ! -f read2org.fq ]]; then
	echo step 3: running samtools bam2fq ...
	$samtools bam2fq -1 read1org.fq -2 read2org.fq  somaticreads.srtnm.bam
else
	echo step 3: read1org.fq and read2org.fq already exist. Go to step 4.
fi

if [[ ! -f read1.fq && ! -f read2.fq ]]; then
    echo step 4: getting intersection of read1org.fq and read2org.fq ...
    python ${handle_fq} intersection read1org.fq read2org.fq read1.fq read2.fq
else
	echo step 4: read1.fq and read2.fq already exist. Go to step 5.
fi

if [[ ! -f bp_reads.txt ]]; then
    echo step 5: grouping reads ...
    perl ${lastdir}/scripts/group_bp_reads.pl ../kmer.stat read1.fq read2.fq > bp_reads.txt
else
	echo step 5: bp_reads.txt already exist.  Go to step 6.
fi

mkdir -p ../ssake
if [[ ! -f ../ssake/ssake.fa ]]; then
	echo step 6: splitting bp_reads.txt into n_cpu files ...
	cls=`tail -1 bp_reads.txt | cut -f1`
	rec=`echo $cls/8 | bc`
	rec=$((rec+1))
	echo "splitting bp_reads.txt has cls="${cls}" groups. Now splitting them into [8] files ..."
	mkdir -p split
	cd split
	rm -rf *
	awk -v rec=$rec '{print > int($1/rec)".txt"}' ../bp_reads.txt

	echo step 7: running ssake to assemble ...
	for file in *.txt
	do
	    perl $lastdir/scripts/run_ssake.pl $file $ssakepath> /dev/null &
		echo ${file} running ssake
	done
	wait
	echo after ssake
	cd ../../ssake/
	#you can split the bp_reads.txt into multiple files to run them together
	#perl $nbbin/run_ssake.pl ../group_reads/bp_reads.txt > /dev/null
	awk 'length($0)>1' ../group_reads/split/*.ssake.asm.out > ssake.fa
	cat ../group_reads/split/*.contig2reads > contig2reads
else
	cd ../ssake/
	echo step 6, 7: ssake.fa already exist. Skip step 6 and step 7. Go to step 8.
fi

if [[ ! -f ssake.sam ]]; then
    echo step 8: running bwa ...
    $bwa mem -t $n_cpus -M $ref ssake.fa > ssake.sam
else
	echo step 8: ssake.sam already exist. Go to step 9.
fi

if [[ ! -f ssake.vcf ]]; then
	echo step 9: running infer_sv.pl ...
    perl ${lastdir}/scripts/infer_sv.pl ssake.sam > ssake.vcf
else
	echo step 9: ssake.vcf already exist. Go to step 10.
fi

if [[ ! -f ssake.pass.vcf ]]; then
	echo step 10: grepping ...
    grep -v '^#' ssake.vcf | sed 's/|/\t/g' | sed 's/read//' |  awk '{if(!x[$1$2]){y[$1$2]=$14;x[$1$2]=$0}else{if($14>y[$1$2]){y[$1$2]=$14; x[$1$2]=$0}}}END{for(i in x){print x[i]}}' | sort -k1,1 -k2,2n  | perl -ne 'if(/TRA/){print}elsif(/SVLEN=(\d+)/){if($1>100){print $_}}elsif(/SVLEN=-(\d+)/){if($1>100){print}}' > ssake.pass.vcf
else
	echo step 10: ssake.pass.vcf already exist. Go to final step.
fi

if [[ ! -f ../novoBreak.pass.flt.vcf ]]; then
	#you can split the ssake.pass.vcf into multiple files to run them together
	num=`wc -l ssake.pass.vcf | cut -f1 -d' '`
	rec=`echo $num/$n_cpus | bc`
	rec=$((rec+1))
	mkdir -p split
	cd split
	split -l $rec ../ssake.pass.vcf # set proper split parameters when needed
	for file in x??
	do
		echo handling $file ...
		perl ${lastdir}/scripts/infer_bp_v4.pl $file $tumor_bam $normal_bam $nbbin > $file.sp.vcf &
	done
	wait
	echo after inter_bp_v4
	cd ..

	##below is a naive filter, pay attention to it
	#perl $nbbin/filter_sv_icgc.pl nbasm.pass.sp.vcf > ../novoBreak.pass.flt.vcf
	grep '^#' ssake.vcf > header.txt	
	perl ${lastdir}/scripts/filter_sv_icgc.pl split/*.sp.vcf | cat header.txt - > ../novoBreak.pass.flt.vcf
	echo after filter_sv_icgc
else
	echo novoBreak.pass.flt.vcf already exist. 
fi

cd $lastdir
echo all done
