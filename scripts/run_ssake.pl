#!/usr/bin/perl -w
#
use strict;
use warnings;
use English;

#my $ssake_path="/data/eg03/sccaller/apps/SSAKE/ssake_v3-8-3/";
my $pre = -1;
my @reads = ();
my $len_e = -1;
my $file = shift or die "$0 <bp_reads.txt>\n";
my $ssake_path = shift;
my $counter = -1;
my $tmp_str = "";
open RD, $file or die $!;

open OUT, ">$file.localreads.fa" or die $!;
my $old_fh = select(OUT);
$| = 1;
open OUT3, '>>', "$file.contig2reads" or die $!;
$| = 1;
select($old_fh);
while (<RD>) {
	chomp;
	my @e = split /\s+/, $_;
	$len_e=@e;

	# meet a different group
	if ($pre != $e[0]) {
		if (@reads > 0) {
			#my $i = 0;
			seek OUT, 0, 0;
			truncate OUT, 0;
			my %seen = ();

			@reads = grep { ! $seen{(split "\n",$_)[1]}++ } @reads; # remove duplicate seq
			foreach my $r (@reads) {
				#print OUT ">read$i\n";
				#print OUT $r, "\n";
				#$i++;
				print OUT $r;
			}
			my $consensus = "";
			my $name = "";
			system("${ssake_path}SSAKE -m 16 -f $file.localreads.fa -w 1 -z 50 -o 1 -b $file.ssake.asm -c 1");
			# contigs -> OUT2
			open OUT2, '>>', "$file.ssake.asm.out" or die $!;
			my $old_fh2 = select(OUT2);
			$| = 1;
			select($old_fh2);
			open IN, "$file.ssake.asm.contigs" or die $!;
			while (<IN>) {
				if (/^>(\S+)/) {
					if ($name ne "") {
						print OUT2 ">", $name, "\n", $consensus, "\n";
					}
					$name = "E".$pre."|".$1;
				} else {
					chomp;
					$consensus = $_;
				}
			}

			print OUT2 ">", $name, "\n", $consensus, "\n";

			close IN;
			close OUT2;

			# contig---read -> OUT3
			open IN2, "$file.ssake.asm.readposition" or die $!;
			while (<IN2>) {
				if (/^>(\S+)/){
					print OUT3 ">E".$pre."|".$1."\n";
				}else{
					print OUT3 $_;
				}
			}
			close IN2;
			
	#		system("cat localreads.fa.cap.contigs >> cap3.asm.out");

			$pre = $e[0];
			@reads = (); # release old reads
		}
		#if($len_e==3){
		#	push @reads, ($e[1], $e[2]);
		#}
		#for($counter = 1; $counter < $len_e-1; $counter = $counter+1){
		#	push @reads, $e[$counter];
		#}
		$tmp_str = $e[3];
		substr($tmp_str, 0, 1, ">");
		push @reads, $tmp_str."+\n".$e[1]."\n";
		push @reads, $tmp_str."-\n".$e[2]."\n";
		$pre = $e[0];
	} else {
		#if($len_e==3){
		#	push @reads, ($e[1], $e[2]);
		#}
		#for($counter = 1; $counter < $len_e-1; $counter = $counter+1){
		#	push @reads, $e[$counter];
		#}
		$tmp_str = $e[3];
		substr($tmp_str, 0, 1, ">");
		push @reads, $tmp_str."+\n".$e[1]."\n";
		push @reads, $tmp_str."-\n".$e[2]."\n";
	}
}

if (@reads > 0) {
	#my $i = 0;
	seek OUT, 0, 0;
	truncate OUT, 0;
	my %seen = ();
	@reads = grep { ! $seen{(split "\n",$_)[1]}++ } @reads;
	foreach my $r (@reads) {
		#print OUT ">read$i\n";
		#print OUT $r, "\n";
		#$i++;
		print OUT $r;
	}
}
system("${ssake_path}SSAKE -m 16 -f $file.localreads.fa -w 1 -z 50 -o 1 -b $file.ssake.asm -c 1");
#contigs -> OUT2
open OUT2, '>>', "$file.ssake.asm.out" or die $!;
open IN, "$file.ssake.asm.contigs" or die $!;
my $consensus = "";
my $name = "";
while (<IN>) {
	if (/^>(\S+)/) {
		if ($name ne "") {
			print OUT2 ">", $name, "\n", $consensus, "\n";
		}
		$name = "E".$pre."|".$1;
	} else {
		chomp;
		$consensus = $_;
	}
}
print OUT2 ">", $name, "\n", $consensus, "\n";
close IN;
close OUT2;
# contig---read -> OUT3
open IN2, "$file.ssake.asm.readposition" or die $!;
while (<IN2>) {
	if (/^>(\S+)/){
		print OUT3 ">E".$pre."|".$1."\n";
	}else{
		print OUT3 $_;
	}
}
close IN2;
		#system("cat localreads.fa.cap.contigs >> cap3.asm.out");
close OUT;
close RD;
close OUT3;
