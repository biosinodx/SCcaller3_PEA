#!/usr/bin/perl -w
#
#Auther: Zechen Chong
#
use strict;
use warnings;

my @dlines = ();
my @org_cigar_list = ();
my $id = 0;
print qq(##fileformat=VCFv4.1
##phasing=none
##INDIVIDUAL=TRUTH
##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="bamsurgeon spike-in">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">
##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description="Germline structural variant">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=UNKNOWN_LEN,Number=0,Type=Flag,Description="Unknown the length of SV">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Breakpoint consensus sequence">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SPIKEIN
);

open FP_OUT_SV2A,  ">sv2alignments" or die $!;
$| = 1;
my $sv2a = "";
my $prime_info = "";
my $pre = "";
my @results = ();
my @vcf_list = ();
my $org_cigar = "";
while (<>) {
	next if /^@/;
	my @e = split /\s+/, $_;
	next if $e[1] == 4 or $e[5]!~/[SH]/;

	# build new cigar
	my $newcigar = "";
	my $match = 0;
	while ($e[5] =~ /(\d+)([SMIDH])/g) {
		my ($l, $t) = ($1, $2);
		if ($t eq 'S' or $t eq 'H') {
			if (!$match) {
				$newcigar .= $l.$t;
			} else {
				$newcigar .= $match."M";
				$newcigar .= $l.$t;
				$match = 0;
			}
		}
		$match += $l if $t eq 'M';
		if ($t eq 'I' or $t eq 'D') {
			if ($l > 20) { # TODO: magic number
				$newcigar .= $match."M";
				$newcigar .= $l.$t;
				$match = 0;
			} else {
				$match += $l if $t eq 'I';
			}
		}
	}
	$newcigar .= $match."M" if $match;
	$org_cigar = $e[5];
	$e[5] = $newcigar;
	#if ($e[5] =~ /(\d+)[SH].+M(\d+)[SH]/) {
	#	next if ($1>5 and $2>5);
	#}

	# screen cigar
	if ($e[5] =~ /(\d+)M.+?(\d+)M/) {
		next if ($1>10 and $2>10);  # skip 10S11M20I11M
	}

	$_ = join("\t", @e);
	if ($e[0] ne $pre) {  # meet a new contig
		if ($pre ne "") {
			if (@dlines >= 2) {
				my $seq = &pick_consensus(\@dlines,\$prime_info, \@org_cigar_list);
				for (my $i = 0; $i < scalar @dlines - 1; $i++) {
					for (my $j = $i+1; $j < scalar @dlines; $j++) {
						@vcf_list = &parse_bp1($dlines[$i], $dlines[$j], $seq, \$sv2a, $org_cigar_list[$i], $org_cigar_list[$j]);
						push @results, @vcf_list if ($sv2a);
						print FP_OUT_SV2A $sv2a, "\t", $prime_info, "\n" if ($sv2a);
					}
				}
			}
			@dlines = ();
			@org_cigar_list = ();
		}
		$pre = $e[0];
		if ($e[5] =~ /(\d+)[SH].+M(\d+)[SH]/) {
			my $tmp = $e[5];
			if ($1 > 20) {
				$e[5] =~ s/$2[SH]//;  # delete
				$_ = join ("\t", @e);
				push @dlines, $_;
				push @org_cigar_list, $org_cigar;
				$e[5] = $tmp;
			}
			if ($2 > 20) {
				$e[5] =~ s/$1[SH]//;
				$_ = join ("\t", @e);
				push @dlines, $_;
				push @org_cigar_list, $org_cigar;
			}
			$e[5] = $tmp;
		} else {
			push @dlines, $_;
			push @org_cigar_list, $org_cigar;
		}
	} else {  # the same contig
		if ($e[5] =~ /(\d+)[SH].+M(\d+)[SH]/) {
			my $tmp = $e[5];
			if ($1 > 20) {
				$e[5] =~ s/$2[SH]//;
				$_ = join ("\t", @e);
				push @dlines, $_;
				push @org_cigar_list, $org_cigar;
				$e[5] = $tmp;
			}
			if ($2 > 20) {
				$e[5] =~ s/$1[SH]//;
				$_ = join ("\t", @e);
				push @dlines, $_;
				push @org_cigar_list, $org_cigar;
			}
			$e[5] = $tmp;
		} else {
			push @dlines, $_;
			push @org_cigar_list, $org_cigar;
		}
	}
}  # end of while

if (@dlines >= 2) {
	my $seq = &pick_consensus(\@dlines,\$prime_info,\@org_cigar_list);
	for (my $i = 0; $i < scalar @dlines - 1; $i++) {
		for (my $j = $i+1; $j < scalar @dlines; $j++) {
			@vcf_list = &parse_bp1($dlines[$i], $dlines[$j], $seq, \$sv2a, $org_cigar_list[$i], $org_cigar_list[$j]);
			push @results, @vcf_list if ($sv2a);
			print FP_OUT_SV2A $sv2a, "\t", $prime_info, "\n" if ($sv2a);
		}
	}
}
close FP_OUT_SV2A;

@results = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @results;
my $prepos = 0;
for (my $i = 0; $i < @results; $i++) {
	print join("\t", @{$results[$i]}), "\n";
}

1;

sub pick_consensus { # same contig
	my $lines = shift;
	my $prime_info = shift;
	my $org_cigar_list = shift;
	my $orgcigar = "";
	if (scalar @$lines != scalar @$org_cigar_list){
		print "not equal len(lines)=".scalar @$lines." and len(org_cigar_list)=".scalar @$org_cigar_list;
	}
	for (my $index=0;$index<scalar @$lines;$index++){
		my @e = split /\s+/, $$lines[$index];
		if ($e[1] & 256){  # not primary alignment
			next;
		}

		$orgcigar = $$org_cigar_list[$index];
		# chr, pos, flag, seq, cigar, org_cigar
		$$prime_info = join("\t", ($e[2],$e[3],$e[1],$e[9],$e[5], $orgcigar));
		return $e[9]; # there should be only 1 primary alignment
	}
	# foreach my $l (@$lines) {
	# 	my @e = split /\s+/, $l;
	# 	if ($e[1] & 256){  # not primary alignment
	# 		$counter++;
	# 		next;
	# 	}
	# 	$orgcigar = $$org_cigar_list[$counter];
	# 	# chr, pos, flag, seq, cigar, org_cigar
	# 	$$prime_info = join("\t", ($e[2],$e[3],$e[1],$e[9],$e[5], $orgcigar));
	# 	return $e[9]; # there should be only 1 primary alignment
	# }
}

sub parse_ins {
	my $l = shift;
	my @lines = @$l;
	my $seq = &pick_consensus(@lines);
	my @bps = ();
	my @ret = ();
	for (my $i = 0; $i < @lines; $i++) {
		my @e = split /\s+/, $lines[$i];
		my ($m1, $s1) = (0, 0);
		while ($e[5] =~ /(\d+)[SH]/g) {
				$s1 = $1 if $1 > $s1;
		}
		while ($e[5] =~ /(\d+)M/g) {
				$m1 = $1 if $1 > $m1;
		}
		if ($e[5] =~ /$m1[M].*?$s1[SH]/) {
			$e[3] = $e[3]+$m1-1;
		}
		push @bps, [ ($e[2],$e[3],$e[1], $e[4],$e[5], $e[0]) ];	
	#	print join("\t", $e[2],$e[3], $e[4], $e[5]), "\n";
	}
	@bps = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @bps;
	#foreach (@bps) {
	#	my @rec = @$_;
	#	print join("\t", @rec), "\n";
	#}
	my @pre = ();
	my @cur = ();
	foreach (@bps) {
		@cur = @$_;
		if (@pre) {
			if ($cur[0] eq $pre[0] and $cur[1]-$pre[1]<3) { #TODO magic number
				#print join("\t", $pre[0],$pre[1], $pre[2], $pre[3], $pre[4]), "\n";
				#print join("\t", $cur[0],$cur[1], $cur[2], $cur[3], $cur[4]), "\n";
				my ($m1, $s1, $m2, $s2) = (0, 0, 0, 0);
				while ($pre[4] =~ /(\d+)[SH]/g) {
					$s1 = $1 if $1 > $s1;
				}
				while ($pre[4] =~ /(\d+)M/g) {
					$m1 = $1 if $1 > $m1;
				}
				while ($cur[4] =~ /(\d+)[SH]/g) {
					$s2 = $1 if $1 > $s2;
				}
				while ($cur[4] =~ /(\d+)M/g) {
					$m2 = $1 if $1 > $m2;
				}
				#print "MS", $m1, "\t", $s1, "\t", $m2, "\t", $s2, "\n";
				if (($pre[4]=~/$m1[M].*?$s1[SH]/ and $cur[4]=~/$s2[SH].*?$m2[M]/) or ($pre[4]=~/$s1[SH].*?$m1[M]/ and $cur[4]=~/$m2[M].*?$s2[SH]/)) {
					if (($pre[2]&16 and $cur[2]&16) or (!($pre[2]&16) and !($cur[2]&16))) { # same strand
						#print join("\t", $cur[0], $pre[1], ".", ".", "<INS>", ($pre[3]+$cur[3])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=INS;END=$cur[1];SVLEN=".($cur[1]-$pre[1]), "GT", "./."), "\n";
						push @ret, [ ($cur[0], $pre[1], "N", ".", "<INS>", ($pre[3]+$cur[3])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=INS;CHR2=$cur[0];END=$cur[1];UNKNOWN_LEN;SVLEN=".($cur[1]-$pre[1]), "GT", "./.", $cur[-1]) ];
					}
				}
				if (($pre[4]=~/$m1[M].*?$s1[SH]/ and $cur[4]=~/$m2[M].*?$s2[SH]/) or ($pre[4]=~/$s1[SH].*?$m1[M]/ and $cur[4]=~/$s2[SH].*?$m2[M]/)) {
					if (($pre[2]&16 and !($cur[2]&16)) or (!($pre[2]&16) and $cur[2]&16)) { # different strand
						#print join("\t", $cur[0], $pre[1], ".", ".", "<INS>", ($pre[3]+$cur[3])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=INS;END=$cur[1];SVLEN=".($cur[1]-$pre[1]), "GT", "./."), "\n";
						push @ret, [ ($cur[0], $pre[1], "N", ".", "<INS>", ($pre[3]+$cur[3])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=INS;CHR2=$cur[0];END=$cur[1];UNKNOWN_LEN;SVLEN=".($cur[1]-$pre[1]), "GT", "./.", $cur[-1]) ];
					}
				}
			} 
			@pre = @cur;
		} else {
			@pre = @cur;
		}
		#print join("\t", $_->[0],$_->[1], $_->[2], $_->[3], $_->[4]), "\n";
	}
	#print @bps;
	return @ret;
}

sub parse_bp1 { # same contig
	my $line1 = shift;
	my $line2 = shift;
	my $seq = shift;
	my $sv2a_line = shift;
	my $orgcigar1 = shift;
	my $orgcigar2 = shift;
	my @ret = ();
	my ($m1, $s1, $m2, $s2) = (0, 0, 0, 0);
	my @e1 = split /\s+/, $line1;
	my @e2 = split /\s+/, $line2;
	$$sv2a_line = "";
	
	if ($e1[5] =~ /[SMH]/ and $e2[5] =~ /[SMH]/) {
			while ($e1[5] =~ /(\d+)[SH]/g) {
					$s1 = $1 if $1 > $s1;
			}
			while ($e1[5] =~ /(\d+)M/g) {
					$m1 = $1 if $1 > $m1;
			}

			while ($e2[5] =~ /(\d+)[SH]/g) {
					$s2 = $1 if $1 > $s2;
			}
			while ($e2[5] =~ /(\d+)M/g) {
					$m2 = $1 if $1 > $m2;
			}
	}
	my ($pos1, $pos2) = (0, 0);
	if ($e1[5] =~ /$s1[SH].*?$m1[M]/) {
		$pos1 = $e1[3];  # 50S241M
	} else {
		$pos1 = $e1[3]+$m1-1;  # 241M50S
	}
	if ($e2[5] =~ /$s2[SH].*?$m2[M]/) {
		$pos2 = $e2[3];
	} else {
		$pos2 = $e2[3]+$m2-1;
	}

	my $chr_1st = "";
	my $pos_1st = "";
	my $chr_2nd = "";
	my $pos_2nd = "";
	my $flag_1st = "";
	my $flag_2nd = "";
	my $seq_1st = "";
	my $seq_2nd = "";
	my $cigar_1st = "";
	my $cigar_2nd = "";
	my $left_num = 1;
	my $org_cigar_1st = "";
	my $org_cigar_2nd = "";
	my $vcf_pos1 = "";
	my $vcf_pos2 = "";
	if (abs($m1-$s2)<=25 or abs($m2-$s1)<=25) {
		if ($e1[2] ne $e2[2]) { # different chr: trans

			if ((($e1[1] & 0x10) ^ ($e2[1] & 0x10)) == 0) { # same strand +,+ or -,-
				if ($e1[1] & 0x10) { #-,-
					if ($pos1 == $e1[3]) {
						push @ret, [ ($e1[2], $pos1, "N", ".", "<TRA>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=5to3;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e2[2];END=$pos2;SVLEN=0", "GT", "./.", $e1[0]) ];
						$left_num = 1;
					}
					if ($pos2 == $e2[3]) {
						push @ret, [ ($e2[2], $pos2, "N", ".", "<TRA>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=5to3;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e1[2];END=$pos1;SVLEN=0", "GT", "./.", $e2[0]) ];
						$left_num = 2;
					}
				} else { #+,+
					if ($pos2 == $e2[3]) {
						push @ret, [ ($e1[2], $pos1, "N", ".", "<TRA>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=3to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e2[2];END=$pos2;SVLEN=0", "GT", "./.", $e1[0]) ];
						$left_num = 1;
					}
					if ($pos1 == $e1[3]) {
						push @ret, [ ($e2[2], $pos2, "N", ".", "<TRA>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=3to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e1[2];END=$pos1;SVLEN=0", "GT", "./.", $e2[0]) ];
						$left_num = 2;
					}
				}
			} else { # different strand +,- or -,+
				if ($e1[1] & 0x10) { # 1-, 2+
					if ($pos1 == $e1[3]) {
						push @ret, [ ($e1[2], $pos1, "N", ".", "<TRA>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=5to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e2[2];END=$pos2;SVLEN=0", "GT", "./.", $e1[0]) ];
						$left_num = 1;
					} else {
						push @ret, [ ($e2[2], $pos2, "N", ".", "<TRA>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=3to3;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e1[2];END=$pos1;SVLEN=0", "GT", "./.", $e2[0]) ];
						$left_num = 2;
					}
				} else { # 1+, 2-
					if ($pos1 == $e1[3]) {
						push @ret, [ ($e2[2], $pos2, "N", ".", "<TRA>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=5to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e1[2];END=$pos1;SVLEN=0", "GT", "./.", $e2[0]) ];
						$left_num = 2;
					} else {
						push @ret, [ ($e1[2], $pos1, "N", ".", "<TRA>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=3to3;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e2[2];END=$pos2;SVLEN=0", "GT", "./.", $e1[0]) ];
						$left_num = 1;
					}
				}
			}
		} else { # same chr

			if ((($e1[1] & 16) and !($e2[1] & 16)) or (!($e1[1] & 16)) and ($e2[1] & 16)) { # +,- or -,+ inv
				if(abs($pos2 - $pos1) < 10) { # TODO magic number controling inv size
				} 
				else {
					if ($pos1 < $pos2) { # TODO wrong CT look clip
						#print join("\t", $e2[2], $pos1, ".", ".", "<INV>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=INV;END=$pos2;SVLEN=0", "GT", "./."), "\n";
						push @ret, [ ($e2[2], $pos1, "N", ".", "<INV>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=5to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=INV;CHR2=$e2[2];END=$pos2;SVLEN=".abs($pos2-$pos1+1), "GT", "./.",$e2[0]) ] ;
						$left_num = 1;
					} else {
						#print join("\t", $e2[2], $pos2, ".", ".", "<INV>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=INV;END=$pos1;SVLEN=0", "GT", "./."), "\n";
						push @ret, [ ($e2[2], $pos2, "N", ".", "<INV>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=3to3;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=INV;CHR2=$e2[2];END=$pos1;SVLEN=".abs($pos2-$pos1+1), "GT", "./.",$e2[0]) ];
						$left_num = 2;
					}
				}
			} else { # +,+ or -,-
				#print $line1, $line2;
				if ($pos1 < $pos2) {
					if ((!($e1[1]&16) and $e1[5]=~/$m1[M].*?$s1[SH]/) or ($e1[1]&16 and $e1[5]=~/$s1[SH].*?$m1[M]/)) { # del
						#print join("\t", $e2[2], $pos1, ".", ".", "<DEL>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=DEL;END=$pos2;SVLEN=".($pos1-$pos2), "GT", "./."), "\n";
						if(abs($pos2 - $pos1) > 10) { # TODO magic number controling del size
							push @ret, [ ($e2[2], $pos1, "N", ".", "<DEL>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=3to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=DEL;CHR2=$e2[2];END=$pos2;SVLEN=".($pos1-$pos2), "GT", "./.",$e2[0])] ;
							$left_num = 1;
						}
					} else { # dup
						#print join("\t", $e2[2], $pos1, ".", ".", "<DUP>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=DUP;END=$pos2;SVLEN=".($pos2-$pos1), "GT", "./."), "\n";
						push @ret, [ ($e2[2], $pos1, "N", ".", "<DUP>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=5to3;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=DUP;CHR2=$e2[2];END=$pos2;SVLEN=".($pos2-$pos1), "GT", "./.",$e2[0]) ];
						$left_num = 1;
					}
				} else {
					if ((!($e2[1]&16) and $e2[5]=~/$m2[M].*?$s2[SH]/) or ($e2[1]&16 and $e2[5]=~/$s2[SH].*?$m2[M]/)) { # del
						#print join("\t", $e2[2], $pos2, ".", ".", "<DEL>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=DEL;END=$pos1;SVLEN=".($pos2-$pos1), "GT", "./."), "\n";
						if(abs($pos2 - $pos1) > 10) { # TODO magic number controling del size
							push @ret, [ ($e2[2], $pos2, "N", ".", "<DEL>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=3to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=DEL;CHR2=$e2[2];END=$pos1;SVLEN=".($pos2-$pos1), "GT", "./.",$e2[0]) ];
							$left_num = 2;
						}
					} else { # dup
						#print join("\t", $e2[2], $pos2, ".", ".", "<DUP>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=DUP;END=$pos1;SVLEN=".($pos1-$pos2), "GT", "./."), "\n";
						push @ret, [ ($e2[2], $pos2, "N", ".", "<DUP>", ($e1[4]+$e2[4])/2, "PASS", "PRECISE;CT=5to3;CIPOS=-10,10;CIEND=-10,10;SOMATIC;CONSENSUS=$seq;SVTYPE=DUP;CHR2=$e2[2];END=$pos1;SVLEN=".($pos1-$pos2), "GT", "./.",$e2[0]) ];
						$left_num = 2;
					}
				}
			}
		}
	}

	
	if (@ret){
		$chr_1st = ($left_num == 1)?$e1[2]:$e2[2];
		$chr_2nd = ($left_num == 1)?$e2[2]:$e1[2];
		$pos_1st = ($left_num == 1)?$e1[3]:$e2[3];
		$pos_2nd = ($left_num == 1)?$e2[3]:$e1[3];
		$flag_1st = ($left_num == 1)?$e1[1]:$e2[1];
		$flag_2nd = ($left_num == 1)?$e2[1]:$e1[1];
		$seq_1st = ($left_num == 1)?$e1[9]:$e2[9];
		$seq_2nd = ($left_num == 1)?$e2[9]:$e1[9];
		$cigar_1st = ($left_num == 1)?$e1[5]:$e2[5];
		$cigar_2nd = ($left_num == 1)?$e2[5]:$e1[5];
		$org_cigar_1st = ($left_num == 1)?$orgcigar1:$orgcigar2;
		$org_cigar_2nd = ($left_num == 1)?$orgcigar2:$orgcigar1;
		$vcf_pos1 = ($left_num == 1)?$pos1:$pos2;
		$vcf_pos2 = ($left_num == 1)?$pos2:$pos1;
		# congig name, chr1, pos1, flag1, seq1, cigar1, orgcigar1, vcf_pos1, chr2, pos2, flag2, seq2, cigar2, orgcigar2, vcf_pos2
		$$sv2a_line = join("\t",($e1[0], $chr_1st, $pos_1st, $flag_1st, $seq_1st, $cigar_1st, $org_cigar_1st, $vcf_pos1, $chr_2nd, $pos_2nd, $flag_2nd, $seq_2nd, $cigar_2nd, $org_cigar_2nd, $vcf_pos2));
	}
	return @ret;
}
