#!/usr/bin/perl -w
#
#Author: Ruan Jue
#
use strict;

# deletion, insertion, invertion, transloaction
my @probs = (0.10, 0.80, 0.10, 0.10);

my $ref  = shift or die("Usage: $0 <ref> <n_sv>\n");
my $n_sv = shift or die("Usage: $0 <ref> <n_sv>\n");

my %refseqs = ();

open(IN, $ref) or die;
my $name = '';
my $seq  = '';
while(<IN>){
	if(/^>(\S+)/){
		$refseqs{$name} = $seq if($seq);
		$name = $1;
		$seq = '';
	} else {
		chomp;
		$seq .= $_;
	}
}
$refseqs{$name} = $seq if($seq);
close IN;
my @genome = ();
my @names = ();
foreach my $chr (sort keys %refseqs){
	$seq = $refseqs{$chr};
	push(@names, $chr);
	push(@genome, [length($seq), [length($seq), $chr, "+", 0, length($seq)]]);
}

$| = 1;

for(my $i=0;$i<$n_sv;$i++){
	my $p = rand();
	my $type = 0;
	while($type < @probs and $p > $probs[$type]){ $p -= $probs[$type]; $type ++; }
	my $chr1 = int(rand(scalar(@genome)));
	my $chr2 = int(rand(scalar(@genome)));
	if($type == 0){ # deletion
		my ($off, $len);
		while(1){
			$off = int(rand($genome[$chr1][0]));
			$len = int(rand($genome[$chr1][0] / 1000));
			last if($off + $len < $genome[$chr1][0]);
		}
		&delete_seq($genome[$chr1], $off, $len);
		#print STDERR "DELETE [$chr1, $off, $len]\n";
	} elsif($type == 1){ #insertion
		my ($off1, $off2, $len);
		while(1){
			$off1 = int(rand($genome[$chr1][0]));
			$len  = int(rand($genome[$chr1][0] / 100));
			last if($off1 + $len < $genome[$chr1][0]);
		}
		my $seg = &delete_seq($genome[$chr1], $off1, $len);
		$off2 = int(rand($genome[$chr2][0]));
		&insert_seq($genome[$chr2], $off2, $seg);
		#print STDERR "INSERT [$chr1, $off1, $len] -> [$chr2, $off2]\n";
	} elsif($type == 2){ #inverse
		my ($off, $len);
		while(1){
			$off = int(rand($genome[$chr1][0]));
			$len = int(rand($genome[$chr1][0] / 1000));
			last if($off + $len < $genome[$chr1][0]);
		}
		my $seg = &delete_seq($genome[$chr1], $off, $len);
		&insert_seq($genome[$chr1], $off, &inverse_seq($seg));
		#print STDERR "INVERSE [$chr1, $off, $len]\n";
	} elsif($type == 3){ #translocation
		if($chr1 == $chr2){ $n_sv ++; }
		my ($off1, $dir1, $off2, $dir2);
		while(1){
			$off1 = int(rand($genome[$chr1][0]));
			$off2 = int(rand($genome[$chr2][0]));
			$dir1 = ($off1 > $genome[$chr1][0] / 2);
			$dir2 = ($off2 > $genome[$chr2][0] / 2);
			last;
		}
		my $seg1 = $dir1? &delete_seq($genome[$chr1], $off1, $genome[$chr1][0] - $off1) : &delete_seq($genome[$chr1], 0, $off1);
		my $seg2 = $dir2? &delete_seq($genome[$chr2], $off2, $genome[$chr2][0] - $off2) : &delete_seq($genome[$chr2], 0, $off2);
		$dir1? &insert_seq($genome[$chr1], $genome[$chr1][0], $seg2) : &insert_seq($genome[$chr1], 0, $seg2);
		$dir2? &insert_seq($genome[$chr2], $genome[$chr2][0], $seg1) : &insert_seq($genome[$chr2], 0, $seg1);
		#print STDERR "TRANSLOCATION [$chr1, $off1, $dir1] [$chr2, $off2, $dir2]\n";
	}
}

for(my $i=0;$i<@genome;$i++){
	my $chr = $genome[$i];
	print ">$names[$i]\tlength=$chr->[0]";
	$seq = '';
	for(my $j=1;$j<@$chr;$j++){
		print "\t" . join("|", @{$chr->[$j]});
		my $s = substr($refseqs{$chr->[$j]->[1]}, $chr->[$j]->[3], $chr->[$j]->[0]);
		if($chr->[$j]->[2] eq '-'){
			$s = reverse $s;
			$s=~tr/ACGTacgt/TGCAtgca/;
		}
		$seq .= $s;
	}
	print "\n";
	while($seq=~/(.{1,100})/g){ print "$1\n"; }
}

1;

sub print_chr {
	my $chr = shift;
	print "$chr->[0]";
	for(my $i=1;$i<@$chr;$i++){
		print "\t" . join("|", @{$chr->[$i]});
	}
	print "\n";
}

sub print_genome {
	my $g = shift;
	print "{\n";
	foreach my $chr (@$g){ &print_chr($chr); }
	print "}\n";
}

sub inverse_seq {
	my $chr = shift;
	my @ret = reverse @$chr;
	my $len = pop(@ret);
	foreach my $e (@ret){
		$e->[2] = ($e->[2] eq '+')? '-' : '+';
	}
	unshift(@ret, $len);
	return \@ret;
}

sub split_seq {
	my ($chr, $off) = @_;
	my $n_seg = @{$chr} - 1;
	my $i = 1;
	my $p = 0;
	while($i <= $n_seg and $p + $chr->[$i][0] <= $off){ $p += $chr->[$i][0]; $i ++; }
	if($off > $p){
		my @a = @{$chr->[$i]};
		my @b = ();
		my @c = ();
		$b[0] = $off - $p; $b[1] = $a[1]; $b[2] = $a[2];
		$c[0] = $a[0] - $b[0]; $c[1] = $a[1]; $c[2] = $a[2];
		if($a[2] eq '+'){
			$b[3] = $a[3]; $b[4] = $a[3] + $b[0];
			$c[3] = $b[4]; $c[4] = $a[4];
		} else {
			$b[4] = $a[4]; $b[3] = $b[4] - $b[0];
			$c[4] = $b[3]; $c[3] = $a[3];
		}
		splice(@{$chr}, $i, 1, \@b, \@c);
		$i ++;
	}
	return $i;
}

sub delete_seq {
	my ($chr, $off, $len) = @_;
	my $bk1 = &split_seq($chr, $off);
	my $bk2 = &split_seq($chr, $off + $len);
	my @del = splice(@{$chr}, $bk1, $bk2 - $bk1);
	my $ret = [$len, @del];
	$chr->[0] -= $len;
	return $ret;
}

sub insert_seq {
	my ($chr, $off, $new) = @_;
	my $bk = &split_seq($chr, $off);
	$chr->[0] += $new->[0];
	my $n_seg = @{$new} - 1;
	splice(@{$chr}, $bk, 0, @{$new}[1..$n_seg]);
}

