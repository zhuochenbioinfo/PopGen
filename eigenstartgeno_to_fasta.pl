use strict;
use warnings;

my($prefix,$out,$misr) = @ARGV;

my $usage = "USAGE:\nperl $0 <geno file prefix> <output file> <missing rate>\n";
$usage .= "<prefix>: eigenstart geno file contains three files, named as [prefix].geno [prefix].snp and [prefix].ind\n";

die $usage unless(@ARGV >= 2);
die $usage unless(-e "$prefix.geno" and -e "$prefix.snp" and -e "$prefix.ind");

my $genoRow = fileRowCount("$prefix.geno");
my $snpRow = fileRowCount("$prefix.snp");
die "#WARNING: geno file and snp file shall have the same number of row.\n" unless($genoRow == $snpRow);

my @samples = ();
open(IN,"<$prefix.ind") or die $!;
while(<IN>){
	chomp;
	push @samples, $_;
}
close IN;

my @seqs = ();
open(IN,"paste $prefix.snp $prefix.geno|");
while(<IN>){
	chomp;
	my($chr,$pos,$ref,$alt,$genos) = split/\t/;
	next unless(length($ref) == 1 and length($alt) == 1);
	my @genos = split//,$genos;
	if($. == 1){
		unless(@genos == @samples){
			die "#WARNING: geno file and sample num conflict.\n";
		}
	}
	for(my $i = 0; $i < @genos; $i++){
		my $base = egeno2base($genos[$i],$ref,$alt);
		$seqs[$i] .= $base;
	}
}
close IN;

open(OUT,">$out");
for(my $i = 0; $i < @samples; $i++){
	my $sample = $samples[$i];
	my $seq = $seqs[$i];
	my $missRate = missingRatio($seq);
	if(defined $misr){
		if($missRate > $misr){
			print "#WARNING: skip sample:$sample for to many missing bases.\n";
			next;
		}
	}
	$seq = makefasta($seq);
	print OUT ">$sample\n$seq";
}
close OUT;

sub fileRowCount{
	my($file) = @_;
	my($num) = `wc -l $file|awk '{printf(\$1)}'`;
	chomp($num);
	return($num);
}

sub makefasta{
	my $seq = $_[0];
	$seq =~ s/(.{50})/$1\n/g;
	unless($seq =~ /\n$/){
		$seq .= "\n";
	}
	return($seq);
}

sub egeno2base{
	my($data,$ref,$alt) = @_;
	my $base = "-";
	if($data == 0){
		$base = $alt;
	}elsif($data == 2){
		$base = $ref;
	}
	return($base);
}

sub missingRatio{
	my($seq) = @_;
	my $count = ($seq =~ s/\-/#/g);
	my $len = length($seq);
	my $ratio = $count/$len;
	return($ratio);
}
