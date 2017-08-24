use strict;
use warnings;
use Getopt::Long;

my($in,$out,$mafCut,$covCut);

my $usage = "USAGE:\nperl $0 --in <hapmap file> --out <outname> --maf <min maf> --cov <min cov>\n";
$usage .= "<min maf>: default=0\n";
$usage .= "<min cov>: default=0\n";

GetOptions(
	"in=s" => \$in,
	"out=s" => \$out,
	"maf=s" => \$mafCut,
	"cov=s" => \$covCut,
) or die $usage;

die $usage unless(defined $in and defined $out);

unless(defined $mafCut){
	$mafCut = 0;
}
unless(defined $covCut){
	$covCut = 0;
}

open(IN,"<$in") or die $!;
open(OUT,">$out.geno");
open(IND,">$out.ind");
open(SNP,">$out.snp");

while(<IN>){
	chomp;
	my($chr,$pos,$refBase,$altBase,$hitNum,$maf,$refNum,$snpNum,@samples) = split/\t/;
	if($. == 1){
		print IND join("\n",@samples)."\n";
		close IND;
		next;
	}
	next unless($hitNum/@samples >= $covCut and $maf >= $mafCut);
	print SNP "$chr\t$pos\t$refBase\t$altBase\n";
	my $outline = "";
	foreach my $geno(@samples){
		my $value = 9;
		if($geno eq $refBase){
			$value = 2;
		}elsif($geno eq $altBase){
			$value = 0;
		}
		$outline .= $value;
	}
	print OUT $outline."\n";
}
close OUT;
close SNP;
