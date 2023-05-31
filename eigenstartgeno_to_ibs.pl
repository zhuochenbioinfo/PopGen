use strict;
use warnings;

my($prefix,$outpath) = @ARGV;

my $usage = "USAGE:\nperl $0 <geno file prefix> <output path>\n";
$usage .= "<prefix>: eigenstart geno file contains three files, named as [prefix].geno [prefix].snp and [prefix].ind\n";

die $usage unless(@ARGV == 2);
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

open(IN,"<$prefix.geno") or die $!;

my %hash_compare;

my $count = 0;

while(<IN>){
	chomp;
	my @genos = split//,$_;
	for(my $i = 0; $i < @samples; $i++){
		my $sample1 = $samples[$i];
		my $geno1 = $genos[$i];
		for(my $j = $i + 1; $j < @samples; $j++){
			my $sample2 = $samples[$j];
			my $geno2 = $genos[$j];
			my $compare = "NA";
			if($geno1 ne "9" and $geno2 ne "9"){
				if($geno1 == 1 or $geno2 == 1){
					if($geno1 == $geno2){
						$compare = "HET_iden";
					}else{
						$compare = "HET_diff";
					}
				}else{
					if($geno1 == $geno2){
						$compare = "HOMO_iden";
					}else{
						$compare = "HOMO_diff";
					}
				}
			}
			$hash_compare{$sample1}{$sample2}{compare}{$compare} ++;
			$hash_compare{$sample2}{$sample1}{compare}{$compare} ++;
		}
	}
	$count++;
	if(($count % 1000) == 0){
		print "$count\n";
	}
}
close IN;

my @types = ("HOMO_iden","HOMO_diff","HET_iden","HET_diff","NA");

foreach my $sample1(sort keys %hash_compare){
	open(OUT,">$outpath/$sample1.ibs");
	print OUT "sample1\tsample2\thomo_iden\thomo_diff\thet_iden\thet_diff\tNA\n";
	foreach my $sample2(sort keys %{$hash_compare{$sample1}}){
		my $line = "$sample1\t$sample2";
		foreach my $compare(@types){
			my $count = 0;
			if(exists $hash_compare{$sample1}{$sample2}{compare}{$compare}){
				$count = $hash_compare{$sample1}{$sample2}{compare}{$compare};
			}
			$line .= "\t$count";
		}
		print OUT "$line\n";
	}
	close OUT;
}

sub fileRowCount{
    my($file) = @_;
    my($num) = `wc -l $file|awk '{printf(\$1)}'`;
    chomp($num);
    return($num);
}
