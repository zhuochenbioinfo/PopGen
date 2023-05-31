use strict;
use warnings;

my $inpath = shift;
my $outfile = shift;

my $suffix = "ibs";

open(IN,"cat $inpath/*/*.$suffix|") or die $!;
open(OUT,">$outfile");

my %hash_pair;

while(<IN>){
	chomp;
	if($. == 1){
		print OUT $_."\n";
	}
	next if($_ =~ /^sample/);
	my($sample1,$sample2,$homo_iden,$homo_diff,$het_iden,$het_diff,$na) = split/\t/;
	#next if(exists $hash_pair{$sample2}{$sample1});
	$hash_pair{$sample1}{$sample2}{homo_iden} += $homo_iden;
	$hash_pair{$sample1}{$sample2}{homo_diff} += $homo_diff;
	$hash_pair{$sample1}{$sample2}{het_iden} += $het_iden;
	$hash_pair{$sample1}{$sample2}{het_diff} += $het_diff;
	$hash_pair{$sample1}{$sample2}{na} += $na;
}
close IN;

foreach my $sample1(sort keys %hash_pair){
	foreach my $sample2(sort keys %{$hash_pair{$sample1}}){
		my $homo_iden = $hash_pair{$sample1}{$sample2}{homo_iden};
		my $homo_diff = $hash_pair{$sample1}{$sample2}{homo_diff};
		my $het_iden = $hash_pair{$sample1}{$sample2}{het_iden};
		my $het_diff = $hash_pair{$sample1}{$sample2}{het_diff};
		my $na = $hash_pair{$sample1}{$sample2}{na};
		print OUT "$sample1\t$sample2\t$homo_iden\t$homo_diff\t$het_iden\t$het_diff\t$na\n";
	}
}
close OUT;
