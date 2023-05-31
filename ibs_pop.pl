use strict;
use warnings;

my $ibs = shift;
my $list = shift;

open(IN,"<$list") or die $!;
my %hash_sample;
while(<IN>){
	chomp;
	my($sample,$group) = split/\t/;
	$hash_sample{$sample}{group} = $group;
}
close IN;


open(IN,"<$ibs") or die $!;

my %hash_pair;

while(<IN>){
	chomp;
	if($. == 1){
		print $_."\tpair\n";
		next;
	}
	my($sample1,$sample2,$v1,$v2,$v3,$v4,$v5) = split/\t/;
	next unless(exists $hash_sample{$sample1} and exists $hash_sample{$sample2});
	next if(exists $hash_pair{$sample2}{$sample1});
	$hash_pair{$sample1}{$sample2}{picked} = 1;
	my $g1 = $hash_sample{$sample1}{group};
	my $g2 = $hash_sample{$sample2}{group};
	my @arr = sort ($g1,$g2);
	my $pair = join("_",@arr);
	if($g1 eq $g2){
		$pair = $g1;
	}
	print "$_\t$pair\n";
}
close IN;
