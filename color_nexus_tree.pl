use strict;
use warnings;

# add color to NEXUS tree
my($in,$out,$list) = @ARGV;

my $usage = "Add color to NEXUS tree.\n";
$usage .= "USAGE:\nperl $0 <input nexus tree> <output> <sample color list>\n";
$usage .= "<input nexus tree> is a tree of NEXUS format with no color annotation.\n";
$usage .= "<output> is the name of output tree.\n";
$usage .= "<sample color list> is a tab-delimited list of samples and their corresponding color code(someting like: #0033ff).\n";

=sample list format
Sample_Y220	#0033ff
Sample_449	#0033ff
Sample_Y248	#0033ff
Sample_329	#ff0000
Sample_L22	#009900
Sample_YCX391	#cc00cc
Sample_Y245	#00cccc
=cut

die $usage unless(@ARGV == 3);

my %hash_sample;
open(IN,"<$list") or die $!;
while(<IN>){
	chomp;
	my($sample,$colorCode) = split/\t/;
	$hash_sample{$sample}{color} = $colorCode;
}
close IN;

open(IN,"<$in") or die $!;
open(OUT,">$out");

my $trigger = 0;

while(<IN>){
	chomp;
	if($_ =~ /^begin/){
		$trigger = 1;
		print OUT $_."\n";
		next;
	}
	if($_ =~ /^end/){
		$trigger = 0;
	}
	if($trigger == 0){
		print OUT $_."\n";
		next;
	}
	my @datas = split/,/,$_;
	my @outs;
	for(my $i = 0; $i < @datas; $i++){
		my $data = $datas[$i];
		my $sample = getNameNexus($data);
		my $out = $data;
		if(exists $hash_sample{$sample}){
			my $colorCode = $hash_sample{$sample}{color};
			$out = addColor($data,$colorCode);
		}
		push @outs, $out;
	}
	print OUT join(",",@outs)."\n";
}

sub getNameNexus{
	my($data) = @_;
	my $sample;
	if($data =~ /\(/){
		($sample) = $data =~ /\S*\((\S+?)\:/;
	}else{
		($sample) = $data =~ /(\S+?)\:/;
	}
	return($sample);
}

sub addColor{
	my($data,$colorCode) = @_;
	my $sample;
	if($data =~ /\(/){
		($sample) = $data =~ /\S*\((\S+?)\:/;
	}else{
		($sample) = $data =~ /(\S+?)\:/;
	}
	$data =~ s/$sample/$sample\[\&\!color\=$colorCode\]/;
	$data =~ s/\)/\)\[\&\!color\=$colorCode\]/g;
	return($data);
}
