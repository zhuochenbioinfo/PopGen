# by Zhuo Chen, IGDB, CAS
# email1: zhuochen@genetics.ac.cn
# email2: zhuochenbioinfo@gmail.com

use strict;
use warnings;
use Getopt::Long;

my($vcf,$sampleinfo,$maf,$outfile,$mincov,$pass,$minqual);
my $usage = "USAGE:\nperl $0 --vcf <vcf> --smp <sample info list> --out <outfile> --maf <maf> --cov <min cov> --qual <min qual>\n";
$usage .= "<vcf> is the input vcf file. [Necessary]\n";
$usage .= "<sample info list>: #SAMPLE\tGROUP [Necessary]\n";
$usage .= "<outfile> is the output file. [Necessary]\n";
$usage .= "<maf> minor allele frequency. [Optional]\n";
$usage .= "<min qual>: minimum QUAL score. [Optional]\n";

GetOptions(
	"vcf=s" => \$vcf,
	"smp=s" => \$sampleinfo,
	"out=s" => \$outfile,
	"maf=s" => \$maf,
	"cov=s" => \$mincov,
	"pass!" => \$pass,
	"qual=s" => \$minqual,
) or die $usage;

die $usage unless(defined $vcf and defined $sampleinfo and defined $outfile);

my %hash_sample;
my %hash_group;

open(IN,"<$sampleinfo") or die $!;
while(<IN>){
	chomp;
	my($sample,$subtype) = split/\t/;
	$hash_sample{$sample}{group} = $subtype;
	push @{$hash_group{$subtype}{samples}}, $sample;
}
close IN;

my @groups = sort keys %hash_group;
my $samplenum = keys %hash_sample;
foreach my $subtype(@groups){
	my $subtypecount = @{$hash_group{$subtype}{samples}};
	my $ratio = $subtypecount/$samplenum;
	$hash_group{$subtype}{ratio} = $ratio;
}

print "# Reading VCF...\n";

open(VCF,"<$vcf") or die $!;
open(OUT,">$outfile");

my @samples = ();
my @subtypes = ();
my %hash_tmp;

while(<VCF>){
	chomp;
	# pick sample names
	if($_ =~ /^#CHROM/){
		my(undef,undef,undef,undef,undef,undef,undef,undef,undef,@samplenames) = split/\t/;
		for(my $i = 0; $i < @samplenames; $i++){
			my $sample = $samplenames[$i];
			$sample =~ s/Sample_//;
			$sample =~ s/\.variant\d+//;
			$sample =~ s/_1\.fastq\.gz//;
			$hash_sample{$sample}{rank} = $i;
			if(exists $hash_sample{$sample}{group}){
				push @samples, $sample;
				my $subtype = $hash_sample{$sample}{group};
				$hash_tmp{$subtype}++;
			}
		}
		@subtypes = sort keys %hash_tmp;
		print OUT "#CHROM\tPOS\tREF\tBASE\t".join("\t",@subtypes)."\n";
	}
	next if($_ =~ /^#/);
	
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@datas) = split/\t/;
	
	if(defined $pass){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}
	if(defined $minqual){
		next if($qual < $minqual);
	}
	
	my @alts = split/,/,$alt;
	my @bases = ($ref, @alts);
	
	# check if the variant is an indel
	my $checkindel = 0;
	foreach(@bases){
		if(length($_) > 1){
			$checkindel = 1;
		}
	}
	next if($checkindel == 1);
	
	# check if the variant exists in the selected samples
	my $covcount = 0;
	my %hash_gt;
	
	foreach my $sample(@samples){
		my $subtype = "-";
		if(exists $hash_sample{$sample}{group}){
			$subtype = $hash_sample{$sample}{group};
		}
		my $i = $hash_sample{$sample}{rank};
		my $spot = $datas[$i];
		my $gt = "-";
		if($spot =~ /^(\d+)\/(\d+)/){
			$gt = $1;
			if($1 ne $2){
				$gt = "h";
			}
		}
		if($gt ne '-' and $gt ne "h"){
			$covcount++;
		}
		$hash_gt{$gt}{subtype}{$subtype} ++;
		$hash_gt{$gt}{all} ++;
	}
	
	foreach my $gt(sort keys %hash_gt){
		next if($gt eq "-" or $gt eq "h");
		if(defined $mincov){
			next unless($covcount/@samples >= $mincov);
		}
		if(defined $maf){
			next unless($hash_gt{$gt}{all}/@samples > $maf and (@samples - $hash_gt{$gt}{all})/@samples > $maf);
		}
		my $base = $bases[$gt];
		
		my @ratios = ();
		foreach my $subtype(@subtypes){
			unless(exists $hash_gt{$gt}{subtype}{$subtype}){
				$hash_gt{$gt}{subtype}{$subtype} = 0;
			}
			my $ratio = $hash_gt{$gt}{subtype}{$subtype}/$hash_tmp{$subtype};
			push @ratios, $ratio;
		}
		print OUT "$chr\t$pos\t$ref\t$base\t".join("\t",@ratios)."\n";
	}
	
}
close VCF;
close OUT;
