use strict;
use warnings;
use Getopt::Long;

my($vcf,$out,$cov,$maf,$pass,$pass1,$keepList,$regList,$reg,$misr);
my $usage = "\nTransform vcf to fasta format sequence.\nBy Zhuo Chen\nEmail: chenomics\@163.com\n\n";
$usage .= "Notice: the program picks only bi-allelic SNPs.\n";
$usage .= "USAGE:\nperl $0 --in <input vcf> --out <output prefix>\n";
$usage .= "<input vcf> is the input vcf file. [Necessary]\n";
$usage .= "<output> is the name of output fasta file. [Necessary]\n";
$usage .= "FILTERING OPTIONS:\n";
$usage .= "--cov <minium coverage ratio>, [0,1], Default=0\n";
$usage .= "--maf <minium allele frequency>, [0,1], Default=0\n";
$usage .= "--keep <keep sample list>\n";
$usage .= "--misr <sample genotype missing rate>, [0,1], Default=1; throw sample with high genotype missing rate\n";
$usage .= "--bed <keep region list> in BED format, shall not be set with --reg\n";
$usage .= "--reg <region>, format: chr:start-end, shall not be set with --bed\n";
$usage .= "--pass to keep SNPs with filter TAG [PASS] or [SnpCluster].\n";
$usage .= "--pass1 to keep SNPs with filter TAG [PASS].\n\n";

GetOptions(
	"in=s" => \$vcf,
	"out=s" => \$out,
	"cov=s" => \$cov,
	"misr=s" => \$misr,
	"maf=s" => \$maf,
	"pass!" => \$pass,
	"pass1!" => \$pass1,
	"keep=s" => \$keepList,
	"bed=s" => \$regList,
	"reg=s" => \$reg,
) or die $usage;

die $usage unless(defined $vcf and defined $out);

unless(defined $cov){
	$cov = 0;
}
unless(defined $maf){
	$maf = 0;
}
unless(defined $misr){
	$misr = 1;
}

my %hash_keep;

if(defined $keepList){
	open(IN,"<$keepList") or die $!;
	while(<IN>){
		chomp;
		my($sample,$other) = split/\t/;
		$hash_keep{$sample} = "";
	}
	close IN;
}

my %hash_bed;
my %hash_chr; # This hash is used to check the end of the region list and stop reading vcf

if(defined $regList){
	die $usage if(defined $reg);
	open(IN,"<$regList") or die $!;
	while(<IN>){
		chomp;
		my($chr,$start,$end,$other) = split/\t/;
		next if(exists $hash_bed{$chr}{$start} and $hash_bed{$chr}{$start} > $end);
		$hash_bed{$chr}{$start} = $end;
	}
	close IN;
	# merge regions
	foreach my $chr(keys %hash_bed){
		$hash_chr{$chr} = "";
		my @starts = sort {$a <=> $b} keys %{$hash_bed{$chr}};
		for(my $i = 0; $i < @starts; $i++){
			my $starti = $starts[$i];
			my $endi = $hash_bed{$chr}{$starti};
			for(my $j = $i+1; $j < @starts; $j++){
				my $startj = $starts[$j];
				my $endj = $hash_bed{$chr}{$startj};
				last if($startj > $endi);
				if($endj > $endi){
					$endi = $endj;
					$hash_bed{$chr}{$starti} = $endj;
				}
				delete($hash_bed{$chr}{$startj});
				splice(@starts,$j,1);
				$j--;
			}
		}
	}
}elsif(defined $reg){
	my($chr,$start,$end) = $reg =~ /(\S+):(\d+)\-(\d+)/;
	die $usage unless(defined $end);
	$hash_chr{$chr} = "";
	$hash_bed{$chr}{$start} = $end;
}


open(IN,"<$vcf") or die $!;

my @allSamples = ();
my @keepRanks = ();
my @keepSamples = ();
my @remained_chrs = ();
my $chr_tmp = "";
my @regions = ();

my %hash_seq;

while(<IN>){
	chomp;
	next if($_ =~ /^##/);
	if($_ =~ /^#CHROM/){
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@allSamples) = split/\t/;
		if(defined $keepList){
			for(my $i = 0; $i < @allSamples; $i++){
				my $sample = $allSamples[$i];
				next unless(exists $hash_keep{$sample});
				push @keepRanks, $i;
				push @keepSamples, $sample;
			}
		}else{
			@keepRanks = (0..@allSamples-1);
			@keepSamples = @allSamples;
		}
		if(@keepRanks == 0){
			die "#ERROR: No sample remain after filtering, please check your files!\n";
		}else{
			print "# Remain ".@keepRanks." samples after filtering.\n";
		}
		next;
	}
	
	#last if($. == 5000);

	# split line
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$datas_join) = split/\t/,$_,10;

#=this function is under construction
	
	# Check if the position locates in the candidate regions
	goto NOBED unless(defined $regList or defined $reg);
	
	next unless(exists $hash_bed{$chr});
	if($chr ne $chr_tmp){
		@regions = ();
		foreach my $start(sort {$a <=> $b} keys %{$hash_bed{$chr}}){
			my $end = $hash_bed{$chr}{$start};	
			push @regions, "$start,$end";
		}
		$chr_tmp = $chr;
		if(exists $hash_chr{$chr}){
			delete($hash_chr{$chr});
		}
		@remained_chrs = keys %hash_chr;
		print "\t# Reading chr:$chr\n";
		if(@regions == 0){
			print "\t# Skipping chr:$chr\n";
		}
	}
	
	if(@remained_chrs == 0 and @regions == 0){
		last;
	}

	if(@regions == 0){
		next;
	}
	
	my $pickIt = 0;
	for(my $i = 0; $i < @regions; $i++){
		my($start,$end) = split/,/,$regions[$i];
		if($pos < $start){
			last;
		}elsif($pos > $end){
			splice(@regions, $i, 1);
			$i--;
			next;
		}else{
			$pickIt = 1;
			last;
		}
	}

	next unless($pickIt == 1);
	NOBED:

#=cut

	# Check if the variant is a bi-allelic SNP
	next unless(length($ref) == 1 and length($alt) == 1);
	
	# Check if the SNP pass the filtration
	if(defined $pass){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}
	if(defined $pass1){
		next unless($filter eq "PASS");
	}
	
	# Check maf and coverage
	my @datas = split/\t/, $datas_join;

	my @keepDatas = ();
	foreach my $rank(@keepRanks){
		push @keepDatas, $datas[$rank];
	}
	my($cov_,$maf_) = get_cov_maf(@keepDatas);
	next unless($cov_ > $cov and $maf_ > $maf);

	# push bases to fasta
	foreach my $rank(@keepRanks){
		my $data = $datas[$rank];
		my $sample = $allSamples[$rank];
		my $base = vcf2base($data,$ref,$alt);
		$hash_seq{$sample}{fasta} .= $base;
	}	
}
close IN;


open(OUT,">$out");
for(my $i = 0; $i < @keepSamples; $i++){
	my $sample = $keepSamples[$i];
	my $seq = $hash_seq{$sample}{fasta};
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

sub vcf2base{
	my($data,@bases) = @_;
	my $base = "-";
	if($data =~ /(\d+)\/(\d+)/){
		if($1 == $2){
			$base = $bases[$1];
		}
	}
	return($base);
}


sub get_cov_maf{
	my @genos = @_;
	my $altcount = 0;
	my $covcount = 0;
	my $lostcount = 0;
	
	foreach my $geno(@genos){
		my $gt = "-";
		if($geno =~ /^(\d+)\/(\d+)/){
			$gt = $1;
			if($1 ne $2){
				$gt = "h";
				$altcount ++;
			}elsif($1 != 0){
				$altcount += 2;
			}
		}
		if($gt ne '-'){
			$covcount += 2;
		}else{
			$lostcount += 2;
		}
	}
	# debug 20171117
	my $cov = 1 - ($lostcount/@genos)/2;
	my $maf = ($altcount/@genos)/2;
	my $maf2 = $cov - $maf;
	if($maf > $maf2){
		$maf = $maf2;
	}
	return($cov,$maf);
}

sub makefasta{
	my $seq = $_[0];
	$seq =~ s/(.{50})/$1\n/g;
	unless($seq =~ /\n$/){
		$seq .= "\n";
	}
	return($seq);
}

sub missingRatio{
	my($seq) = @_;
	my $count = ($seq =~ s/\-/#/g);
	my $len = length($seq);
	my $ratio = $count/$len;
	return($ratio);
}
