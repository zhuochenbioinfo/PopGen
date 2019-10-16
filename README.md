# PopGen
Small scripts for population genetics analysis. 

Most of the scripts contain detailed usage in the codes.

## SNP_group_ratio.pl
Calculate SNP allele frequency for selected groups of samples.

## color_nexus_tree.pl
Add color annotation to samples in a nexus format tree. The colors are added to the branches.

## dissect_admixture_Q.pl
Dissect the output Q file of [ADMIXTURE](http://software.genetics.ucla.edu/admixture/) population structure inference.

## eigenstartgeno_to_fasta.pl
Transform eigenstart geno file to fasta sequence file.

## eigenstartgeno_to_fasta.withRefSeq.pl
Transform eigenstart geno file to fasta sequence file, with a pseudo reference sequence.

## hapmap_to_eigenstartgeno.pl
Transform hapmap format file to eigenstart geno file. May not adapt to all types of hapmap file.

## hetPi.pl
Nucleotide diversity parameter Pi, with heterozygous genotype treated as an independent genotype.

## vcf2treemix.pl
Transform VCF file to group allele frequency table. Prepare input file for [Treemix](https://bitbucket.org/nygcresearch/treemix/wiki/Home). TreeMix is a method for inferring the patterns of population splits and mixtures in the history of a set of populations.

## vcf_to_eigenstartgeno.pl
Transform VCF file to eigenstart geno format. Typically for [ADMIXTURE](http://software.genetics.ucla.edu/admixture/) inputs.

## vcf_to_fasta.pl
Pick SNPs from VCF file and transform to fasta format sequences for selected samples.

