#!/usr/bin/perl

use strict;
use warnings;

# 从含有gene 那么的文件中逐行注释基因的信息，包括染色体位置，基因家族
#注释GeneCode文件是，选择最长转录本来返回基因组位置


#my @inFileName = ("merged_table.txt","GeneCodeCompV19_refformat.txt","HGNC_completeset.txt");
my @inFileName = ("merged_table_defined_CAR.txt","GeneCodeCompV19_refformat.txt","HGNC_completeset.txt");
		
my @inFileHandle = ();

for (my $i = 0; $i < scalar @inFileName; $i++) {
	open($inFileHandle[$i], "<", $inFileName[$i]) or die "Couldn't open $inFileName[$i] for appending: $!";
}

my @outFileName = ("anno.txt",);
		
my @outFileHandle = ();

for (my $i = 0; $i < scalar @outFileName; $i++) {
	open($outFileHandle[$i], ">", $outFileName[$i]) or die "Couldn't open $outFileName[$i] for appending: $!";
}

my(%locusType,%locusFam);


#注释GeneCode文件是，选择最长转录本来返回基因组位置


sub max{
	if ($_[0] > $_[1]){
		return $_[0];
	}
	else{
		return $_[1];
	}
}

sub min {
		if ($_[0] > $_[1]){
		return $_[1];
	}
	else{
		return $_[0];
	}
}


sub start {
	if ($_[0] eq "+"){
		return &min($_[1], $_[2]);
	}
	else {
		return &max ($_[1], $_[2]);
	}
}

sub end {
	if ($_[0] eq "+"){
		return &max ($_[1], $_[2]);
	}
	else {
		return &min ($_[1], $_[2]);
	}
}





my (%Chr_gene,%strand_gene,%txStart_gene,%txEnd_gene);
readline($inFileHandle[1]);
	while (defined(my $line=readline($inFileHandle[1]))) {
		chomp $line;
		$line=~s/\r//;
		my @content=split("\t", $line);
		if(!$Chr_gene{$content[0]}){
		$Chr_gene{$content[0]}=$content[2];
		$strand_gene{$content[0]}=$content[3];
		$txStart_gene{$content[0]}= $content[4];
		$txEnd_gene{$content[0]}=$content[5];
	}

}


# readline($inFileHandle[1]);
# 	while (defined(my $line=readline($inFileHandle[1]))) {
# 		chomp $line;
# 		$line=~s/\r//;
# 		my @content=split("\t", $line);
# 		$Chr_gene{$content[5]}=$content[1];
# 		$strand_gene{$content[5]}=$content[2];
# 		$txStart_gene{$content[5]}=$txStart_gene{$content[5]}? &min($txStart_gene{$content[5]},$content[3]): $content[3];
# 		$txEnd_gene{$content[5]}=$txEnd_gene{$content[5]}? &max($txEnd_gene{$content[5]},$content[4]): $content[4];
# }

readline($inFileHandle[-1]);
	while (defined(my $line=readline($inFileHandle[-1]))) {
		chomp $line;
		$line=~s/\r//;
	my	@content=split("\t",$line);
		$locusType{$content[1]}=$content[4];
		$locusFam{$content[1]}=$content[12];
		
	}
	
my $temp=readline($inFileHandle[0]);
chomp $temp;
print "please specify the column where the geneId is(0 based)\n";
my $i=<>;
chomp $i;


	$temp=~s/\r//;
	print {$outFileHandle[0]} $temp, "\t","Chrom","\t","strand","\t","txStart\t","txEnd\t","gene_type","\t","Gene_fam","\n";
		while (defined(my $line=readline($inFileHandle[0]))) {
		chomp $line;
		$line=~s/[\r"]//;
		my	@content=split("\t",$line);
		#$content[$i]=~s/"//;

		my $Chrom=$Chr_gene{$content[$i]}?$Chr_gene{$content[$i]}:"Nd";
		my $strand=$strand_gene{$content[$i]}?$strand_gene{$content[$i]}:"Nd";
		my $txStart=$txStart_gene{$content[$i]}?$txStart_gene{$content[$i]}:"Nd";
		my $txEnd=$txEnd_gene{$content[$i]}?$txEnd_gene{$content[$i]}:"Nd";
		my $gene_type=$locusType{$content[$i]}?$locusType{$content[$i]}:"Nd";
		my $Gene_fam=$locusFam{$content[$i]}?$locusFam{$content[$i]}:"Nd";


		
		print {$outFileHandle[0]} $line, "\t","$Chrom","\t","$strand","\t","$txStart\t","$txEnd\t","$gene_type","\t","$Gene_fam","\n";
		
	}