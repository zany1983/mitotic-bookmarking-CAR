#!/usr/bin/perl
use strict;
use warnings;


# find distal and inter-chromosomal MARGI interactions

my @lines;
my %hash;

my @inFileName = ("/Users/yan/Documents/scientific/CAR/downloaded_data/MARGI_zhongsheng/diMARGI_H9_associations.txt");

my @inFileHandle = ();
my @outFileHandle = ();

for (my $i = 0; $i < scalar @inFileName; $i++) {
	open($inFileHandle[$i], "<", $inFileName[$i]) or die "Couldn't open $inFileName[$i] for appending: $!";
}

#my @outFileName = ("A549_merged_GeneCount.txt",);
my @outFileName = ("diMARGI_H9_distal_associations.txt","diMARGI_H9_distal_RNA_bed.txt","diMARGI_H9_distal_DNA_bed.txt","diMARGI_H9_RNA_bed.txt","diMARGI_H9_DNA_bed.txt");
for (my $i = 0; $i < scalar @outFileName; $i++) {
	open($outFileHandle[$i], ">", $outFileName[$i]) or die "Couldn't open $outFileName[$i] for appending: $!";
}



 readline($inFileHandle[0]);
print {$outFileHandle[1]} "RNA_chr","\t","RNA_start","\t","RNA_end","\t","RNA_name","\t","RNA_strand","\n";
print {$outFileHandle[2]} "DNA_chr","\t","DNA_start","\t","DNA_end","\t","DNA_strand","\n";
print {$outFileHandle[3]} "RNA_chr","\t","RNA_start","\t","RNA_end","\t","RNA_name","\t","RNA_strand","\n";
print {$outFileHandle[4]} "DNA_chr","\t","DNA_start","\t","DNA_end","\t","DNA_strand","\n";

	while(defined( my $line=readline($inFileHandle[0]))){ 
		chomp $line;
	#	$line=~s/\r"//g;
		my @content=split("\t",$line);	
			print {$outFileHandle[3]} $content[0],"\t",$content[1],"\t",$content[2],"\t",$content[4],"\t",$content[3],"\n";
			print {$outFileHandle[4]} $content[5],"\t",$content[6],"\t",$content[7],"\t",$content[8],"\n";

		if($content[10] !~/proximal/) {
			print {$outFileHandle[0]} $line,"\n";
			print {$outFileHandle[1]} $content[0],"\t",$content[1],"\t",$content[2],"\t",$content[4],"\t",$content[3],"\n";
			print {$outFileHandle[2]} $content[5],"\t",$content[6],"\t",$content[7],"\t",$content[8],"\n";
		}	
	
	}




