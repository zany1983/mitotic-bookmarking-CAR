#!/usr/bin/perl

use strict;
use warnings;

# pick TSS information from gene.bed


#my @inFileName = ("merged_table.txt","GeneCodeCompV19_refformat.txt","HGNC_completeset.txt");
my @inFileName = ("gene_table.bed");
my @inFileHandle = ();

for (my $i = 0; $i < scalar @inFileName; $i++) {
	open($inFileHandle[$i], "<", $inFileName[$i]) or die "Couldn't open $inFileName[$i] for appending: $!";
}

my @outFileName = ("gene_TSS.bed",);
		
my @outFileHandle = ();

for (my $i = 0; $i < scalar @outFileName; $i++) {
	open($outFileHandle[$i], ">", $outFileName[$i]) or die "Couldn't open $outFileName[$i] for appending: $!";
}

my(%locusType,%locusFam);

my $upwindow=500;
my $downwindow=500;


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
my $t=readline($inFileHandle[0]);

print {$outFileHandle[0]} "#TSS_chr","\t","TSS_start","\t","TSS_end","\t",$t;
	while (defined(my $line=readline($inFileHandle[0]))) {
		chomp $line;
		$line=~s/[\r"]//;
		my @content=split("\t", $line);
		if($content[5] eq "+"){
			print {$outFileHandle[0]} $content[0],"\t",$content[1]-$upwindow,"\t",$content[1]+$downwindow,"\t",$line,"\n";
		}
		else{
			print {$outFileHandle[0]} $content[0],"\t",$content[2]-$upwindow,"\t",$content[2]+$downwindow,"\t",$line,"\n";

		}
}
