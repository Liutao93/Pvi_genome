#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path getcwd);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Config::General;

use File::Path qw(make_path remove_tree);
use File::Copy;
use File::Spec;
my $BEGIN_TIME=time();
my $version="1.0";
use PerlIO::gzip;
use Bio::SeqIO;

###################

my %opts;
GetOptions(\%opts,"in=s","anno=s" ,"o=s", "h" );

#&help()if(defined $opts{h});
if(!defined($opts{in}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";


	Usage:

		-in          genome.tblout.final.xls file <infile>                                 must be given
		-o           output file <infile>                                                  must be given
		-h           Help document

	example:
		perl /share/work/wangq/piplines/bin/genome_assembly/nocoding/Rfam_anno.pl 
		-in /share/work/wangq/piplines/bin/genome_assembly/nocoding/ncRNA.final.gff 
		-o /share/work/wangq/piplines/bin/genome_assembly/nocoding/noncoding.stat.xls

	Usage End.


	exit;
}





my%stat;

open (IN, "$opts{in}") || die "can't open $opts{in}, $!\n" ;
open (OUT, ">$opts{o}") || die "can't open $opts{o}, $!\n" ;
while(<IN>){
	chomp;
	next if(/^#/);
	my@line = split("\t",$_);
	next if($line[2] eq "gene");
	my$len =abs($line[4] - $line[3]) + 1 ;

	
	if($line[2] eq "tRNA"){
		$stat{"tRNA"}[0] ++;
		$stat{"tRNA"}[1] += $len;
		
	}elsif($line[2] eq "miRNA"){
	  $stat{"miRNA"}[0] ++;
		$stat{"miRNA"}[1] += $len;
		
	}elsif($line[2] eq "rRNA"){
		$stat{"rRNA"}[0] ++;
		$stat{"rRNA"}[1] += $len;
		
		my($subtype) = $line[8] =~ /ID=(.+?)_rRNA/;
		
		$stat{$subtype}[0] ++;
		$stat{$subtype}[1] += $len;
		
		
	}elsif($line[2] eq "snRNA" || $line[2] eq "snoRNA"){
		$stat{"snRNA"}[0] ++;
		$stat{"snRNA"}[1] += $len;
		
		my($subtype) = $line[8] =~ /subtype=(.*)/;
		
		$stat{$subtype}[0] ++;
		$stat{$subtype}[1] += $len;
		
		
	}
}
close(IN);


my$p;
while(my($key,$value) = each %stat){
	$p = $stat{$key}[1] / $stat{$key}[0] ;
	$stat{$key}[2] = sprintf("%.2f",$p);
}

my@order = ("miRNA","tRNA", "rRNA", "5s", "8s", "18s", "28s", "snRNA", "CD-box", "HACA-box", "splicing", "scaRNA");
print OUT "Type\tClass\tNumber\tTotal length(bp)\tAverage length(bp)\n";
my$a;
foreach my$i (@order){
	if(!exists $stat{$i}){
		#$stat{$i} = [0,0,0];
		next;	
	}
	if($i eq "miRNA"){
	  $a = "miRNA\tmiRNA";
	}elsif($i eq "tRNA"){
		$a = "tRNA\ttRNA";
	}elsif($i eq "rRNA"){
		$a = "rRNA\trRNA";
	}elsif($i eq "snRNA"){
		$a = "snRNA\tsnRNA";
	}else{
		$a = " \t$i";
	}

	$stat{$i}[0] =~ s/(?<=\d)(?=(\d{3})+\b)/,/g;
	$stat{$i}[1] =~ s/(?<=\d)(?=(\d{3})+\b)/,/g;
	$stat{$i}[2] =~ s/(?<=\d)(?=(\d{3})+\b)/,/g;


	print OUT join("\t",$a,$stat{$i}[0],$stat{$i}[1],$stat{$i}[2])."\n";
}
close(OUT);


	

