#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://study.163.com/course/introduction/1006448023.htm?share=1&shareId=1030291076
#perl高级：https://study.163.com/course/introduction/1004833023.htm?share=1&shareId=1030291076

use strict ;

die "perl $0 <in.ks.tsv>  <in.collinearity>  > out.ks.tsv\n" unless @ARGV==2;
my $ks = shift;
my $syn = shift;


my %syn;
open IN, $syn or die $!;
while(<IN>){
	chomp;
	next if (/^#/) ;
	#276-  0:        SIN_1010716     SIN_1024252       3e-82
	if($_ =~ /^\s*\d+-\s*\d+:\s*(\S+)\s+(\S+)\s+\S+/){
		my $id1 = $1;
		my $id2 = $2;
		$syn{$id2}{$id1} = 1;
		$syn{$id1}{$id2} = 1;
	}
}
close IN;

open KS, $ks or die $!;
my $hd = <KS>;
print $hd;
while(<KS>){
	chomp;
	my @t = split;
	#SIN_1014290__SIN_1023532   0.94118 0.84115 2448.0  2223.0  0.12698 GF_002203
	my @id = split /__/ ,$t[0];
	my $id1 = $id[0];
	my $id2 = $id[1];
	unless(!$syn{$id1}{$id2} && !$syn{$id2}{$id1}){
		print "$_\n";
	}
}
close KS;
