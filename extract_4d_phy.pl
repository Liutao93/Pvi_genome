#!/usr/bin/env perl
use strict;

die "perl $0 <input.phy>  <out_4d.phy>\n"  unless @ARGV==2;

my $file = shift;
my $phy = shift;

my %codons = (
	'CTT'=>'L', 'CTC'=>'L', 'CTA'=>'L', 'CTG'=>'L',
	'GTT'=>'V', 'GTC'=>'V', 'GTA'=>'V', 'GTG'=>'V',
	'TCT'=>'S', 'TCC'=>'S', 'TCA'=>'S', 'TCG'=>'S',
	'CCT'=>'P', 'CCC'=>'P', 'CCA'=>'P', 'CCG'=>'P',
	'ACT'=>'T', 'ACC'=>'T', 'ACA'=>'T', 'ACG'=>'T',
	'GCT'=>'A', 'GCC'=>'A', 'GCA'=>'A', 'GCG'=>'A',
	'CGT'=>'R', 'CGC'=>'R', 'CGA'=>'R', 'CGG'=>'R',
	'GGT'=>'G', 'GGC'=>'G', 'GGA'=>'G', 'GGG'=>'G'
);


my $i=0;

my ($num_species,$length_seq);
my (@seq,@name);

open IN,$file or die "fail to open $file\n";
while(<IN>){
        chomp;
	next if (/^\s*$/);
        if(/(\d+)\s+(\d+)/){
                $num_species = $1;
                $length_seq = $2;
                next;
        }
        my @temp=split;
        $seq[$i] = $temp[1];
        $name[$i] = $temp[0];
        $i++;
}
close IN;

die "error: check seq number in header\n" unless ($i   ==  $num_species);

my @out;
for(my $j = 0; $j < $length_seq; $j += 3){
        my @codon = ();
        my @site = ();
        my @first2 = ();
        my $flag = 1;
        for(my $i = 0; $i < $num_species ; $i++){
                $codon[$i] = uc(substr($seq[$i],$j,3));
                $site[$i] = uc(substr($seq[$i],$j+2,1));
                $first2[$i] = uc(substr($seq[$i],$j,2));
                if( $i > 0 and $first2[$i] ne $first2[$i-1]){
                        $flag = 0;
                        last;
                }
                if(! exists $codons{$codon[$i]}){
                        $flag = 0;
                        last;
                }


        }
        if($flag == 1){
                for( $i = 0; $i < $num_species; $i++){
                        $out[$i] .= $site[$i];
                }
        }

}



open OUT, ">$phy" or die $! ;

my $length = length($out[0]);

print OUT " $num_species  $length\n";
for( my $i = 0; $i < $num_species; $i++){
        printf OUT "$name[$i]  $out[$i]\n";
}
close OUT;

