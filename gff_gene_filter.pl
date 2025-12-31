#!/usr/bin/perl
use strict;
use List::Util qw(sum);
use Data::Dumper;
use Cwd qw(abs_path getcwd);
my $usage = <<USAGE;
Usage:
    perl $0 all.gff3  Min_exon_number_of*s_length Min_cds_ratio_vs_cDNA Max_intron_length_ratio Max_cds_length_ratio
For Example:
    perl $0 all.gff3  2 150  > final.gff3
    1. 按照先后顺序对assemblies.fasta.transdecoder.genome.gff3中的基因模型进行过滤。
    2. 再过滤cds数目<3的基因模型。
    3. 再过滤cds长度<900的基因模型。
    4. 再过滤cds region占exon region比例<0.6的基因模型。
    5. 再过滤含有过长intron的基因模型。根据全部数据的cds(而不是exon)信息得到的intron长度信息；再将intron长度按从小到大进行排序，选取95%分位数为阈值。
    6. 最后过滤cds长度过长的基因模型。根据全部数据得到的cds长度信息；将cds长度长度从小到大排序，选取95%分为数为阈值。
USAGE
if (@ARGV==0){die $usage}


open IN, $ARGV[0] or die $!;
my (@intron_length, @cds_length, %validate);
my ($gene_number_total, $gene_number_validate,  $gene_number_exonLack,$gene_number_short)=(0,0,0,0);
my $cds_length_threshold=$ARGV[2];

my%gff;
my%mRNA2gene=();
while(<IN>){
        chomp;
        my@a=split(/\t/);
#	print join(",",@a)."\n";
#        if($a[2] eq "gene" ){
#		
 #               my($geneid)=($a[8]=~ m/ID=([^;]*)/);
#		$gff{$geneid}{chr}=$a[0];
#		$gff{$geneid}{start}=$a[3];
#		$gff{$geneid}{end}=$a[4];
#		
#	}
        if($a[2] eq "mRNA" || $a[2] eq "transcript"){
	#	print $_; 
               my($mRNAid)=($a[8]=~ m/ID=([^;]*)/);
                my($geneid)=($a[8]=~ m/Parent=([^;]*)/);
                $mRNA2gene{$mRNAid}=$geneid;
		$gff{$geneid}{$mRNAid}{"exon"}=[];
		$gff{$geneid}{$mRNAid}{"CDS"}=[];
		$gff{$geneid}{$mRNAid}{"chr"}=$a[0];
		$gff{$geneid}{$mRNAid}{"start"}=$a[3];
		$gff{$geneid}{$mRNAid}{"end"}=$a[4];
		$gff{$geneid}{$mRNAid}{"strand"}=$a[6];
		 
       }
	if($a[2] eq "CDS"){
		 my($p)=($a[8]=~ m/Parent=([^;]*)/);
		push   @{$gff{$mRNA2gene{$p}}{$p}{"CDS"}},"CDS\t$a[3]\t$a[4]";   
		
	}

        if($a[2] eq "exon"){
                 my($p)=($a[8]=~ m/Parent=([^;]*)/);
		push   @{$gff{$mRNA2gene{$p}}{$p}{"exon"}},"exon\t$a[3]\t$a[4]";   

        }
}
close IN;


#一个基因保留最长的转录本用于统计分析
#print Dumper(\%gff);
my%gff1=();
my %longest=();

for my $g (keys %gff){
  
  for my $m (keys %{$gff{$g}}){
	next if ($m eq "chr");
	next if ($m eq "end");
	next if ($m eq "start");

        my @cds_line = @{$gff{$g}{$m}{"CDS"}};
        my @exon_line = @{$gff{$g}{$m}{"exon"}};
	
        #my @intron_len=&get_intron_length(@cds_line);

        my ($cds_length, $exon_length);
        foreach (@cds_line) {
            m/CDS\t(\d+)\t(\d+)/;
            $cds_length += ($2 - $1 + 1);
        }
        #push @cds_length, $cds_length;
        foreach (@exon_line) {
            m/exon\t(\d+)\t(\d+)/;
            $exon_length += ($2 - $1 + 1);
        }
	
	if(exists $gff1{$g} and $cds_length >$gff1{$g}{cds_len}){
                $gff1{$g}{exon_len}= $exon_length;
                $gff1{$g}{exon_num}= \@exon_line;
                $gff1{$g}{cds_num}= \@cds_line;
                $gff1{$g}{cds_len}= $cds_length;
		$gff1{$g}{chr}=$gff{$g}{$m}{chr};
		$gff1{$g}{start}=$gff{$g}{$m}{start};
		$gff1{$g}{end}=$gff{$g}{$m}{end};
		$gff1{$g}{strand}=$gff{$g}{$m}{strand};
                #$gff1{$g}{intron_len}=$intron_len ;
		$longest{$g}=$m;
	}else{

	
        	$gff1{$g}{exon_len}= $exon_length;
                $gff1{$g}{exon_num}= \@exon_line;
        	$gff1{$g}{cds_num}= \@cds_line;
        	$gff1{$g}{cds_len}= $cds_length;
        	#$gff1{$g}{intron_len}=$intron_len ;
		$gff1{$g}{chr}=$gff{$g}{$m}{chr};
		$gff1{$g}{start}=$gff{$g}{$m}{start};
		$gff1{$g}{end}=$gff{$g}{$m}{end};
		$gff1{$g}{strand}=$gff{$g}{$m}{strand};
		$longest{$g}=$m;
	}
  }
}

#汇总统计结果
#print Dumper(\%gff1);
my %keepID=();

my $wd=getcwd;
my $keepfile="$wd/keeped.gene.info.tsv";
open OUT ,">$keepfile" or die "$! $keepfile\n";
print OUT "geneID\tmRNAID\tchr\tstart\tend\tstrand\tcds_lengths\tcds_length_sum\t\tcds_num\texon_lengths\texon_length_sum\texon_num\tintron_lengths\tintron_length_sum\tintron_num\n";
for my $g ( sort {$gff1{$a}{chr} cmp $gff1{$b}{chr} || $gff1{$a}{start} <=> $gff1{$b}{start} }  keys %gff1) {
	$gene_number_total++;
	my @cds_region=@{$gff1{$g}{cds_num}};
	my @intron_length = &get_intron_length(@cds_region);
	
	@intron_length = sort {$b <=> $a} @intron_length;
	my $cds_keep=0;
	my $exon_keep=0;
	my@exon_lengths=&get_length(@{$gff1{$g}{exon_num}});
	my@cds_lengths=&get_length(@{$gff1{$g}{cds_num}});
	my$cds_length=sum(@cds_lengths);
        
	if ($cds_length <= $cds_length_threshold) {
            $gene_number_short++;
                
        }else{
		$cds_keep=1;
	}
         
	
        if (@{$gff1{$g}{exon_num}} < $ARGV[1]){
            $gene_number_exonLack++;
        }else{
		$exon_keep=1;	
	}
	
	if($cds_keep and $exon_keep){
		$gene_number_validate++;
		print OUT "$g\t$longest{$g}\t$gff1{$g}{chr}\t$gff1{$g}{start}\t$gff1{$g}{end}\t".
				join(",",@cds_lengths)."\t"."$gff1{$g}{cds_len}\t".@cds_lengths."\t".
				join(",",@exon_lengths)."\t"."$gff1{$g}{exon_len}\t".@exon_lengths."\t".
				join(",",@intron_length)."\t".sum(@intron_length)."\t".@intron_length."\n";
		$keepID{$g}=1;
		$keepID{$longest{$g}}=1;	
	}

}
close OUT;


print STDERR "
Total gene number:                  $gene_number_total
Keeped gene number:               $gene_number_validate
Exon < $ARGV[1] gene number:               $gene_number_exonLack
CDS length < $ARGV[2] gene number:       $gene_number_short \n";
#######################################
print STDERR "save gene info: $keepfile\n";

open IN, $ARGV[0] or die $!;

while(<IN>){
        chomp;
	next if /^\s*$/;
        my@a=split(/\t/);
	if($a[8]){
		if($a[8]=~ m/ID=([^;]*)/ and exists $keepID{$1}){
			print "$_\n";
			next;
		}

		if($a[8]=~ m/Parent=([^;]*)/ and exists $keepID{$1}){
			next if $a[2] eq "mRNA" or $a[2] eq "transcript";
			print "$_\n";
			next;
		}
	}else{
		print "$_\n"
	}	 
}


close(IN);




#######################################
sub get_intron_length {
    my @exon;
    foreach (@_) {
        m/(\d+)\t(\d+)$/;
        push @exon, "$1\t$2";
    }
    @exon = sort {$a <=> $b} @exon;
    my $first = shift @exon;
    my $end = $1 if $first =~ m/(\d+)$/;
    my @intron_length;
    foreach (@exon) {
        m/^(\d+)\t(\d+)$/;
        my $intron_length = $1 - $end-1;
        $end = $2;
        push @intron_length, $intron_length;
    }
    return @intron_length;
}

sub get_length{
	my @region;
	foreach(@_){
        	m/(\d+)\t(\d+)$/;
        	push @region, abs($1-$2+1);
		
	}
	

	return @region;

}


