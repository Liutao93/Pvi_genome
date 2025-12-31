#!/usr/bin/env perl

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory


=head1 Name

repeat_to_gff.pl  --  convert repeat raw formats to gff format

=head1 Description

This program can read repeatmasker .out file, repeatproteinmask .annot file, or trf .dat file,
and convert it to gff format.

=head1 Version

  Author: huangls@biomics.com.cn
  Version: 1.0,  Date: 2022.12.07
  Note: suffix of input flie name  must be  .dat, .out or .annot
    
=head1 Usage

  --in <str>  set input file path
  --out <str>  set output file name
  --prefix <str>  set a prefix before repeat element ID
  --verbose   output running progress information to screen  
  --help      output help information to screen  
  
=head1 Example

 perl $Bin/repeat_to_gff.pl -in ./trf.dat -out trf.gff
 perl $Bin/repeat_to_gff.pl -in ./RepeatMasker.out -out RepeatMasker.gff
 perl $Bin/repeat_to_gff.pl -in ./Proteinmask.annot -out Proteinmask.gff
 
=cut


##get options from command line into variables and set default values
my ($Verbose,$Help,$Prefix,$in,$out);
GetOptions(
	"prefix:s"=>\$Prefix,
	"in:s"=>\$in,
	"out:s"=>\$out,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if ( !$in || $Help);

my $repeat_file = $in;
my $outfile=$out;
open OUT,">$outfile" || die "fail creat $outfile";
#print OUT "##gff-version 3\n";

dat_to_gff3($repeat_file,$Prefix) if($repeat_file =~ /\.dat$/);
out_to_gff3($repeat_file,$Prefix) if($repeat_file =~ /\.out$/);
annot_to_gff3($repeat_file,$Prefix) if($repeat_file =~ /\.annot$/);
close(OUT);
####################################################
################### Sub Routines ###################
####################################################


##facilitate to creat marks
####################################################
sub create_marker {
	my $number = shift || 100000;
	$number =~ s/\d/0/g;
	$number++;
	return $number;
}


##Start	End	PeriodSize 	CopyNumber	ConsensusSize	PercentMatches	PercentIndels	Score	A	C	G	T	Entropy(0-2)	consensus	repeatSequences
##19670039 19670073 4     8.8         4              93                0              61    22  0   28 48     1.51           TGTA       TGTATGTATGTATGTATGTATGTAGGTATGTATGT
####################################################
sub dat_to_gff3 {
	my $file = shift;
	my $pre_tag = shift;
	my $output;
	$pre_tag .= "_" if($pre_tag); 
	my $line_num = `wc -l $file`;
	$line_num = $1 if($line_num =~ /^(\d+)/);
	my $mark = create_marker($line_num);
	my $chr;
	open IN,$file || die "fail open $file";
	while (<IN>) {
		$chr = $1 if(/^Sequence:\s+(\S+)/);
		my @t = split /\s+/;
		next if(@t != 15);
		my $start = $t[0];
		my $end = $t[1];
		my $id = $pre_tag."TR".$mark;
		
		my $score = $t[7];
		my $strand = "+";
		$output = "$chr\tTRF\tTandemRepeat\t$start\t$end\t$score\t$strand\t.\tID=$id;ConsensusSize=$t[4];CopyNumber=$t[3];PercentMatches=$t[5];PercentIndels=$t[6];Consensus=$t[13];RepeatSeq=$t[14]\n";
		print OUT $output;
		$mark++;
	}
	close IN;


	
	close OUT;	

}



##  SW   perc perc perc  query     position in query              matching       repeat       position in repeat
##score   div. del. ins.  sequence  begin    end          (left)   repeat         class/family begin   end   (left)  ID
#245   35.2  2.5  0.6  chr1       1001400  1001556 (18753082) + TS2            SINE            80   239  (416)  12  
####################################################
sub out_to_gff3 {
	my $file = shift;
	my $pre_tag = shift;
	my $output;
	
	$pre_tag .= "_" if($pre_tag); 
	my $line_num = `wc -l $file`;
	$line_num = $1 if($line_num =~ /^(\d+)/);
	my $mark = create_marker($line_num);

	open IN,$file || die "fail open $file";
	while (<IN>) {
		s/^\s+//;
		my @t = split /\s+/;
		next if($t[0] =~ /\D/ || !$t[0]);
		next if $t[10]=~/Low/i || $t[10] =~/Simple/i  || $t[10]=~ /Satellite/i ; #不统计简单重复序列
		my $start = $t[5];
		my $end = $t[6];
		my $id = $pre_tag."TE".$mark;
		my $chr = $t[4];
		my $score = $t[0];
		my $strand = ($t[8] eq '+') ? "+" : "-";
		my $target = $t[9];
		my $class = $t[10];
		my @ary;
		push @ary,$t[11] if($t[11] !~ /[\(\)]/);
		push @ary,$t[12] if($t[12] !~ /[\(\)]/);
		push @ary,$t[13] if($t[13] !~ /[\(\)]/);
		@ary = sort {$a<=>$b} @ary;
		my ($target_start,$target_end) = ($ary[0],$ary[1]);

		$output = "$chr\tRepeatMasker\tTransposon\t$start\t$end\t$score\t$strand\t.\tID=$id;Target=$target $target_start $target_end;Class=$class;PercDiv=$t[1];PercDel=$t[2];PercIns=$t[3];\n";
		print OUT $output;
		$mark++;
	}
	close IN;


}


##6.00e-22      117 WUBlastX Chr07frag1M       26352    26816    - RETRO1_pol      LTR/Gypsy            360      533
####################################################
sub annot_to_gff3 {
	my $file = shift;
	my $pre_tag = shift;
	my $output;

	$pre_tag .= "_" if($pre_tag); 
	my $line_num = `wc -l $file`;
	$line_num = $1 if($line_num =~ /^(\d+)/);
	my $mark = create_marker($line_num);

	open IN,$file || die "fail open $file";
	while (<IN>) {
		s/^\s+//;
		my @t = split /\s+/;
		next if($t[0] =~ /^pValue/);
		my $start = $t[4];
		my $end = $t[5];
		my $id = $pre_tag."TP".$mark;
		my $chr = $t[3];
		my $score = $t[1]; ##use pvalue here
		my $strand = $t[6];
		my $target = $t[7];
		my $class = $t[8];
		my @ary = ($t[9],$t[10]);
		@ary = sort {$a<=>$b} @ary;
		my ($target_start,$target_end) = ($ary[0],$ary[1]);

		$output = "$chr\tRepeatProteinMask\tTEprotein\t$start\t$end\t$score\t$strand\t.\tID=$id;Target=$target $target_start $target_end;Class=$class;pValue=$t[0];\n";
		print OUT $output;
		$mark++;
	}
	close IN;

}