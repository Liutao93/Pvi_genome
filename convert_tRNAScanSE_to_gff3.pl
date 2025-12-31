#!/usr/bin/env perl

=head1 NAME

convert_RNAmmer_to_gff3.pl - convert bad gff2 format RNAmmer output to good GFF3

=head1 SYNOPSIS

USAGE: template.pl 
            --input=/path/to/some_file.out 
            
=head1 OPTIONS

B<--input,-i>
The raw output from RNAmmer:
##gff-version2
##source-version RNAmmer-1.2
##date 2014-03-26
##Type DNA
# seqname           source                      feature     start      end   score   +/-  frame  attribute
# ---------------------------------------------------------------------------------------------------------
tp.assembly.567497685.1 RNAmmer-1.2     rRNA    308022  312190  3003.7  -       .       28s_rRNA
tp.assembly.567468735.1 RNAmmer-1.2     rRNA    916705  922401  2985.9  -       .       28s_rRNA
tp.assembly.567492885.1 RNAmmer-1.2     rRNA    1034463 1034568 5.8     -       .       8s_rRNA
tp.assembly.567478335.1 RNAmmer-1.2     rRNA    762933  763050  69.8    -       .       8s_rRNA
tp.assembly.567478335.1 RNAmmer-1.2     rRNA    769881  769998  69.8    -       .       8s_rRNA
tp.assembly.567478335.1 RNAmmer-1.2     rRNA    775140  775257  69.8    -       .       8s_rRNA
tp.assembly.567492885.1 RNAmmer-1.2     rRNA    1036777 1036886 5.7     -       .       8s_rRNA
tp.assembly.567497685.1 RNAmmer-1.2     rRNA    312768  314511  1441.7  -       .       18s_rRNA
tp.assembly.567468735.1 RNAmmer-1.2     rRNA    922979  924722  1441.7  -       .       18s_rRNA
# ---------------------------------------------------------------------------------------------------------
B<--log,-l> 
    Log file
B<--help,-h>
    This help message
    
=head1  DESCRIPTION

File converter

=head1  INPUT

Input above.

=head1  OUTPUT

GFF3 to STDOUT

=head1  CONTACT

    huangls
    huangls@biomics.com.cn
    
=cut


use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage qw(pod2usage);


my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'log|l=s',
                          'ignoreIntrons|g',
                          'help|h') || pod2usage(1);

## display documentation


if( $options{'help'}){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
    #pod2usage();
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## set ignoreIntrons flag
my $ignoreIntrons = '';
if (defined $options{ignoreIntrons}) {
    $ignoreIntrons = 1;
}


## open the input file
my $ifh;
open($ifh, "<$options{input}") || die "can't open input file: $!";

# all output needs the gff header
print "##gff-version 3\n";

## globals
my $tRNAGenNum=1;
my $tRNAExonNum=1;

## parse the file
foreach my $line (<$ifh>){
	my @cols = split /[\t]/, $line;
	chomp @cols;
	my $contig = $cols[0];

    if ($contig =~ /^(.+?)\s+$/) {
        $contig = $1;
    }

    ## skip the header lines
    next if $contig eq 'Sequence' || $contig eq 'Name' || $contig eq '--------';
    
	my $start = $cols[2] + 0;
	my $stop = $cols[3] + 0;
	my $trna = $cols[4];
	my $score = $cols[8];
    my $pseudo = "";
	my $anticodon = $cols[5];
    my $intron_start = $cols[6]-1; ## offset intron boundaries
	my $intron_end = $cols[7]+1;

    $pseudo = ";pseudogene=unknown" if $cols[9] =~ m/pseudo/;
    # Suppressor (sup) tRNAs are not recognized by NCBI's table2asn
    $trna = "Xxx" if $trna =~ m/Sup/ or $trna =~ m/Undet/;
	
	## Check if tRNAScan-SE predicted a possible intron (value in column 6 > 0)
	if ($cols[6] > 0 and not $ignoreIntrons) {
   		if ($start > $stop){
			my $intron_end = $cols[6]+1; ## offset intron boundaries
			my $intron_start = $cols[7]-1;
            ($start, $stop) = ($stop, $start);
    	}
        print "$contig\ttRNAScan-SE\tgene\t$start\t$stop\t$score\t-\t.\tID=tRNA-$trna$tRNAGenNum\_gene$pseudo\n";
		print "$contig\ttRNAScan-SE\ttRNA\t$start\t$intron_start\t$score\t-\t.\tID=tRNA-$trna$tRNAExonNum\_tRNA;Parent=tRNA-$trna$tRNAGenNum\_gene;product=tRNA-$trna;anticodon=$anticodon;\n";
        $tRNAExonNum++;
		print "$contig\ttRNAScan-SE\ttRNA\t$intron_end\t$stop\t$score\t-\t.\tID=tRNA-$trna$tRNAExonNum\_tRNA;Parent=tRNA-$trna$tRNAGenNum\_gene;product=tRNA-$trna;anticodon=$anticodon;\n";
		$tRNAGenNum++;
        $tRNAExonNum++;
	} else{
		($start, $stop) = ($stop, $start) if ($start > $stop);
        print "$contig\ttRNAScan-SE\tgene\t$start\t$stop\t$score\t-\t.\tID=tRNA-$trna$tRNAGenNum\_gene$pseudo\n";
		print "$contig\ttRNAScan-SE\ttRNA\t$start\t$stop\t$score\t-\t.\tID=tRNA-$trna$tRNAExonNum\_tRNA;product=tRNA-$trna;anticodon=$anticodon;Parent=tRNA-$trna$tRNAGenNum\_gene\n";
		$tRNAGenNum++;
        $tRNAExonNum++;
	}
}

exit(0);


sub _log {
    my $msg = shift;
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;
    ## make sure required arguments were passed
    my @required = qw( input );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    ## handle some defaults
    $options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}
