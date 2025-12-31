#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://study.163.com/course/introduction/1006448023.htm?share=1&shareId=1030291076
#perl高级：https://study.163.com/course/introduction/1004833023.htm?share=1&shareId=1030291076

use Getopt::Long;
use strict;
use Cwd qw(abs_path getcwd);

my %opts;

GetOptions (\%opts,"od=s","colline=s","gff=s","name=s"); 

my $od=$opts{od};
$od||=getcwd;
$od=abs_path($od);
unless(-d $od){	mkdir $od;}

my$name=$opts{name};
$name||="genome";
##############################





open (IN,"$opts{gff}") || die "open $opts{gff} failed\n";
my %gff;
my @info;
my $chr;
my $start;
my $end;
my $gene;
while(<IN>){
	chomp;
	next if /^#/;
	
	@info=split(/\t/,$_);
	
	#next unless($info[2]=~/gene/);
	#($gene)=($info[8]=~/ID=([^;]+)/);
	
	$gene=$info[1];
	$chr=$info[0];
	$start=$info[2];
	$end=$info[3];
	$gff{$gene}=$chr."\t".$start."\t".$end;
}

close(IN);



####################### list ##############


my %list;
my $pair;
my $Len;
my $Agene;
my $Bgene;



######### collinearity for genome block colline #####

open (IN,"$opts{colline}") || die "open $opts{colline} failed\n";
open (OUT,">$od/$name.blocklink.txt") || die "open $od/$name.blocklink.txt failed\n";
open (OUTA,">$od/$name.align.blocklink.txt") || die "open $od/$name.align.blocklink.txt failed\n";
my $n;
my $align;
my $colline;
my %block;
my $Agene1S;
my $AgeneNE;
my $Bgene1S;
my $BgeneNE;
my $Achr;
my $Bchr;
while(<IN>){
	chomp;
	if(/^#/){
		if(/Alignment/){
			$n=1;
			$_=~/Alignment ([^:]*)/;
			$align="Alignment".$1;
		}
		next;
	}
	
	$colline=$_;
	@info=split("\t",$colline);
	$Agene=$info[1];
	$Bgene=$info[2];
	
	if(exists $gff{$Agene} && exists $gff{$Bgene} ){
	
	    if($n ==1 ){
		
		    ($chr,$start,$end)=split(/\t/,$gff{$Agene});
	    	$Agene1S=$start;
		    $Achr=$chr;
		
		    ($chr,$start,$end)=split(/\t/,$gff{$Bgene});
		    $Bgene1S=$start;
		    $Bchr=$chr;
		
	    }else{
		
				
		    ($chr,$start,$end)=split(/\t/,$gff{$Agene});
		    $AgeneNE=$end;
		
		    ($chr,$start,$end)=split(/\t/,$gff{$Bgene});
		    $BgeneNE=$end;
		
	    }
	}
	$n=$n+1;
	$block{$align}=$Achr."\t".$Agene1S."\t".$AgeneNE."\t".$Bchr."\t".$Bgene1S."\t".$BgeneNE;	
		
	
}

close(IN);

my $block_info;

while(($align,$block_info)=each %block){
	print OUT $block_info."\n";
	print OUTA $align."\t".$block_info."\n";
}
close(OUT);
close(OUTA);






