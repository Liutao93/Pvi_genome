#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://study.163.com/course/introduction/1006448023.htm?share=1&shareId=1030291076
#perl高级：https://study.163.com/course/introduction/1004833023.htm?share=1&shareId=1030291076


use strict;
use Bio::SeqIO;
die "perl $0 <gff> <genome.fa>  <cds.hint>" unless(@ARGV==3);

my $gff_file = shift;
my $genome_file = shift;
#my $op1 = shift; #fasta序列
my $op2 = shift; #cds位置信息的文件
my @scaffold_name;
open (G,"<$genome_file") || die ("$! $genome_file");
open (F,"<$gff_file") || die ("$! $gff_file");
#open (O1,">$op1") || die ("!/n");
open (O2,">$op2") || die ("!/n");
#my $fa = Bio::SeqIO->new(-file=>$genome_file, -format=>'fasta');
my %seq;
my %hash;
while(my$line=<G>){
	chomp $line;
    if($line=~/^>/){
	    my ($id )= ($line=~/^>(\S+)/);
    	push(@scaffold_name,$id);
	# print O1 "$line\n";
    }else{
	#print O1 "$line\n";
    }
}
close G;

my %strand;

while (my $eve = <F>){
    chomp ($eve);
    my @infor = split/\s+/,$eve;
    my $scaffold_name = $infor[0];
    my $type = $infor[2];
    
    my $strand = $infor[6];
    if ($type eq "exon"){
      my( $gene_name) = $infor[8]=~/Parent=([^;]+)/;
      	 $strand{$gene_name}=$strand;
      	if ($strand eq "-"){
      	    
      	    push @{$hash{$scaffold_name}{$gene_name} }, [$scaffold_name,$infor[4],$infor[3]];
      	}elsif($strand eq "+"){
      	   push @{$hash{$scaffold_name}{$gene_name} }, [$scaffold_name,$infor[3],$infor[4] ];
      	}
    }
}

foreach my $scaffold_name(@scaffold_name){
    if (exists $hash{$scaffold_name}){
    	foreach my $gene_name(sort {$hash{$scaffold_name}{$a}->[0]->[1] <=> $hash{$scaffold_name}{$b}->[0]->[1] } keys %{$hash{$scaffold_name}}){
    	  if($strand{$gene_name} eq "-"){
    	    for my $cds (sort {$b->[1] <=> $a->[1] } @{$hash{$scaffold_name}{$gene_name}}){
    	      print O2 join("\t",@{$cds})."\n";
    	    }
    	     
    	  }elsif($strand{$gene_name} eq "+"){
    	    for my $cds (sort {$a->[1] <=> $b->[1]} @{$hash{$scaffold_name}{$gene_name}}){
    	      print O2 join("\t",@{$cds})."\n";
    	    }
    	  }else{
    	    die "gene : $gene_name don't have strand info in gff file \n";
    	  }
    	  
    	 
    	 print O2 "\n";
    	}
    }
}
close F;
close O1;
close O2;
