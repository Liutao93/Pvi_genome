die "perl $0 <fa> <gff> <out>" unless @ARGV==3;
use Bio::SeqIO;
use Bio::Seq;
 
open IN,"$ARGV[1]" or die "$!; can't open file $ARGV[1]\n";
my%gff;
while(<IN>){
	chomp;
	my@a=split(/\t/);	
	if($a[2] eq "mRNA" || $a[2] eq "transcript"){
		my($mRNAid)=($a[8]=~ m/ID=([^;]*)/);
		my($geneid)=($a[8]=~ m/Parent=([^;]*)/);
		$gff{$mRNAid}=$geneid;
		}	
}
close IN;

$in  = Bio::SeqIO->new(-file => "$ARGV[0]" ,
                               -format => 'Fasta');
$out = Bio::SeqIO->new(-file => ">$ARGV[2]" ,
                               -format => 'Fasta');
my%fa=();
while ( my $seq = $in->next_seq() ) {
	my($id,$s,$desc,$len)=($seq->id,$seq->seq,$seq->desc,$seq->length);
	my $gene_id="";
	if(exists $gff{$id}){
		$gene_id=$gff{$id};
	}else{
		print STDERR "mRNA seq $id not in gff file, please check\n";
	next;

	}
	if(exists $fa{$gene_id}){
		if($fa{$gene_id}->length < $len){
			$fa{$gene_id}=$seq;
		}
	}else{
		$fa{$gene_id}=$seq;
	}
}

for my$id(keys %fa){
	
	my$seq_obj=Bio::Seq->new(-seq =>$fa{$id}->seq,
							-desc=>$id."=".$fa{$id}->length,
							-id=>$id
							 );
	$out->write_seq($seq_obj);
}
$in->close();
$out->close();
