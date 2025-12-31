die "perl $0 <fa>  <outdir> [cut length 50000000]" if @ARGV==0;
use Bio::SeqIO;
use Bio::Seq;

mkdir "$ARGV[1]" unless(-d $ARGV[1]);



my $in  = Bio::SeqIO->new(-file => "$ARGV[0]" ,
                               -format => 'Fasta');

my$i=1;
my$first=1;
my$out="";
my$s=0;
my $cut_len=50000000;
if (@ARGV==3) {$cut_len=$ARGV[2];}
while ( my $seq = $in->next_seq() ) {
	my($id,$desc,$len)=($seq->id,$seq->desc,$seq->length);
	
        if($len>$cut_len ){
                $out = Bio::SeqIO->new(-file => ">$ARGV[1]/seq$i.fa" ,
                               -format => 'Fasta');

                $i++;
                $s=0;
        }elsif($s>$cut_len ){
                $out = Bio::SeqIO->new(-file => ">$ARGV[1]/seq$i.fa" ,
                               -format => 'Fasta');

                $i++;
                $s=0;
        }elsif($first and $len<=$cut_len and $i==1){
                $out = Bio::SeqIO->new(-file => ">$ARGV[1]/seq$i.fa" ,
                               -format => 'Fasta');

                $i++;
                $s=0;
                $first=0;

        }
        $s+=$len;
        $out->write_seq($seq);

}

$out->close();
$in->close();
