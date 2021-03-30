use strict;
use warnings;
use Bio::Perl;

my $InputFA = $ARGV[0];

my $GMind=$ENV{'GMind'};

my @b=split(/[\/\.]/,$InputFA); my $Target=$b[$#b-1];
my @spe=split(/\_/,$Target);
my $Abbr        = $spe[0];
my %AList=();

my $TrainPath = "/public/home/chengj/Programs/GlimmerHMM/trained_dir/arabidopsis";

my %BList=();
if ( !(-e "$Target\_glimmerHMM.gff") ) {
	`glimmhmm.pl glimmerhmm_linux_x86_64 $InputFA $TrainPath -g >$Target\_glimmerHMM.gff`;
	my $GeneNum=1;
	my $Str="10000000";
	my $NameID=();
	open INPUT,"$Target\_glimmerHMM.gff";
	open OUTPUT,">$Target\_glimmerHMM1.gff";
	while (<INPUT>) {
		chomp($_);
		my @a=split(/\t/,$_);
		my @b=split(/\;/,$a[$#a]);
		if ($a[2] and $a[2] eq "mRNA") {
			$NameID=$Abbr."\_glm_".substr($Str,0,length($Str)-length($GeneNum)).$GeneNum; $GeneNum++;
			$BList{$NameID}=$a[0]." ".$a[3]."_".$a[4]." ".$a[6];
			$b[0]="ID=".$NameID;
			$a[$#a]=join("\;",@b);
			print OUTPUT join("\t",@a),"\n";
		}
		if ($a[2] and $a[2] eq "CDS") {
			$b[1]="Parent=".$NameID;
			$b[2]="Name=".$NameID;
			$a[$#a]=join("\;",@b);
			print OUTPUT join("\t",@a),"\n";
		}
	}
	close INPUT;
	close OUTPUT;
	`mv $Target\_glimmerHMM1.gff $Target\_glimmerHMM.gff`;
}

my %CList = ();
my $seqdata = $InputFA;
my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {   
	$CList{$seq->display_id}=$seq->seq();
}

my $Tem_Gff1 = "CDS1_".$$."\.gff";
my $Tem_Gff2 = "CDS2_".$$."\.gff";

open INPUT,"$Target\_glimmerHMM.gff";
open OUTPUT,">$Tem_Gff1";
while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ( $a[2] and $a[2] eq "CDS" ) {
		my $ID="";
		if ( $a[8] =~ /gene_id \"(.*?)\"/ ) {
			$ID=$1;
		}
		elsif ( $a[8] =~ /Parent\=(.*?)\./ ) {
			$ID=$1;
		}
		elsif ( $a[8] =~ /Parent\=(.*?)\;/ ) {
			$ID=$1;
		}
		$a[8]=$ID;
		if ( $a[8] ) {
			print OUTPUT join("\t",@a),"\n";
		}
		else {
			print "No gene ID.....\n";
		}
	}
}
close(INPUT);

`perl $GMind/Tab_Rank.pl $Tem_Gff1 $Tem_Gff2 1 9 4`;

my @Gene  = ();
my %GList = ();
my %SList = ();

open INPUT,$Tem_Gff2;
while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ( $a[2] and $a[2] eq "CDS" ) {
		my $ID=$a[8];
		if ( !$GList{$ID} ) {
			$GList{$ID}=FetchLocSeq($CList{$a[0]},$a[3],$a[4],"+");
			$SList{$ID}=$a[6];
			push @Gene,$ID;
		}
		else {
			$GList{$ID}.=FetchLocSeq($CList{$a[0]},$a[3],$a[4],"+");
		}
	}
}
close(INPUT);

open OUTPUT,">$Target\_glimmerHMM.cds";

foreach my $ID (@Gene) {
	my $GeneSeq=$GList{$ID};
	if ( $SList{$ID} =~ /\-/ ) {
		my $seqobj = Bio::Seq->new( -display_id => 'my_id',-seq => $GeneSeq);
		$GeneSeq=$seqobj->revcom->seq();
	}
	if ( $AList{$Abbr} ) {
		print OUTPUT ">",$ID," ",$BList{$ID}," ","[",$AList{$Abbr},"]","\n",$GeneSeq,"\n";
	}
	else {
		print OUTPUT ">",$ID," ",$BList{$ID}," ","[",$Target,"]","\n",$GeneSeq,"\n";
	}
}

close(OUTPUT);

open OUTPUT,">$Target\_glimmerHMM_Adj.cds";

my $seqio2 = new Bio::SeqIO(-file   => "$Target\_glimmerHMM.cds",-format => 'fasta');
while ( my $seq = $seqio2->next_seq ) {
	my $Len=length($seq->seq());  my $MySeq = $seq->seq();
	my @SeqDNA=();
	$SeqDNA[0]=$seq->seq();
	$SeqDNA[1]=substr($SeqDNA[0],1,$Len-1);
	$SeqDNA[2]=substr($SeqDNA[0],2,$Len-2);
	
	if ( $Len % 3 ne 0 ) {
		my $k = 0;
		for (my $i=0;$i<@SeqDNA ;$i++) {
			my $seqobj1 = Bio::Seq->new( -display_id => 'my_id',-seq => $SeqDNA[$i]);
			my $protein=$seqobj1->translate()->seq();
			$protein =~ s/\*$//;
			if ( $protein =~ /\*/ ) {				
				$k++;
			}
			else {
				$MySeq = $SeqDNA[$i];
			}
		}
		print $seq->display_id,"\t",$k,"\n";
	}
	print OUTPUT ">",$seq->display_id," ",$seq->description,"\n";
	print OUTPUT $MySeq,"\n";
}

close(OUTPUT);

open OUTPUT,">$Target\_glimmerHMM.pep";
################################################# cds2pep #################################################
my $seqio1 = new Bio::SeqIO(-file   => "$Target\_glimmerHMM_Adj.cds",-format => 'fasta');
while ( my $seq = $seqio1->next_seq ) {
	print OUTPUT ">",$seq->display_id," ",$seq->description,"\n";
	print OUTPUT $seq->translate()->seq(),"\n";
}
close(OUTPUT);

unlink("$Tem_Gff1");
unlink("$Tem_Gff2");
`rm "$Target\_glimmerHMM.cds"`;
`mv $Target\_glimmerHMM_Adj.cds $Target\_glimmerHMM.cds`;

#####################################################################################################################################################

sub FetchLocSeq{

	my ($MySeq,$a1,$a2,$Strand) = @_;
	my $str1="";
	if ( $a1 > 0 ) {
		$str1=substr($MySeq,$a1-1,$a2-$a1+1);
		if ( $Strand eq "-" ) {
			my $seqobj = Bio::Seq->new( -display_id => 'my_id',-seq => $str1);
			$str1=$seqobj->revcom->seq();
		}
	}
	return uc($str1);
}
