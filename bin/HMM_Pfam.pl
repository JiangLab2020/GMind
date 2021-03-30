use strict;
use warnings;
use Bio::SeqIO;

my $InputPep = $ARGV[0];
my $Pfam     = $ARGV[1]; 

my $CPU		 = 20;
if ( $ARGV[2] ) {
	$CPU	 = $ARGV[2]
}

my $GMind=$ENV{'GMind'};
my @b=split(/[\/\.]/,$InputPep); my $Target=$b[$#b-1];
######################################## 1.HMM $Pfam #################################
print "################# 1.HMM $Pfam #################","\n";
`perl $GMind/Hmmer_MultiThreads.pl $InputPep $CPU $Pfam`;
my %AList=();
open INPUT,"$Target\_Hmm.txt";
while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ($_=~ /$Pfam/) {
		$AList{$a[0]}=1;
	}
}
close INPUT;
print "1.HMM $Pfam(Fetch PEP)","\n";  
my $Target1="$Target\_1.fasta";
open OUTPUT,">$Target1";
my $seqdata = $InputPep;
my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) { 
	my $id=$seq->display_id;
	if ($AList{$id}) {
		print OUTPUT ">",$seq->display_id," ",$seq->description,"\n";
		print OUTPUT $seq->seq(),"\n";
	}
}
close OUTPUT;
######################################## 2.HMM ALL #################################
print "################# 2.HMM ALL #################","\n";
`perl $GMind/Hmmer_MultiThreads.pl $Target1 $CPU`;
my %BList=();
open INPUT,"$Target\_1_Hmm.txt";
while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ($_=~ /$Pfam/) {
		$BList{$a[0]}=join(" ",@a[1..$#a]);
	}
}
close INPUT;
print "2.HMM ALL(Fetch PEP)","\n";
open OUTPUT,">$Target\_$Pfam\.fasta";
my $seqdata1 = $Target1;
my $seqio1 = new Bio::SeqIO(-file => $seqdata1, -format => 'fasta');
while ( my $seq = $seqio1->next_seq ) { 
	my $id=$seq->display_id;
	if ( $BList{$id} ) {
		print OUTPUT ">",$seq->display_id," ",$seq->description," ",$BList{$id},"\n";
		print OUTPUT $seq->seq(),"\n";
	}
}
close OUTPUT;

unlink("$Target1");
unlink("$Target\_1_Hmm.txt");
unlink("$Target\_Hmm.txt");

