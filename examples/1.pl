use strict;
use warnings;
use Bio::Perl;

my $k=0;
open OUTPUT,">Eryas_drummond.fna";
my $seqio=new Bio::SeqIO(-file=>"/public/home/chengj/3_LiuYuQian/Genome_Annotation/genome/PUGTA/bin/Dryas_drummondii.fna", -format=>'fasta');
while (my $seq =$seqio->next_seq) {
	if ($k>601 and $k<1101) {
		print OUTPUT ">",$seq->display_id," ",$seq->description,"\n",$seq->seq(),"\n";
	}
	$k++;
}
print $k,"\n";
close OUTPUT;