use Bio::Perl;
use strict;
use warnings;

my $input  = $ARGV[0];          
my $output = $ARGV[1];          
my $Cut    = $ARGV[2];          

my $seqdata = "$input";

my $seqio = new Bio::SeqIO(-file   => $seqdata,-format => 'fasta');
my @Ami=(); my @Gap=(); my @b=();
while ( my $seq = $seqio->next_seq ) {
	my $FaSeq=$seq->seq();
	$FaSeq=~s/[A-Z]/\+/g;
	@b=split(//,$FaSeq);
	for (my $i=0;$i<@b ;$i++) {
		if ($b[$i] eq "+") {
			$Ami[$i]+=1;
		}
		else {
			$Gap[$i]+=1;
		}
	}
}

my @Com=();
for (my $i=0;$i<@b ;$i++) {
	if (!$Ami[$i]) {
		$Ami[$i]=0;
	}
	if (!$Gap[$i]) {
		$Gap[$i]=0;
	}
	if ( $Ami[$i] > $Gap[$i] ) {
		$Com[$i]="+";
	}
	else {
		$Com[$i]="-";
	}
}

open OUTPUT,">$output";

$seqio = new Bio::SeqIO(-file => $seqdata,-format => 'fasta');
my $Dif=0;
while ( my $seq = $seqio->next_seq ) {  
	my $FaSeq=$seq->seq();
	$FaSeq=~s/[A-Z]/\+/g;
	my @b=split(//,$FaSeq);
	for (my $i=0;$i<@b ;$i++) {
		if ($b[$i] ne $Com[$i]) {
			$Dif++;
		}
	}
	$Dif=Float($Dif/length($FaSeq));
	if ( $Dif >= $Cut ) {
		print OUTPUT $seq->display_id,"\t",$Dif,"\t",$input,"\n";
	}
}

close(OUTPUT);

############ 
sub Float{
	my $str=shift;
	if ($str =~ /\./) {
		$str = sprintf "%0.3f",$str;
	}
	return $str;
}
