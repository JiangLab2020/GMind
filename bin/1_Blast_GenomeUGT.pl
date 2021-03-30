use strict;
use warnings;
use Bio::Perl;

my $FamilyPath = $ARGV[0];
$FamilyPath=~s/\/$//;

my $CPU_Num = 20;
if ( $ARGV[1] ) {
	$CPU_Num = $ARGV[1];
}
my $Genome_dir = 0;
if ( $ARGV[2] ) {
	$Genome_dir = $ARGV[2];
}

my $Query_DB2 = "PUGTDB60_Cut.fasta";

my $GMind=$ENV{'GMind'};
my $GMind_lib=$ENV{'GMind_lib'};

opendir(DIR, $FamilyPath);
my @dots = readdir(DIR); 
closedir DIR; 

chdir($FamilyPath);

for (my $j=2;$j<@dots;$j++) {
	if ( $dots[$j] =~ /(.*)\.genome$/ or $dots[$j] =~ /(.*)\.fna$/) {
		my $Name = $1;
		if ( !(-e "$Name\_$Name\.fasta") ) {

			print $dots[$j],"\n";

			my @args = stat ($dots[$j]);
			my $size = $args[7]/1e+6;
			`perl $GMind/1_BlastSpeGene.pl $dots[$j] $GMind_lib/$Query_DB2 $Name $CPU_Num $Genome_dir`;
			
			unlink("$Name\.txt");
			unlink("$Name\.sort");
			unlink("$Name\.sort.num");
			unlink("$Name\.sort.combine");
			unlink("$Name\.sort.combine.short");			
		}
	}
}

