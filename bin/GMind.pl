use strict;
use warnings;

my $GenomeFile = $ARGV[0];
my $CPU		   = $ARGV[1];

my $Pfam	   = "PF00201";
if ( $ARGV[2] ) {
	$Pfam = $ARGV[2];
}
my $Query_DB   = "PUGTDB60_Cut.fasta";
if ( $ARGV[3] ) {
	$Query_DB = $ARGV[3];
}

my $GMind=$ENV{'GMind'};
my $GMind_lib=$ENV{'GMind_lib'};

$GenomeFile=~s/\/$//;
my $Path = `pwd`;  $Path=~s/\n$//;

print '##############################    ',"1: ","Blast ProteinDB to Genome\n";

my $TargetDir  = $GenomeFile;							
if ( !(-e $TargetDir) ) {
	`mkdir $TargetDir`;
}
my $ResultDir  = "$GenomeFile/Resultdir";
if ( !(-e $ResultDir) ) {
	`mkdir $ResultDir`;
}
else {
	`rm -rf $ResultDir/*`;
}

my $Genome_dir="$GenomeFile/GenomeAnnotation";
if (! (-e $Genome_dir)) {
	`mkdir $Genome_dir`;
}
else {
	`rm -rf $Genome_dir/*`;
}

my $TargetDirBlast  = "$Genome_dir/_Blast";                   
if ( !(-e $TargetDirBlast) ) {
	`mkdir $TargetDirBlast`;
}
my $Blast_MapDIR  = "$Genome_dir/Genome_Map"; 
if ( !(-e $Blast_MapDIR) ) {
	`mkdir $Blast_MapDIR`;
}

opendir(DIR, $TargetDir);
my @dots = readdir(DIR); 
closedir(DIR);
chdir($TargetDir); 

open OUTPUT,">rename_genome.txt";
print OUTPUT "Old_name","\t","New_name","\n";
for (my $j=2;$j<@dots ;$j++) {
	if ($dots[$j]=~ /(.*)\.fasta$/ or $dots[$j]=~ /(.*)\.fna$/) {
		my $species=$1;
		my $Abbr        = "";           
		if ( $species =~ /(^[A-Z]?).*\_([a-z][a-z][a-z]?).*/ ) {
			$Abbr=$1.$2;
		}
		else {
			$Abbr=$species;
		}
		`cp $dots[$j] $Blast_MapDIR/$Abbr\.genome`;
		`cp $dots[$j] $Genome_dir`;
		print OUTPUT $dots[$j],"\t",$Abbr,".genome","\n";
	}
}
close OUTPUT;
`perl $GMind/1_Blast_GenomeUGT.pl $Blast_MapDIR $CPU $Genome_dir`;

print '##############################    ',"2: ","Three methods Annotate homo-Genomes\n";
chdir("$Path/Genomedir");

`perl $GMind/Genome_Annoation.pl $Genome_dir $CPU $Pfam $Query_DB`;

print '##############################    ',"3: ","Filtering good proteins in homo-Genomes\n";

`perl $GMind/GeneAnn_TrimRightSeq.pl $Genome_dir`;

`mv $Genome_dir/Trim_All.fasta $TargetDirBlast/`;

print '##############################    ',"4: ","Blast bad genome to ProteinDB, fetch protein region\n";

chdir $TargetDirBlast;

`perl $GMind/GeneAnn_BlastGene.pl Trim_All.fasta $GMind_lib/$Query_DB 0 $CPU`;

my $Blast_Result  = "Blast_Result";                 
if ( !(-e $Blast_Result) ) {
	`mkdir $Blast_Result`;
}

`mv 0.* $Blast_Result/`;   `mv 0_Uni.txt $Blast_Result/`;  `mv Trim_All.fasta $Blast_Result/`; 

print '##############################    ',"5: ","Three methods Annotate Filtering Genomes\n";
chdir $Genome_dir;

`perl $GMind/Genome_Annoation.pl $TargetDirBlast $CPU $Pfam $Query_DB`;

print '##############################    ',"6: ","Combining the annotations of homo-Genomes and Filtering Genomes\n";

`cat All_uni_DB_adj_div.pep _Blast/All_uni_DB_adj_div.pep >All_uni_DB_adj_div_Com.pep`;
`cat All_uni_DB_adj_div.txt _Blast/All_uni_DB_adj_div.txt >All_uni_DB_adj_div_Com.txt`;


`perl $GMind/Uni_GeneAnn3.pl All_uni_DB_adj_div_Com.txt All_uni_DB_adj_div_Com.pep 0.2`;
`cp All_uni_DB_adj_div_Com_Final.pep $ResultDir/GMind.pep`;
`cp All_uni_DB_adj_div_Com_Final.txt $ResultDir/GMind.txt`;

unlink("All_uni_DB_adj_div_good.pep");
unlink("All_uni_DB_adj_div_bed.pep");
