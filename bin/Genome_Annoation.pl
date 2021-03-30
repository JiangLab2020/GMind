use strict;
use warnings;
use Bio::SeqIO;

my $TargetDir  = $ARGV[0];
my $CPU		   = $ARGV[1];
my $Pfam	   = $ARGV[2];
my $Query_DB   = $ARGV[3];

$TargetDir=~s/\/$//;

opendir DIR,$TargetDir;
my @dots = readdir(DIR);
close(DIR);

my $GMind=$ENV{'GMind'};
my $GMind_lib=$ENV{'GMind_lib'};
`cp $GMind_lib/$Query_DB $TargetDir/Query_DB.pep`;


chdir $TargetDir;

my $Auga_dir="AugAth_dir";
if (! (-e $Auga_dir)) {
	`mkdir $Auga_dir`;
}
my $Augr_dir="AugPUGT_dir";
if (! (-e $Augr_dir)) {
	`mkdir $Augr_dir`;
}
my $Glm_dir="Glimmhmm_dir";
if (! (-e $Glm_dir)) {
	`mkdir $Glm_dir`;
}
my $Merge_dir="Merge_dir";
if (! (-e $Merge_dir)) {
	`mkdir $Merge_dir`;
}

open OUTPUT,">Genome_size.error";

for (my $j=2;$j<@dots ;$j++) {
	my @b=split(/[\/\.]/,$dots[$j]); my $species=$b[$#b-1];
	my @args=stat($dots[$j]);
	my $size=$args[7];
	my @Files=();
	if ( $dots[$j] =~ /\.fasta/ and $size > 0 and $dots[$j] !~ /Trim_All.fasta/ ) {
		if ( $dots[$j] !~ /\_$Pfam\.fasta/ ) {
			if (! -e "$Auga_dir/$species\_augustus.pep") {
				`perl $GMind/Genome_augustus.pl $dots[$j] $CPU`;
				`perl $GMind/HMM_Pfam.pl $species\_augustus.pep $Pfam $CPU`;
				`mv $species\_augustus.pep $Auga_dir`;
				`mv $species\_augustus.cds $Auga_dir`;
				`mv $species\_augustus.gff $Auga_dir`;
			}
			if (! -e "$Augr_dir/$species\_augPU.pep") {
				`perl $GMind/Genome_augPUGT.pl $dots[$j] $CPU`;
				`perl $GMind/HMM_Pfam.pl $species\_augPU.pep $Pfam $CPU`;
				`mv $species\_augPU.pep $Augr_dir`;
				`mv $species\_augPU.cds $Augr_dir`;
				`mv $species\_augPU.gff $Augr_dir`;
			}
			if (! -e "$Glm_dir/$species\_glimmerHMM.pep") {
				`perl $GMind/Genome_glimmhmm.pl $dots[$j]`;
				`perl $GMind/HMM_Pfam.pl $species\_glimmerHMM.pep $Pfam $CPU`;
				`mv $species\_glimmerHMM.cds $Glm_dir`;
				`mv $species\_glimmerHMM.pep $Glm_dir`;
				`mv $species\_glimmerHMM.gff $Glm_dir`;
			}
		}

		if ( !(-e "$Merge_dir/$species\_uni_DB.pep") ) {
			if ( !(-e "$species\_augustus\_$Pfam\.fasta") and -e "$Auga_dir/$species\_augustus\_$Pfam\.fasta" ) {
				`cp $Auga_dir/$species\_augustus\_$Pfam\.fasta .`;
				`cp $Augr_dir/$species\_augPU\_$Pfam\.fasta .`;
				`cp $Glm_dir/$species\_glimmerHMM\_$Pfam\.fasta .`;
			}

			if (-e "$species\_augustus\_$Pfam\.fasta" ) {
				my @arg=stat("$species\_augustus\_$Pfam\.fasta");
				my $size1=$arg[7];
				if ($size1 > 0) {
					push @Files,"$species\_augustus\_$Pfam\.fasta";
				}
			}
			if (-e "$species\_augPU\_$Pfam\.fasta" ) {
				my @arg=stat("$species\_augPU\_$Pfam\.fasta");
				my $size1=$arg[7];
				if ($size1 > 0) {
					push @Files,"$species\_augPU\_$Pfam\.fasta";
				}
			}
			if (-e "$species\_glimmerHMM\_$Pfam\.fasta" ) {
				my @arg=stat("$species\_glimmerHMM\_$Pfam\.fasta");
				my $size1=$arg[7];
				if ($size1 > 0) {
					push @Files,"$species\_glimmerHMM\_$Pfam\.fasta";
				}
			}

			if ( @Files ) {
				my $file=join(" ",@Files);
				print "Merge: ",$dots[$j],"\t",$file,"\n";
				`perl $GMind/Uni_GeneAnn1.pl $file`;
				`mv GeneLoc_Sort_uni_DB.pep $Merge_dir/$species\_uni_DB.pep`;
				`mv GeneLoc_Sort_class.txt $Merge_dir/$species\_class.txt`;
				`mv GeneLoc_div.txt $Merge_dir/$species\_div.txt`;
				`mv $species\_augustus\_$Pfam\.fasta $Auga_dir`;
				`mv $species\_augPU\_$Pfam\.fasta $Augr_dir`;
				`mv $species\_glimmerHMM\_$Pfam\.fasta $Glm_dir`;
			}
			if ( -e "$species\_augustus\_$Pfam\.fasta" and !(-e "$Auga_dir/$species\_augustus\_$Pfam\.fasta" )) {
				`mv $species\_augustus\_$Pfam\.fasta $Auga_dir/`;
				`mv $species\_augPU\_$Pfam\.fasta $Augr_dir/`;
				`mv $species\_glimmerHMM\_$Pfam\.fasta $Glm_dir/`;
			}
		}

	}
	if ($dots[$j]=~ /\.fasta/ and $size == 0) {
		print OUTPUT $dots[$j],"\n";
	}
}

close OUTPUT;

unlink("Tem_Tem.fas");
unlink("Tem_Tem.pep");
unlink("Tem_Tem.txt");

if ( !(-e "All_uni_DB_adj_div.pep") ) {
	`cat AugAth_dir/*.fasta AugPUGT_dir/*.fasta Glimmhmm_dir/*.fasta >All_Methods.pep`;
	`perl $GMind/Fasta_div.pl All_Methods.pep Query_DB.pep $CPU`;
	`cat Merge_dir/*_uni_DB.pep >All_uni_DB.pep`;

	`perl $GMind/Adjust_div.pl`;
	`perl $GMind/Uni_GeneAnn2.pl All_uni_DB_adj_div.txt All_uni_DB_adj.pep 0.2`;
}
