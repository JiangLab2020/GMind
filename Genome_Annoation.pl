use strict;
use warnings;
use Bio::SeqIO;

my $TargetDir  = $ARGV[0];                   ####  基因组所在文件夹
my $CPU		   = $ARGV[1];
my $Pfam	   = $ARGV[2];
my $Query_DB   = $ARGV[3];                   ####  注释依赖的蛋白库文件

$TargetDir=~s/\/$//;

opendir DIR,$TargetDir;
my @dots = readdir(DIR);
close(DIR);

`cp /public/home/chengj/3_LiuYuQian/Genome_Annotation/genome/Genome_augustus.pl $TargetDir/`;
`cp /public/home/chengj/3_LiuYuQian/Genome_Annotation/genome/Genome_augrice.pl $TargetDir/`;
`cp /public/home/chengj/3_LiuYuQian/Genome_Annotation/genome/Genome_glimmhmm.pl $TargetDir/`;
`cp /public/home/chengj/3_LiuYuQian/Genome_Annotation/genome/LocRegionUni.pl $TargetDir/`;

`cp $Query_DB $TargetDir/Query_DB.pep`;

`cp /public/home/chengj/Database/Genomes_Map/Viridiplantae/Uni_GeneAnn1.pl $TargetDir/`;
`cp /public/home/chengj/Database/Genomes_Map/Viridiplantae/HMM_Pfam.pl $TargetDir/`;
`cp /public/home/chengj/Database/Genomes_Map/Viridiplantae/Fasta_div.pl $TargetDir/`;
`cp /public/home/chengj/Database/Genomes_Map/Viridiplantae/Adjust_div.pl $TargetDir/`;

chdir $TargetDir;
my $Auga_dir="AugAth_dir";
if (! (-e $Auga_dir)) {
	`mkdir $Auga_dir`;
}
my $Augr_dir="AugRice_dir";
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
	my @b=split(/[\/\.]/,$dots[$j]); my $species=$b[$#b-1]; my $len_species=length($species);
	my @args=stat($dots[$j]);
	my $size=$args[7];
	my @Files=();
	if ($dots[$j]=~ /\.fasta/ and $size>0) {
		if ($len_species < 10) {
			if (! -e "$Auga_dir/$species\_augustus.pep") {
				`perl Genome_augustus.pl $dots[$j] $CPU`;
				`perl HMM_Pfam.pl $species\_augustus.pep $Pfam $CPU`;
				`mv $species\_augustus.pep $Auga_dir`;
				`mv $species\_augustus.cds $Auga_dir`;
				`mv $species\_augustus.gff $Auga_dir`;
			}
			if (! -e "$Augr_dir/$species\_augrice.pep") {
				`perl Genome_augrice.pl $dots[$j] $CPU`;
				`perl HMM_Pfam.pl $species\_augrice.pep $Pfam $CPU`;
				`mv $species\_augrice.pep $Augr_dir`;
				`mv $species\_augrice.cds $Augr_dir`;
				`mv $species\_augrice.gff $Augr_dir`;
			}
			if (! -e "$Glm_dir/$species\_glimmerHMM.pep") {
				`perl Genome_glimmhmm.pl $dots[$j]`;
				`perl HMM_Pfam.pl $species\_glimmerHMM.pep $Pfam $CPU`;
				`mv $species\_glimmerHMM.cds $Glm_dir`;
				`mv $species\_glimmerHMM.pep $Glm_dir`;
				`mv $species\_glimmerHMM.gff $Glm_dir`;
			}
		}

		if (-e "$species\_augustus\_$Pfam\.fasta" ) {
			my @arg=stat("$species\_augustus\_$Pfam\.fasta");
			my $size1=$arg[7];
			if ($size1 > 0) {
				push @Files,"$species\_augustus\_$Pfam\.fasta";
			}
		}
		if (-e "$species\_augrice\_$Pfam\.fasta" ) {
			my @arg=stat("$species\_augrice\_$Pfam\.fasta");
			my $size1=$arg[7];
			if ($size1 > 0) {
				push @Files,"$species\_augrice\_$Pfam\.fasta";
			}
		}
		if (-e "$species\_glimmerHMM\_$Pfam\.fasta" ) {
			my @arg=stat("$species\_glimmerHMM\_$Pfam\.fasta");
			my $size1=$arg[7];
			if ($size1 > 0) {
				push @Files,"$species\_glimmerHMM\_$Pfam\.fasta";
			}
		}
		if (! -e "$Merge_dir/$species\_uni_DB.pep" and @Files) {
			my $file=join(" ",@Files);
			print "Merge: ",$dots[$j],"\t",$file,"\n";
			`perl Uni_GeneAnn1.pl $file`;
			`mv GeneLoc_Sort_uni_DB.pep $Merge_dir/$species\_uni_DB.pep`;
			`mv GeneLoc_Sort_class.txt $Merge_dir/$species\_class.txt`;
			`mv GeneLoc_Sort_uni.txt $Merge_dir/$species\_uni.txt`;
			`mv GeneLoc_div.txt $Merge_dir/$species\_div.txt`;
			`mv $species\_augustus\_$Pfam\.fasta $Auga_dir`;
			`mv $species\_augrice\_$Pfam\.fasta $Augr_dir`;
			`mv $species\_glimmerHMM\_$Pfam\.fasta $Glm_dir`;
		}
		if (! @Files) {
			print "Merge: ",$dots[$j],"\t","Nothing to merge!","\n";
		}
	}
	if ($dots[$j]=~ /\.fasta/ and $size==0) {
		print OUTPUT $dots[$j],"\n";
	}
}

close OUTPUT;

unlink("Tem_Tem.fas");
unlink("Tem_Tem.pep");
unlink("Tem_Tem.txt");

`cat AugAth_dir/*.fasta AugRice_dir/*.fasta Glimmhmm_dir/*.fasta >All_Methods.pep`;
`perl Fasta_div.pl All_Methods.pep Query_DB.pep $CPU`;
`cat Merge_dir/*_uni_DB.pep >All_uni_DB.pep`;

`perl Adjust_div.pl`;
