# Gene_Modeller
A novel plant UGTs annotated process for the mining of UGTs in plant genomes

# Overview

Gene_Modeller is a command line program (in Perl) based on Augustus, GlimmerHMM, 
which is used for sufficiently mining of UDP-glycosyltransferases from plant genomes.


# Installation

Running Gene_Modeller requires a Linux system with Bash. The following dependencies
need to be satisfied.

## 1 Install the following dependencies

### Perl dependencies

Perl 5.10 or higher is required. Download the Core Package from
http://www.bioperl.org

The following non-Core Perl modules are required:

* `Bio::Perl`
* `Bio::SeqIO`
* `List::Util`

These modules are available at CPAN and can be installed/updated with

    cpan Bio::Perl Bio::SeqIO List::Util

### HMMER

HMMER3 (available at http://hmmer.org/)

### Augustus

Augustus3 (available at http://bioinf.uni-greifswald.de/augustus/binaries/)

### GlimmerHMM

GlimmerHMM (available at http://ccb.jhu.edu/software/glimmerhmm/)

### Clustal Omega

Clustal Omega (available at http://www.clustal.org/omega/)

### Blast+
ncbi-blast+ (available at https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

## 2 Download the Pfam
    cd your_path_to_Gene_Modeller
    perl lib.pl

## 3 Add the path to environment:
    export Gene_Modeller=/your_path_to_Gene_Modeller/bin
    export Gene_Modeller_lib=/your_path_to_Gene_Modeller/lib

# Usage

## Input

* Prepare the Genomic sequence file(s) and put it/them into a folder.
* Genomic sequence from the target species in FASTA format and named with its specie name linked with "_".
(For example Arabidopsis_thaliana.fna, Zea_mays.fasta)

The tool is applicable to complete as well as draft genome assemblies. Every sequence needs to have a unique ID 
(first word of a FASTA header is used for ID). Examples of valid FASTA headers:

    >contig10
    sequence(ATGCT...)
    > seq3  genome Z
    sequence((ATGCT...))
    >IV contig 25
    sequence((ATGCT...))

## Running Gene_Modeller

To run Gene_Modeller, use the following command:

    perl $Gene_Modeller/Gene_Modeller.pl your_path_to_genomic_Folder CPU


## Output
* The annotated proteins: your_path_to_genomic_Folder/Resultdir/Gene_Modeller.pep
* The divergence score of each protein compared with PUGTDB: your_path_to_genomic_Folder/Resultdir/Gene_Modeller.txt 

## About
Yuanqian Liu, Xiaoping Liao, Qian Wang, Zihe Li, Zhigang Li, Jian Cheng, Huifeng Jiang.

Tianjin Institute of Industrial Biotechnology, Chinese Academy of Sciences, China.

Reference: []

Contact: jianglab_admin@tib.cas.cn.
