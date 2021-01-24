# GMind
A novel plant UGTs annotated process for the mining of UGTs in plant genomes

# Overview

GMind is a command line program (in Perl) based on Augustus, GlimmerHMM, which is used for sufficiently mining of UDP-glycosyltransferases from plant genomes. Firstly, the protein sequences of the known plant UGTs from an initial plant UGT database will be mapped into the un-annotated plant genomes, and the mapped regions will be extracted for the next annotation. Secondly, the extracted UGTs regions were independently de nove annotated by: 1) Augustus with the model training from all annotated UGTs in the 195 plant genomes, 2) Augustus with the model species of Arabidopsis thaliana and 3) GlimmerHMM with model species of A. thaliana. Thirdly, HMMER and UGT domain (PF00201) from Pfam were used to filter UGTs from all these annotated proteins. Fourthly, the filtered UGTs were evaluated by calculating accuracy score, the optimal UGTs were selected as candidates, while the UGTs with bad sores would be re-annotated by removing the good UGTs from the mapped regions. Finally, all these candidate UGTs would be output as the result of genomic annotation. 


# Installation

Running GMind requires a Linux system with Bash. The following dependencies
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
    cd your_path_to_GMind
    perl lib.pl

## 3 Add the trained plant UGT model to Augustus path
    cp -r /your_path_to_GMind/lib/PUGT /your_path_to_Augustus/config/species/

## 4 Add the path to environment:
    export GMind=/your_path_to_GMind/bin
    export GMind_lib=/your_path_to_GMind/lib

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

## Running GMind

To run GMind, use the following command:

    perl $GMind/GMind.pl your_path_to_genomic_Folder CPU


## Output
* The annotated proteins: your_path_to_genomic_Folder/Resultdir/GMind.pep
* The divergence score of each protein compared with PUGTDB: your_path_to_genomic_Folder/Resultdir/GMind.txt 

## About
Yuanqian Liu, Xiaoping Liao, Qian Wang, Zihe Li, Zhigang Li, Jian Cheng, Huifeng Jiang.

Tianjin Institute of Industrial Biotechnology, Chinese Academy of Sciences, China.

Reference: []

Contact: jianglab_admin@tib.cas.cn.
