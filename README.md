# PUGTA
A tool for mining of plant UGTs.

# Overview

PUGTA is a command line program (in Perl) based on Augustus, GlimmerHMM, 
which is used for sufficiently mining of UDP-glycosyltransferases from plant genomes.


# Installation

Running PUGTA requires a Linux system with Bash. The following dependencies
need to be satisfied.

### Perl dependencies

Perl 5.10 or higher is required.

The following non-Core Perl modules are required:

* `Bio::Perl`
* `Bio::SeqIO`
* `List::Util`

These modules are available at CPAN and can be installed/updated with

    cpan Bio::Perl Bio::SeqIO List::Util

### Augutus

Augustus (available at http://bioinf.uni-greifswald.de/augustus/binaries/)

### GlimmerHMM

GlimmerHMM (available at http://ccb.jhu.edu/software/glimmerhmm/)

### Clustal Omega

Clustal Omega (available at http://www.clustal.org/omega/)

### Blast+
ncbi-blast+ (available at https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)


# Usage

## Input

* Genomic sequence from the target species in FASTA format

The tool is applicable to complete as well as draft genome assemblies. Every sequence needs to have a unique ID 
(first word of a FASTA header is used for ID). Examples of valid FASTA headers:

    >contig10
    sequence: contig10
    > seq3  genome Z
    sequence: seq3
    >IV contig 25
    sequence: IV

## Running PUGTA

To run PUGTA, use the following command:

    Genome_Annoation.pl TargetDir CPU Pfam Query_DB


## Output

## About
Yuanqian Liu, Jian Cheng, Huifeng Jiang

Tianjin Institute of Industrial Biotechnology, Chinese Academy of Sciences, China

Reference: []
