use strict;
use warnings;

`wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`;
`gunzip -f Pfam-A.hmm.gz`;
`mv Pfam-A.hmm lib/`;

`wget http://pfam.xfam.org/family/PF00201/hmm`;
`mv hmm lib/HMM/PF00201.hmm`;