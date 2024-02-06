#!/bin/bash

set -uex

# ModelFinder indicated that K3Pu+F+ASC was the best-fit model for our example dataset

# Now run iqtree with bootstrapping
iqtree -s SNPs_alignments_95percent-filter.fasta -m K3Pu+F+ASC -alrt 1000 -B 1000 -T AUTO
