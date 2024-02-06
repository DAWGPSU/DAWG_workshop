#!/bin/bash

set -uex

# Install pip commands
# conda install anaconda::pip
# pip install snptoolkit

REF2=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/SterlingLWright/REFERENCES/Anaerolineacea/GCF_001717545.1_ASM171754v1_genomic.gbff
snptoolkit annotate -i vcf -g $REF2 -q 30 -d 3 -r 0.9 -f 3

# annotate: Annotate one or more vcf files
# -i: Provide a specific identifier to recognize the file(s) to be analyzed
# -q: Quality score to consider as a cutoff for variant calling. default value [20]
# -d: Minimum depth caverage. default value [3]
# -r: Minimum ratio that correspond to the number of reads that has the mutated allele; total depth in that particular position. default value [0]
# -f: Exclude SNPs if the distance between them is lower then the specified window size in bp
