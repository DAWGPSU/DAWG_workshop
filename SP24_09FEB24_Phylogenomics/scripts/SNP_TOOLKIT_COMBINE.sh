#!/bin/bash

set -uex

snptoolkit combine --loc NZ_CP017039.1 -r 0.9

# --loc the name of chromosome or plasmid you want to concatenate for SNPs; this can be found in the last column of the output file of the annotate command  
# --bam this option takes three parameters in the following order: depth ratio patho_to_bam_files
# --snps type of SNPs to be concatenated (default=all)
# -r = M/M+R where M is the number of reads that carry the mutated allele and R is the number of reads that carry the reference allele; if not specified, all SNPs will be taken into account 
# -e specifies a yaml file with two key arguments KEYWORDS and COORDINATES
