#!/bin/bash

# This shell script takes in fastq files within a directory and aligns them to an oral microbial reference
# It generates a filtered, sorted, and uniq reads bam file

set -uex

module load bwa
module load samtools
module load anaconda3

# Change path to where you saved the Reference genome
REF=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/SterlingLWright/DAWG_Phylogenomics/reference/Anaerolineacea_bacterium_oral_taxon_439.fasta

# Change path to where you have the fastq files
PATH=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/SterlingLWright

# Index the reference
bwa index $REF
samtools faidx $REF

##########################################################################################################################
#
#
# No edits below necessary #
#
#########################################################################################################################

ls *.fastq | while read samples
do
bwa aln -l 1000 -n 0.01 $REF $samples > $samples.sai;
ls *.fastq|while read samples; do bwa samse $REF ${samples}.sai ${samples} > ${samples}_bwa_all.sam; done

