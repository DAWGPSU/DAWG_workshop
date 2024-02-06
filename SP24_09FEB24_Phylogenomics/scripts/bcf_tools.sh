#!/bin/bash

set -uex

module load bwa anaconda3 gcc/8.3.1 samtools/1.13

for i in *pmds1filter.bam; do
REF=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/SterlingLWright/DAWG_Phylogenomics/reference/Anaerolineacea_bacterium_oral_taxon_439.fasta ;
bcftools mpileup -Ou -f $REF $i | bcftools call -mv -Ob -o ${i}.bcf; bcftools view ${i}.bcf |vcfutils.pl varFilter -D400 > ${i}.vcf ;
done
