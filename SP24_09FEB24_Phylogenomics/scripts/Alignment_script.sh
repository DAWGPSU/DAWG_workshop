#!/bin/bash

# This shell script takes in fastq files within a directory and aligns them to an oral microbial reference
# It generates a filtered, sorted, and uniq reads bam file

set -uex

module load bwa
module load samtools
module load anaconda3

# Change path to where you saved the Reference genome
REF=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/SterlingLWright/DAWG_Phylogenomics/reference/Anaerolineacea_bacterium_oral_taxon_439.fasta

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
bwa aln -l 1000 -n 0.01 $REF $samples > $samples.sai
bwa samse $REF $samples.sai $samples > $samples_bwa_all.sam
samtools view -bSh $samples_bwa_all.sam > $samples_bwa_all.bam
echo " " >> $samples_MappingReport.txt
echo "${NAME}" >> $samples_MappingReport.txt
echo "---------------------------------------" >> $samples_MappingReport.txt
echo "Analysis-ready reads" >> $samples_MappingReport.txt
samtools view -c $samples_bwa_all.bam >> $samples_MappingReport.txt
samtools view -bh -F4 $samples_bwa_all.bam > $samples_bwa_mapped.bam
samtools view -bh -q 37 $samples_bwa_mapped.bam > $samples_bwa_mapped_q37.bam
samtools sort -o $samples_bwa_mapped_q37_sorted.bam $samples_bwa_mapped_q37.bam
echo "Mapped reads" >> $PATHS/Sample_"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
samtools view -c $samples_bwa_mapped.bam >> $samples_MappingReport.txt
echo "Q37 mapped reads" >> $samples_MappingReport.txt
samtools view -c $samples_bwa_mapped_q37_sorted.bam >> $samples_MappingReport.txt
samtools rmdup -s $samples_bwa_mapped_q37_sorted.bam $samples_bwa_mapped_q37_sorted_rmdup.bam
echo "After removing duplicates" >> $samples_MappingReport.txt
samtools view -c $samples_bwa_mapped_q37_sorted_rmdup.bam >> $samples_MappingReport.txt
samtools view -h $samples_bwa_mapped_q37_sorted_rmdup.bam |grep -v 'XT:A:R'|grep -v 'XA:Z' |grep -v 'XT:A:M' |awk '{if($0~/X1:i:0/||$0~/^@/)print $0}' | samtools view -bS - > $samples_bwa_mapped_q37_sorted_rmdup_uniq.bam
echo "Unique reads" >> $samples_MappingReport.txt
samtools view -c $samples_bwa_mapped_q37_sorted_rmdup_uniq.bam >> $samples_MappingReport.txt
echo "Average length of mapped reads" >> $samples_MappingReport.txt
samtools view $samples_bwa_mapped_q37_sorted_rmdup_uniq.bam |awk '{SUM+=length($10);DIV++}END{print SUM/DIV}' >> $samples_MappingReport.txt
no_unique_reads=`samtools view -c $samples_bwa_mapped_q37_sorted_rmdup_uniq.bam`
samtools index $samples_bwa_mapped_q37_sorted_rmdup_uniq.bam
samtools flagstat $samples_bwa_mapped_q37_sorted_rmdup_uniq.bam >> $samples_flagstat.txt
done

