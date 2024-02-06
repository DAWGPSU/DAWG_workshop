#!/bin/bash

set -uex

module load samtools bwa anaconda3

ls *.fastq|while read samples;
do NAME=${samples};
samtools view -bSh ${samples}_bwa_all.sam > ${samples}_bwa_all.bam;
echo " " >> ${samples}_MappingReport.txt;
echo "${NAME}" >> ${samples}_MappingReport.txt;
echo "---------------------------------------" >> ${samples}_MappingReport.txt;
echo "Analysis-ready reads" >> ${samples}_MappingReport.txt;
samtools view -c ${samples}_bwa_all.bam >> ${samples}_MappingReport.txt;
samtools view -bh -F4 ${samples}_bwa_all.bam > ${samples}_bwa_mapped.bam;
samtools view -bh -q 37 ${samples}_bwa_mapped.bam > ${samples}_bwa_mapped_q37.bam;
samtools sort -o ${samples}_bwa_mapped_q37_sorted.bam ${samples}_bwa_mapped_q37.bam;
echo "Mapped reads" >> "${samples}"_MappingReport.txt;
samtools view -c ${samples}_bwa_mapped.bam >> ${samples}_MappingReport.txt;
echo "Q37 mapped reads" >> ${samples}_MappingReport.txt;
samtools view -c ${samples}_bwa_mapped_q37_sorted.bam >> ${samples}_MappingReport.txt;
samtools rmdup -s ${samples}_bwa_mapped_q37_sorted.bam ${samples}_bwa_mapped_q37_sorted_rmdup.bam;
echo "After removing duplicates" >> ${samples}_MappingReport.txt;
samtools view -c ${samples}_bwa_mapped_q37_sorted_rmdup.bam >> ${samples}_MappingReport.txt;
samtools view -h ${samples}_bwa_mapped_q37_sorted_rmdup.bam |grep -v 'XT:A:R'|grep -v 'XA:Z' |grep -v 'XT:A:M' |awk '{if($0~/X1:i:0/||$0~/^@/)print $0}' | samtools view -bS - > ${samples}_bwa_mapped_q37_sorted_rmdup_uniq.bam;
echo "Unique reads" >> ${samples}_MappingReport.txt;
samtools view -c ${samples}_bwa_mapped_q37_sorted_rmdup_uniq.bam >> ${samples}_MappingReport.txt;
echo "Average length of mapped reads" >> ${samples}_MappingReport.txt;
samtools view ${samples}_bwa_mapped_q37_sorted_rmdup_uniq.bam |awk '{SUM+=length($10);DIV++}END{print SUM/DIV}' >> ${samples}_MappingReport.txt;
no_unique_reads=`samtools view -c ${samples}_bwa_mapped_q37_sorted_rmdup_uniq.bam`;
samtools index ${samples}_bwa_mapped_q37_sorted_rmdup_uniq.bam;
samtools flagstat ${samples}_bwa_mapped_q37_sorted_rmdup_uniq.bam >> ${samples}_flagstat.txt;
done
