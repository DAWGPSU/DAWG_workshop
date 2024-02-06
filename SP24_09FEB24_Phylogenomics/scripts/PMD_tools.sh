
# Install PMD tools
# This tool includes only reads that exhibit post mortem damage
# conda install bioconda::pmdtools
for i in *.bam; do samtools view -h $i |pmdtools --threshold 3 --header |samtools view -Sb > ${i}.pmds3filter.bam;
echo ${i}.pmds3filter.bam >> PMD3Filter_Counts.txt;
samtools view -c ${i}.pmds3filter.bam >> PMD3Filter_Counts.txt;
echo "______________________________________________________________________________" >> PMD3Filter_Counts.txt;
done
