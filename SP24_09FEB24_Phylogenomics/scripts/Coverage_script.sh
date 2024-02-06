#!/bin/bash

set -uex

ls *uniq.bam|while read samples; do echo "----------------------------------------------------------------------" >> Anaero_samples_coverage.txt;
echo $samples >> Anaero_samples_coverage.txt; samtools coverage $samples >> Anaero_samples_coverage.txt;
echo "----------------------------------------------------------------------" >> Anaero_samples_coverage.txt;
echo "----------------------------------------------------------------------" >> Anaero_histo_coverage.txt;
echo $samples >> Anaero_histo_coverage.txt;
samtools coverage -A -w 32 $samples >> Anaero_histo_coverage.txt; echo "----------------------------------------------------------------------" >> Anaero_histo_coverage.txt; done;
