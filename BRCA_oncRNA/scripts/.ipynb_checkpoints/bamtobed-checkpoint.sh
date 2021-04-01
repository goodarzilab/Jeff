#!/usr/bin/bash
#To parallelize: ls /rumi/shams/jwang/BRCA_oncRNA/data/TCGA/*.dust.bam | parallel -j 30 bash scripts/bamtobed.sh {} &> log/bamtobed.out
out=${1/.dust.bam/.bed}
echo "bamToBed -i $1 > $out"
bamToBed -i $1 > $out
