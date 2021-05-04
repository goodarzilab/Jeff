#!/usr/bin/bash
#To parallelize: ls /rumi/shams/jwang/TCGA_oncRNA/data/*.dust.bam | parallel -j 30 bash scripts/bamtobed.sh {} &> log/bamtobed.out
pos_out=${1/.dust.bam/.pos.bg.txt}
neg_out=${1/.dust.bam/.neg.bg.txt}
echo "bedtools genomecov -ibam $1 -bg -strand - > $neg_out"
bedtools genomecov -ibam $1 -bg -strand - > $neg_out

echo "bedtools genomecov -ibam $1 -bg -strand + > $pos_out"
bedtools genomecov -ibam $1 -bg -strand + > $pos_out