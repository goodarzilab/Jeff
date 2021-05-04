#!/usr/bin/bash
out=${1/.filter.bed/.loci.bed}
echo "intersectBed -s -wo -a $1 -b  /rumi/shams/jwang/TCGA_oncRNA/data/tcga_smRNAs_loci_map_exRNA_filtered.bed > $out"
intersectBed -s -wo -a $1 -b  /rumi/shams/jwang/TCGA_oncRNA/data/tcga_smRNAs_loci_map_exRNA_filtered.bed > $out