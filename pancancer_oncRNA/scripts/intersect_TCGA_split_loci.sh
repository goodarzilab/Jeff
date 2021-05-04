#!/usr/bin/bash
out=${1/.filter.bed/.split.bed}
echo "intersectBed -s -wo -a $1 -b  /rumi/shams/jwang/TCGA_oncRNA/data/tcga_merged_loci_to_split.bed > $out"
intersectBed -s -wo -a $1 -b  /rumi/shams/jwang/TCGA_oncRNA/data/tcga_merged_loci_to_split.bed > $out;