#!/usr/bin/bash
out=${1/.bed/.loci.bed}
echo "intersectBed -s -wo -a $1 -b  /rumi/shams/jwang/BRCA_oncRNA/results/thresholded_oncRNAs.bed > $out"
intersectBed -s -wo -a $1 -b  /rumi/shams/jwang/BRCA_oncRNA/results/thresholded_oncRNAs.bed > $out