#!/usr/bin/bash
out=${1/.filter.bed/.test.bed}
echo "intersectBed -s -wo -a $1 -b /rumi/shams/jwang/BRCA_oncRNA/results/oldall_cell_lines_sig_loci.bed > $out"
intersectBed -s -wo -a $1 -b /rumi/shams/jwang/BRCA_oncRNA/results/oldall_cell_lines_sig_loci.bed > $out