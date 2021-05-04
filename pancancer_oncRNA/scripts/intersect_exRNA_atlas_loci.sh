#!/usr/bin/bash
out=$(basename $1 .srt.filter.bed).exRNA.intersect.bed
echo "intersectBed -s -u -a data/tcga_smRNAs_loci_map.bed -b $1 > data/exRNA/$out"
intersectBed -s -u -a data/tcga_smRNAs_loci_map.bed -b $1 > data/exRNA/$out
