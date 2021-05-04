#!/usr/bin/bash
out=${1/.dust.bam/.bed}
echo "bamToBed -i $1 > $out"
bamToBed -i $1 > $out
