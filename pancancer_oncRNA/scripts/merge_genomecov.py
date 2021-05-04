#!/usr/bin/python3
#Script to filter all TCGA bed data
import os
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import pandas as pd
import re
import json as js
from operator import add

project_dir = "/rumi/shams/jwang/TCGA_oncRNA"
total_inputs = [f.path for f in os.scandir(f"{project_dir}/data/TCGA/") if f.name.endswith(".split.genomecov.json")]
print("Processing", len(total_inputs), "files.")

def processInput(inputs): 
    loci_to_split = pd.read_csv("data/tcga_merged_loci_to_split.bed", header=None, sep="\t")
    loci_coverage_dict = {}
    for locus, start, end in zip(loci_to_split[3],loci_to_split[1], loci_to_split[2]):
        length = end-start
        loci_coverage_dict[locus] = [0 for i in range(length)]
    
    for i in inputs:
        with open(i, "r") as f:
            sample_dict = js.load(f)
            
        for locus in loci_to_split[3]:
            loci_coverage_dict[locus] = list(map(add, loci_coverage_dict[locus], sample_dict[locus]))
            
    print("Processed:", len(inputs),"files")        
    return loci_coverage_dict


def batch(ls, n):
    for i in range(0, len(ls), n):
        yield ls[i:i+n]
        
batched_inputs = [i for i in batch(total_inputs, 555)]
results = Parallel(n_jobs=20)(delayed(processInput)(i) for i in batched_inputs)

print("Finished Processing", len(results), "Batches")

loci_to_split = pd.read_csv("data/tcga_merged_loci_to_split.bed", header=None, sep="\t")

final_coverage_dict = {}
for locus, start, end in zip(loci_to_split[3],loci_to_split[1], loci_to_split[2]):
    length = end-start
    final_coverage_dict[locus] = [0 for i in range(length)]  
    
for merged_dict in results: #Final merge from 20 results.
    for locus in loci_to_split[3]:
        final_coverage_dict[locus] = list(map(add, final_coverage_dict[locus], merged_dict[locus]))  

with open("data/merged_loci_to_split_coverage.json", "w") as out:
    js.dump(final_coverage_dict, out)

print("Finished Merging Genome Coverage")
