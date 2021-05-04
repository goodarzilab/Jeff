#!/usr/bin/python3
#Script to filter all TCGA bed data
import os
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import pandas as pd
import re
import json as js

project_dir = "/rumi/shams/jwang/TCGA_oncRNA"
tcga_files = [f.path for f in os.scandir(f"{project_dir}/data/TCGA/") if f.name.endswith(".split.bed")]
print("Processing", len(tcga_files), "files.")


def processInput(f): 
    name = f.split("/")[-1].split(".")[0]
    loci_to_split = pd.read_csv(f"{project_dir}/data/tcga_merged_loci_to_split.bed", header=None, sep="\t")
    loci_coverage_dict = {}
    for locus, start, end in zip(loci_to_split[3],loci_to_split[1], loci_to_split[2]):
        length = end-start
        loci_coverage_dict[locus] = [0 for i in range(length)]

    with open(f, "r") as read:
        for line in read:
            splits = line.split("\t")
            locus = splits[9]
            if locus in loci_coverage_dict:
                read_start = int(splits[1])
                read_end = int(splits[2])
                locus_start = int(splits[7])
                locus_end = int(splits[8])
                #This maps to the list index for the read. Range works well with bed format as it is end exclusive.
                for i in range(read_start - locus_start, read_end - locus_start): 
                    loci_coverage_dict[locus][i] += 1
                    
    with open(f'{project_dir}/data/TCGA/{name}.split.genomecov.json', 'w') as out:
        js.dump(loci_coverage_dict, out)
        out.close() 
        
    print(f"{name} Done")
    return f"{name} Done"
        

results = Parallel(n_jobs=20)(delayed(processInput)(i) for i in tcga_files)

print("Finished Processing", len(results), "files")