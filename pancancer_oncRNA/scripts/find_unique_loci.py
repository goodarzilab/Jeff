#!/usr/bin/python3
import os
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import pandas as pd

def processInput(f):
    print("Start", f)
    unique_loci = set()
    with open(f, "r") as read:
        for line in read:
            line = line.strip() #Get rid of trailing "\n" since we are 
            splits = line.split("\t")
            ref = splits[0]
            start = splits[1]
            end = splits[2]
            strand = splits[5]
            bed_locus = f"{ref}\t{start}\t{end}\t{ref}:{start}-{end}:{strand}\t.\t{strand}"
            unique_loci.add(bed_locus)
    print("Finished", f)
    return unique_loci

def batch(ls, n):
    for i in range(0, len(ls), n):
        yield ls[i:i+n]

total_inputs = [f.path for f in os.scandir("data/TCGA") if f.name.endswith(".filter.bed")]
print("Processing", len(total_inputs), "files.")

all_unique_loci = set()
num_batch = 1
num_processed = 0
for inputs in batch(total_inputs, 25):
    results = Parallel(n_jobs=25)(delayed(processInput)(i) for i in inputs)
    for unique_loci in results:
        all_unique_loci.update(unique_loci)
    
    print("Complete batch:", num_batch)
    num_batch += 1
    num_processed += len(results)
    
with open("data/tcga_all_unique_loci.bed", "w") as out:
    for locus in all_unique_loci:
        out.write(locus + "\n")
    
print("Finished Processing", num_processed, "files")

