#!/usr/bin/python3
#Script to filter all TCGA bed data
import os
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import pandas as pd
import re

project_dir = "/rumi/shams/jwang/TCGA_oncRNA"
tcga_files = [f.path for f in os.scandir(f"{project_dir}/data/") if f.name.endswith(".bed")]
print("Processing", len(tcga_files), "files.")

def processInput(f):
    name = f.split("/")[-1].split(".")[0]
    out = f"{project_dir}/data/{name}.filter.bed"
    print("Start:", name)
    with open(out, "wt") as out, open(f, "rt") as file:
        for line in file:
            if re.match("chr[\d+,X,Y]", line) and "None" not in line:
                out.write(line) 
    print("Finished:", name)
    return f"{name} done"

results = Parallel(n_jobs=30)(delayed(processInput)(i) for i in tcga_files)

print("Finished Processing", len(results), "files")