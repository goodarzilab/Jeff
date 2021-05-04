#!/usr/bin/python3
import os
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import pandas as pd
import re
import pysam

project_dir = "/rumi/shams/jwang/TCGA_oncRNA"
tcga_metadata = pd.read_csv(f"{project_dir}/TCGA_metadata.csv").set_index("UUID")
tcga_files = [f.path for f in os.scandir("/rumi/shams/jwang/TCGA") if f.is_dir()]
print("Processing", len(tcga_files), "files.")

"""Simple dust score adapted from Dust Algorithm Score: A fast and symmetric DUST implementation to mask
low-complexity DNA sequences (2006, Morgulis et al.)"""
def simpleDustScore(seq):
    assert len(seq) > 2
    if len(seq) == 3:
        return 0
    else:
        triplets = {}
        num_trip = len(seq) - 2
        for i in range(num_trip):
            subseq = seq[i:i+3]
            if subseq in triplets:
                triplets[subseq] += 1
            else:
                triplets[subseq] = 1
        sum_triplet = 0
        for triplet, count in triplets.items():
            sum_triplet += count * (count - 1) / 2
        return sum_triplet/(num_trip - 1)

def processInput(f):
    UUID = re.split("/", f)[-1]
    file = tcga_metadata.loc[UUID]
    study_code = file["TCGAbarcode"]
    subf = [f.path for f in os.scandir(f) if f.name.endswith(".bam")]
    subf = subf[0]
    print("Process:", study_code)
    out = f"{project_dir}/data/{study_code}.dust.bam"
    infile = pysam.AlignmentFile(subf, "rb")
    outfile = pysam.AlignmentFile(out, "wb", template=infile)
    for read in infile.fetch():
        if simpleDustScore(read.get_forward_sequence()) < 3 and len(read.get_forward_sequence()) >= 15:
            outfile.write(read)
    outfile.close()
    infile.close()
    print("Finished:", study_code)
    return f"{study_code} done"

results = Parallel(n_jobs=30)(delayed(processInput)(i) for i in tcga_files)
print("Finished Processing", len(results), "files")