#!/usr/bin/python3
import pysam
import glob
import argparse
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import pandas as pd
import scipy.stats as stats
import statsmodels
import seaborn as sns
import re
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
import time

def processInput(path):
    samfile = pysam.AlignmentFile(path, "rb")
    data = samfile.fetch()
    sample_name = re.split("\.",re.split("/", path)[-1])[0]
    rna = {}
    total_count = samfile.count()
    cpm = 1000000/total_count
    for read in data:
        seq = read.get_forward_sequence()
        ref_id = read.reference_name
        ref_start = read.reference_start
        ref_end = read.reference_end
        strand = "+"
        if read.is_reverse:
            strand = "-"

        seq_strand = seq + "_" + strand
        if seq_strand in rna:
            rna[seq_strand][1] += cpm
        else:
            bed = f"{ref_id}\t{ref_start}\t{ref_end}"
            rna[seq_strand] = [sample_name, cpm, bed]
    return rna

def aggregateRNA(rna_results):
    agg_rna = {}
    loc_rna = {}
    for rna in rna_results:
        for seq_strand, sample in rna.items():
            sample_name = sample[0]
            cpm = sample[1]
            if seq_strand in agg_rna:
                seq_dict = agg_rna[seq_strand]
                seq_dict[sample_name] = cpm
            else:
                agg_rna[seq_strand] = {sample_name: cpm}
                loc_rna[seq_strand] = sample[2]
    return agg_rna, loc_rna

def filterRNA(agg_rna, dust, thresh):
    to_delete=[]
    for seq_strand, rna in agg_rna.items():
        num_ctrl = 0
        total = 0
        seq = seq_strand[:-2]
        if simpleDustScore(seq) >= dust:
            to_delete.append(seq_strand)
            continue
        for study, v in rna.items():
            total += 1
            if study in ctrl_labels:
                num_ctrl += 1
        if num_ctrl >= thresh or total < 2:
            to_delete.append(seq_strand)
    for seq_strand in to_delete:
        agg_rna.pop(seq_strand)
    return agg_rna

def simpleDustScore(seq):
    assert len (seq) > 2
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
        sum_trip = 0
        for triplet, count in triplets.items():
            sum_trip += count*(count-1)/2
        return sum_trip/(num_trip - 1)


def fisher(row, ctrl_total, test_total):
    """Run one-sided fisher exact test on a row of data. Test vs Control.
    """
    test_count= row["test_count"]
    ctrl_count = row["ctrl_count"]
    oddsratio, pvalue = stats.fisher_exact([[test_count, ctrl_count],[test_total-test_count, ctrl_total-ctrl_count]], "greater")
    return pvalue

def addBed(row):
    return loc_rna[row]



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fisher Test")
    parser.add_argument("--p", help="Path to Bamfiles")
    parser.add_argument("--ctrl", help="Path to labels of control")
    parser.add_argument("--test", help="Path to labels of test")
    parser.add_argument("--dust", help="Dust score threshold, default 2.5")
    parser.add_argument("--t", help="Filter control threshold, default 3")

    args = parser.parse_args()

    data_path, ctrl_label_path, test_label_path = None, None, None
    dust, filter, thresh = 2.5, "count", 3

    if args.p:
        data_path = args.p
    if args.ctrl:
        ctrl_label_path = args.ctrl
    if args.test:
        test_label_path = args.test
    if args.dust:
        dust = float(args.dust)
    if args.t:
        thresh = float(args.t)

    bamfiles =[f for f in glob.glob(data_path + "**/*.bam", recursive=True)]

    print("Begin Process Data", len(bamfiles), "Files.")
    start = time.perf_counter()
    num_cores = cpu_count()
    rna_results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in bamfiles)
    end = time.perf_counter()
    print("Time:", (end-start)/60)

    print("Begin Aggregating Data")
    start = time.perf_counter()
    agg_rna, loc_rna = aggregateRNA(rna_results)
    print(len(agg_rna), "sequences")
    end = time.perf_counter()
    print("Time:", (end-start)/60)

    ctrl_file = open(ctrl_label_path,'r')
    test_file = open(test_label_path,'r')
    ctrl_labels = [line for line in ctrl_file.read().splitlines() if line]
    test_labels = [line for line in test_file.read().splitlines() if line]
    ctrl_file.close()
    test_file.close()
    ctrl_total = len(ctrl_labels)
    test_total = len(test_labels)


    print("Begin Filtering Data")
    start = time.perf_counter()
    agg_rna = filterRNA(agg_rna, dust, thresh)
    print(len(agg_rna), "rnas post filter")
    end = time.perf_counter()
    print("Time:", (end-start)/60)

    print("Begin Creating Count Data Matrix")
    start = time.perf_counter()

    pres_df = np.zeros((len(agg_rna), 2))
    seq_strand_index = list(agg_rna.keys())
    for i in range(len(seq_strand_index)):
        seq_strand = seq_strand_index[i]
        samples = agg_rna[seq_strand]
        for s in samples:
            if s in ctrl_labels:
                pres_df[i, 1] += 1
            elif s in test_labels:
                pres_df[i, 0] += 1
            else:
                assert False
    pres_df = pd.DataFrame(pres_df, index=seq_strand_index, columns=["test_count", "ctrl_count"])
    print(pres_df.shape)
    end = time.perf_counter()
    print("Time:", (end-start)/60)

    #Fisher Part
    print("Begin Fisher Exact Test")
    start = time.perf_counter()
    pres_df["pval"] = pres_df.apply(fisher, args=(ctrl_total, test_total),axis=1)
    end = time.perf_counter()
    print("Time:", (end-start)/60)

    #Fdr Part
    print("Begin FDR correction. Alpha=0.1")
    rej, test_fdr = fdrcorrection(pres_df["pval"], alpha=0.1)
    pres_df["fdr"] = test_fdr
    if not any(rej):
        print("No significant RNAs after FDR correction, plotting RNAs with unadjusted pval <= 0.05")
        sig_df = pres_df[pres_df["pval"] <= 0.05]
    else:
        sig_df = pres_df[rej]
    sig_df = sig_df.sort_values(by="test_count", ascending = False)
    sig_index = sig_df.index
    print(len(sig_index), "significant RNAs found.")

    print("Creating Data Matrix for Plot")
    ctrl_df = np.zeros((len(sig_index), ctrl_total))
    test_df = np.zeros((len(sig_index), test_total))
    for i in range(len(sig_index)):
        for j in range(test_total + ctrl_total):
            if j < ctrl_total:
                if ctrl_labels[j] in agg_rna[sig_index[i]]:
                    ctrl_df[i,j] = agg_rna[sig_index[i]][ctrl_labels[j]]
            else:
                j -= ctrl_total
                if test_labels[j] in agg_rna[sig_index[i]]:
                    test_df[i,j] = agg_rna[sig_index[i]][test_labels[j]]
    ctrl_df = pd.DataFrame(ctrl_df, index=sig_index, columns=ctrl_labels)
    test_df = pd.DataFrame(test_df, index=sig_index, columns=test_labels)

    ctrl_binary = ctrl_df.clip(upper=1)
    ctrl_binary = ctrl_binary.apply(np.ceil)
    test_binary = test_df.clip(upper=1)
    test_binary = test_binary.apply(np.ceil)

    #Visualization Part: Binary
    print("Begin Binary visualization.")
    fig, (ax, ax2) = plt.subplots(ncols=2, gridspec_kw={"width_ratios":[ctrl_total, test_total]})
    fig.subplots_adjust(wspace=0.005)
    sns.heatmap(ctrl_binary.values, yticklabels=False, xticklabels=False, cmap="rocket_r", ax=ax, cbar=False)
    sns.heatmap(test_binary.values, cmap="rocket_r", ax=ax2, yticklabels=False, xticklabels=False)
    ax.set(xlabel="Control", ylabel="RNAs")
    ax2.set(xlabel="Test")
    ax.set_fc("w")
    plt.savefig("Binary Test Pipeline RNA Heatmap.png")


    #Normalize by row:
    max_row = pd.DataFrame([test_df.max(axis=1), ctrl_df.max(axis=1)]).max()
    norm_ctrl_df = ctrl_df.div(max_row, axis=0)
    norm_test_df = test_df.div(max_row, axis=0)

    # Visualization Part: Normalized
    print("Begin Normalized visualization.")
    fig, (ax, ax2) = plt.subplots(ncols=2, gridspec_kw={"width_ratios":[ctrl_total, test_total]})
    fig.subplots_adjust(wspace=0.005)
    sns.heatmap(norm_ctrl_df.values, yticklabels=False, xticklabels=False, cmap="rocket_r", ax=ax, cbar=False)
    sns.heatmap(norm_test_df.values, cmap="rocket_r", ax=ax2, yticklabels=False, xticklabels=False)
    ax.set(xlabel="Control", ylabel="RNAs")
    ax2.set(xlabel="Test")
    ax.set_fc("w")
    plt.savefig("Normalized Test Pipeline RNA Heatmap.png")
    #Bed + Annotation Part
    bed_data = pd.Series(sig_index).apply(addBed)
    bed_data.to_csv("sig_rna.bed", index=False, header=False)

    #Export Data
    test_df["pval"] = sig_df["pval"]
    test_df["fdr"] = sig_df["fdr"]
    test_df.to_csv("sig_rna_test_expression.csv")
    ctrl_df.to_csv("sig_rna_ctrl_expression.csv")
    print("Complete Run.")
