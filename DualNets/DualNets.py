#!/usr/bin/env python
# coding: utf-8

#Input files must be .txt and tab-delimited with an empty line at the end. 
#Columns are the samples' ID and rows are the values (gene expression / stability). So the first row of two files are the same (sample IDs). Next lines begin with the gene names (ENSG format). The gene names in mRNA stability file must be preceded by 'mRNAStab_'. 

import numpy as np
from scipy import stats
import pandas as pd
import json as js
import MI
import random
import subprocess
import sys
from multiprocessing import Pool
import matplotlib.pyplot as plt
import time
from datetime import datetime
import argparse
import os

subprocess.call('clear', shell=True)

def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--cpu", help="Number of cpus to use for DualNets", type=int)
    parser.add_argument("--netA", help="Expression file name for networkA", type=str)
    parser.add_argument("--netB", help="Expression file name for networkB", type=str)
    parser.add_argument("--output", help="Output directory path", type=str)

    parser.set_defaults(cpus=1)
    args = parser.parse_args()
    return args


def check_args(args):
    if args.cpu and args.netA and args.netB and args.output:
        return True
    return False

def preprocess_data(netA_path, netB_path, CPU_num):
    with open(netA_path, 'r') as i:
        lines = i.readlines()
    with open(netB_path, 'r') as i:
        lines2 = i.readlines()[1:]
    alllines = lines + lines2

    with open("combined_networks_data.txt", 'w') as output:
        for i in range(0, len(alllines)):
            output.write(alllines[i])
    del alllines

    data = pd.read_csv("combined_networks_data.txt", index_col=0, sep='\t')

    ### Filtering out data based on threshold std (here 0 for zero variation; see notebook for more explanation on this)
    threshold_std = 0
    names = data.T.std()
    names = names[names > threshold_std].index
    data = data.loc[names,]
    print('STD threshold applied to data @', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    ### Winsorization
    for i in range(0, len(data.index)):
        data.iloc[i,] = stats.mstats.winsorize(data.iloc[i,], limits=0.01)

    ## Add noise and rank transform & save data into txt file
    for i in range(0, len(data.index)):
        data.iloc[i,] = data.iloc[i,] + np.random.normal(loc=0, scale=.00001, size=len(data.columns))
        data.iloc[i,] = data.iloc[i,].argsort().argsort()
    print('Data winsorized, noised, and rank transformed')
    data.to_csv('combined_networks_ranked_filtered_data.txt', sep='\t')

    # Now, we create a dictionary with names as keys and expression data sequence as values. This dictionary will be used over and over again in this notebook, so we write it down in `json` format into `'exp_dic_ranked.json'`. Later, this file will be used by other subprocesses in the notebook; so if you want to change the name of this file make sure to  change it in other places, inculding `MI_X.py` files
    dic_genes_exp = data.T.to_dict(orient='list')
    with open('exp_dic_ranked.json', 'w') as f:
        js.dump(dic_genes_exp, f)
        f.close()
    print('Dictionary created and saved @', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    # Here, we extract names of genes in gene-expression and mRNA stability data and write them down into two separate files. Then we generate ALL possible pairs for MI calculation. The order of the pairs in the `names_pairs.txt` is as: gene.expression-mRNA.Stability, gene.expression-gene.expression, and mRNA.Stability-mRNA.Stability. Generally speaking, if you comapre network `X` with network `Y`, the arrangement is `X-Y`, `X-X`, and `Y-Y`.
    # As the number of pairs is too high, we consequently break down the `names_pairs.txt` into `CPU_percent*NCPU+1` parts, where the `NCPU` is number of CPUs in the running machine. This is required for the following step which is MI calculation.

    names = list(data.index)
    position = [i for i in range(len(names)) if names[i].startswith("NetB")][0]
    netA_names = names[0:position]
    netB_names = names[position:]

    with open('netA_names.txt', 'w') as networkA:
        for each in netA_names:
            networkA.write(each + '\n')
        networkA.close()
    with open('netB_names.txt', 'w') as networkB:
        for each in netB_names:
            networkB.write(each + '\n')
        networkB.close()

    # Generate an index file with one pair of names of genes at each row
    with open("name_pairs.txt", "w") as names:
        total_pairs = 0
        for A_name in netA_names:
            for B_name in netB_names:
                names.write(A_name + "\t" + B_name + "\n")
                total_pairs += 1

        for i in range(len(netA_names)):
            for j in range(len(netA_names)):
                if i > j:
                    names.write(netA_names[i] + "\t" + netA_names[j] + "\n")
                    total_pairs += 1

        for i in range(len(netB_names)):
            for j in range(len(netB_names)):
                if i > j:
                    names.write(netB_names[i] + "\t" + netB_names[j] + "\n")
                    total_pairs += 1
    print("All possible name pairs(", total_pairs, ")were generated and saved @",
              datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    with open('name_pairs.txt', 'r') as names:
        bins = CPU_num
        width = total_pairs // bins
        line = names.readline()
        part = 1
        lines_written = 0
        write_part = open(f"part{part}.txt", "w")
        while line:
            if lines_written <= width or part == bins:  # This forces the last [part.txt] to take the remaining lines so that we have exactly [bins] amount of part files.
                write_part.write(line)
                lines_written += 1
            else:
                write_part.close()
                part += 1
                write_part = open(f"part{part}.txt", "w")
                write_part.write(line)
                lines_written = 1
            line = names.readline()
        write_part.close()
        names.close()

    print('Name pairs split into %s parts and saved @' % (bins), datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    #Create shuffled expression dictionary for background MI Calculation
    with open("exp_dic_ranked.json", "r") as f:
        shuffle_dic_genes_exp = js.load(f)
        f.close()

    for name, exprank_list in shuffle_dic_genes_exp.items():
        random.shuffle(exprank_list)

    with open("shuffle_exp_dic_ranked.json", "w") as f:
        js.dump(shuffle_dic_genes_exp, f)
        f.close()


def MI_calc(n):
    print('MI calculation %s began' %n)
    with open('exp_dic_ranked.json', 'r') as f:
        dic_genes_exp = js.load(f)
        f.close()
    partfile = f"part{n}.txt"
    outfile = f"MI_pairs{n}.txt"
    with open(partfile, 'r') as index, open(outfile, "w") as out:
        start_time = time.time()
        header = "nameA\tnameB\tMI\n"
        out.write(header)
        for each in index:
            names = each.strip("\n").split("\t")
            nameA = names[0]
            nameB = names[1]
            exp_A = dic_genes_exp[nameA]
            exp_B = dic_genes_exp[nameB]
            mi = MI.computeMI_outer(exp_A, exp_B)
            out.write(f"{nameA}\t{nameB}\t{mi}\n")
        end_time = time.time()
        exe_time = end_time - start_time
        print(f"MI calculation {n} done:{exe_time}")
        index.close()
        out.close()


#Calculate Background MI
def background_MI_calc(n):
    print('Background MI calculation %s began' % n)
    with open('shuffle_exp_dic_ranked.json', 'r') as f:
        dic_genes_exp = js.load(f)
        f.close()

    partfile = f"part{n}.txt"
    outfile = f"MI_random_pairs{n}.txt"
    with open(partfile, 'r') as index, open(outfile, "w") as out:
        start_time = time.time()
        header = "nameA\tnameB\tMI\n"
        out.write(header)
        for each in index:
            names = each.strip("\n").split("\t")
            nameA = names[0]
            nameB = names[1]
            exp_A = dic_genes_exp[nameA]
            exp_B = dic_genes_exp[nameB]
            mi = MI.computeMI_outer(exp_A, exp_B)
            out.write(f"{nameA}\t{nameB}\t{mi}\n")
        end_time = time.time()
        exe_time = end_time - start_time
        print(f"Background MI calculation {n} done:{exe_time}")
        index.close()
        out.close()

def process_FDR(bins):
    #TODO can precalculate how many pairs for each and then preallocate space for using np array so we don't have to waste time coverting to np.array later.
    MI_netA_A = []
    MI_random_netA_A = []
    MI_netB_B = []
    MI_random_netB_B = []
    MI_netA_B = []
    MI_random_netA_B = []

    for i in range(1, bins + 1):
        print("Start part", i)
        with open(f"MI_pairs{i}.txt", "r") as pairs, open(f"MI_random_pairs{i}.txt", "r") as random_pairs:
            line = pairs.readline()  # This reads the header
            line = pairs.readline()
            random_line = random_pairs.readline()  # This reads the header
            random_line = random_pairs.readline()
            while line:
                parts = line.strip("\n").split("\t")
                random_parts = random_line.strip("\n").split("\t")
                if "NetB" in parts[0] and "NetB" in parts[1]:
                    MI_netB_B.append(float(parts[2]))
                    MI_random_netB_B.append(float(random_parts[2]))
                elif "NetB" in parts[0] or 'NetB' in parts[1]:
                    MI_netA_B.append(float(parts[2]))
                    MI_random_netA_B.append(float(random_parts[2]))
                else:
                    MI_netA_A.append(float(parts[2]))
                    MI_random_netA_A.append(float(random_parts[2]))
                line = pairs.readline()
                random_line = random_pairs.readline()
            pairs.close()
            random_pairs.close()

    MI_netA_A_array = np.array(MI_netA_A)
    del MI_netA_A
    MI_random_netA_A_array = np.array(MI_random_netA_A)
    del MI_random_netA_A
    MI_netA_B_array = np.array(MI_netA_B)
    del MI_netA_B
    MI_random_netA_B_array = np.array(MI_random_netA_B)
    del MI_random_netA_B
    MI_netB_B_array = np.array(MI_netB_B)
    del MI_netB_B
    MI_random_netB_B_array = np.array(MI_random_netB_B)
    del MI_random_netB_B

    # Sort MIs
    MI_netA_A_array = np.sort(MI_netA_A_array)
    MI_random_netA_A_array = np.sort(MI_random_netA_A_array)
    netA_A_results = calculate_FDR(MI_netA_A_array, MI_random_netA_A_array, "NetA_A")


    MI_netA_B_array = np.sort(MI_netA_B_array)
    MI_random_netA_B_array = np.sort(MI_random_netA_B_array)
    netA_B_results = calculate_FDR(MI_netA_B_array, MI_random_netA_B_array, "NetA_B")

    MI_netB_B_array = np.sort(MI_netB_B_array)
    MI_random_netB_B_array = np.sort(MI_random_netB_B_array)
    netB_B_results = calculate_FDR(MI_netB_B_array, MI_random_netB_B_array, "NetB_B")

    return netA_A_results, netA_B_results, netB_B_results


def calculate_FDR(MI_network_array, MI_random_array, network_name):
    list_FDR = []
    list_N = []
    list_Th = []
    list_percent = []
    percentiles = np.linspace(0.95, 1, 1000)
    array_length = MI_network_array.shape[0]
    indexes = np.floor(percentiles * array_length)

    for i in indexes:
        Th = MI_random_array[int(i - 1)]
        m = binary_search(Th, MI_network_array, array_length)

        if MI_network_array[m] > Th:
            if MI_network_array[m - 1] <= Th:
                j = m
            else:
                assert False  # Sanity Check
        elif MI_network_array[m] < Th:
            if MI_network_array[m + 1] >= Th:
                j = m + 1
            else:
                assert False  # Sanity Check
        else:
            j = m

        N = array_length - j
        FP = array_length - i
        if N > 0:
            FDR = FP / N
            percent = N / array_length
            list_FDR.append(FDR)
            list_N.append(N)
            list_Th.append(Th)
            list_percent.append(percent)

    network_FDR_results = pd.DataFrame({"FDR": list_FDR, "N": list_N, "Th": list_Th, "Proportion of data" : list_percent})
    network_FDR_results = network_FDR_results.sort_values("FDR", ascending=False).reset_index(drop=True)
    network_FDR_results.to_csv(f"{network_name}_FDR.csv")
    return network_FDR_results

def plot_FDR(dataframe, network_name):
    # plt.plot(dataframe["FDR"], np.log10(dataframe["N"]))
    plt.plot(dataframe["FDR"], np.log10(dataframe["N"]))

    plt.xlim(0, 0.1)
    plt.ylabel("log10(N)")
    plt.xlabel("FDR")
    plt.title(f"ROC: {network_name}")
    plt.savefig(f"ROC_{network_name}.png")
    plt.clf()

def binary_search(Th, array, length):
    L = 0
    R = length - 1
    while L <= R:
        m = (L + R) // 2
        if array[m] < Th:
            L = m + 1
        elif array[m] > Th:
            R = m - 1
        else:
            break
    return m


def main():
    args = handler()

    if not check_args(args):
        sys.exit("Missing or incompatible arguments")

    os.chdir(args.output)
    print('Script launched at:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'), ' with the following parameters:')
    print('Number of cpus to use=', args.cpu)

#Preprocess Data
    preprocess_data(args.netA, args.netB, args.cpu)

#MI Calculations
    print('MI calculation initiated @', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    pool = Pool(processes=args.cpu)
    pool.map(MI_calc, [i+1 for i in range(args.cpu)])
    print(' ALL MIs CALCULATIONS DONE @', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    print('Background MI calculation initiated @', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    pool = Pool(processes=args.cpu)
    pool.map(background_MI_calc, [i + 1 for i in range(args.cpu)])
    print(' ALL BACKGROUND MIs CALCULATIONS DONE @', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

#FDR Calculation

    print('FDR calculation initiated')
    netA_A_results, netA_B_results, netB_B_results = process_FDR(args.cpu)

#Plot ROC Curves
    print("Plot ROC Curves")
    plot_FDR(netA_A_results, "NetA_A")
    plot_FDR(netA_B_results, "NetA_B")
    plot_FDR(netB_B_results, "NetB_B")
    print("DONE")

if __name__ == "__main__":
    main()
