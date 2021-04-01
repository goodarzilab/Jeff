#!/usr/bin/env python
import argparse
import sys
import os
import pandas as pd
import numpy as np

def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Th_AA", help="FDR threshold for networkAA", type=float)
    parser.add_argument("--Th_BB", help="FDR threshold for networkBB", type=float)
    parser.add_argument("--Th_AB", help="FDR threshold for networrkAB", type=float)
    parser.add_argument("--dpi", help="DPI epsilon to use", type=float)
    parser.add_argument("--dir", help="Directory path with all the MI calculations", type=str)
    parser.set_defaults(dpi=0.15)
    args = parser.parse_args()
    return args

def check_args(args):
    if args.Th_AA and args.Th_BB and args.Th_AB and args.dpi and args.dir:
        return True
    return False


def filter_MI(Th_AA, Th_AB, Th_BB, num_parts):
    with open("filtered_MI_netAA.txt", "w") as AA_out, open("filtered_MI_netAB.txt", "w") as AB_out, open(
            "filtered_MI_netBB.txt", "w") as BB_out:
        header = "NameA\tNameB\tMI\n"
        AA_out.write(header)
        AB_out.write(header)
        BB_out.write(header)
        for i in range(1, num_parts + 1):
            print("Start part", i)
            with open(f"MI_pairs{i}.txt", "r") as pairs:
                line = pairs.readline()  # This reads the header
                line = pairs.readline()  # First line with data
                while line:
                    parts = line.strip("\n").split("\t")
                    if "NetB" in parts[0] and "NetB" in parts[1]:
                        if float(parts[2]) >= Th_BB:
                            BB_out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")
                    elif "NetB" in parts[0] or 'NetB' in parts[1]:
                        if float(parts[2]) >= Th_AB:
                            AB_out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")
                    else:
                        if float(parts[2]) >= Th_AA:
                            AA_out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")
                    line = pairs.readline()
                pairs.close()

        AA_out.close()
        AB_out.close()
        BB_out.close()

def build_network_dict():
    MI_AA = pd.read_csv("filtered_MI_netAA.txt", sep="\t")
    MI_AB = pd.read_csv("filtered_MI_netAB.txt", sep="\t")
    matrix = pd.concat([MI_AA, MI_AB])

    nameA_list = np.array(matrix["NameA"])
    nameB_list = np.array(matrix["NameB"])
    MI_list = np.array(matrix["MI"])

    network_dict = {}
    mydict = {}
    nameA_set = set()
    for i in range(len(nameA_list)):
        nameA = nameA_list[i]
        nameB = nameB_list[i]
        netMI = MI_list[i]
        if "NetB" in nameA and "NetB" in nameB: #Here we exlcude all netB-netB pairs
            continue
        network_dict[(nameA, nameB)] = netMI
        mydict[(nameA, nameB)] = netMI
        nameA_set.add(nameA)
    return network_dict, mydict, nameA_set


def DPI(network_dict, mydict, newrows, DPI_epsilon):
    network_series = pd.Series(mydict)
    network_series.index.names = ['nameA', None]

    for each in newrows:
        m = network_series[each].sort_values(ascending=False)
        for i in range(0, len(m)):
            valueAB = mydict[each, m.index[i]]
            minMI = valueAB/(1.00 - DPI_epsilon)
            for j in range(0, i):
                if not("NetB" in m.index[i] and "NetB" in m.index[j]): #no DPI on BC when it is stability-stability
                    valueAC = mydict[each, m.index[j]]
                    if valueAC < minMI:
                        break
                    elif (m.index[i], m.index[j]) in mydict:
                        valueBC = mydict[m.index[i], m.index[j]]
                    elif (m.index[j], m.index[i]) in mydict:
                        valueBC = mydict[m.index[j], m.index[i]]
                    else:
                        continue

                    if valueBC > minMI:
                        del network_dict[each, m.index[i]]
                        break
    return network_dict

def split_df(final_df):
    final_df["nameA"] = final_df["pairs"].apply(lambda x: x[0])
    final_df["nameB"] = final_df["pairs"].apply(lambda x: x[1])
    final_df = final_df[["nameA", "nameB", "MI"]]
    netAA = final_df[(final_df["nameA"].str.contains("NetA")) & (final_df["nameB"].str.contains("NetA"))]
    netAB = final_df[(final_df["nameA"].str.contains("NetA")) & (final_df["nameB"].str.contains("NetB"))]
    netBB = final_df[(final_df["nameA"].str.contains("NetB")) & (final_df["nameB"].str.contains("NetB"))]
    assert netBB.shape[0] == 0

    return netAA, netAB



def main():
    args = handler()
    if not check_args(args):
        sys.exit("Missing or incompatible arguments")

    num_parts = len([f for f in os.scandir(args.dir) if "MI_pairs" in f.name]) #Number of MI pair files

    os.chdir(args.dir)

    filter_MI(args.Th_AA, args.Th_AB, args.Th_BB, num_parts)

    network_dict, mydict, nameA_set = build_network_dict()

    final_dict = DPI(network_dict, mydict, nameA_set, args.dpi)

    final_df = pd.DataFrame.from_dict(final_dict, orient='index', columns=["MI"]).reset_index().rename(columns={"index":"pairs"})

    fileout = 'all_thresholded_DPI_removed_epsilon_' + str(args.dpi) + '.txt'
    final_df.to_csv(fileout, sep='\t')

    netAA_df, netAB_df = split_df(final_df)

    netAA_df.to_csv("final_netAA_results.txt", sep="\t")
    netAB_df.to_csv("final_netAB_results.txt", sep="\t")

    print("Done")

if __name__ == "__main__":
    main()
