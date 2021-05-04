#!/usr/bin/python3
import pymongo
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import pandas as pd
import time
import os

def processInput(f):
    myclient = pymongo.MongoClient(port=27027)
    mydb = myclient["TCGA_loci"]
    rnacol = mydb["smRNAs"]
    tcga_barcode = f.split("/")[-1].split(".")[0]
    
    metadata = pd.read_csv("metadata/TCGA_metadata.csv").set_index("TCGAbarcode", drop=False)
    file = metadata.loc[tcga_barcode]
    sample = file["Sample"]
    tissue = file["TSS"]
    abb = file["TCGA_abb"]
    
    sample_loci = {}
    cpm_map = {}
    
    bed_file = pd.read_csv(f"data/TCGA/{tcga_barcode}.filter.bed", header=None, sep="\t")
    #cpm calculations: total reads within pre-exRNA filtered bedfile. Note reads that had "chrUn" or low-complexity seqeunces where filered out and not included in cpm normalization.           
    cpm = 1000000/bed_file.shape[0]
    cpm_map[sample] = cpm
    del bed_file
    
    #Counts of loci 
    loci_bed = pd.read_csv(f, header=None, sep="\t")
    loci_features = loci_bed[9] #Locus annotation as a feature.
    loci_bp_overlaps = loci_bed[12] #Number of bp overlaps between read and locus annotation.
    read_ids = loci_bed[3] #Query ID from original bamfile for each read.
    
    #First create read_id map
    read_id_locus_map = {}
    for i in range(len(read_ids)):
        _id = read_ids.iloc[i]
        num_bp_match = loci_bp_overlaps.iloc[i]
        locus = loci_features.iloc[i]
        
        if _id in read_id_locus_map: #Indicates multiple hits/overlaps for one read.
            if num_bp_match > read_id_locus_map[_id]["bp"]: #Update locus count if more number of bp matched. This is our simple binning procedure.
                read_id_locus_map[_id] = {"locus":locus, "bp":num_bp_match}
        else:
            read_id_locus_map[_id] = {"locus":locus, "bp":num_bp_match}
            
    assert len(read_id_locus_map) == len(loci_bed[3].unique()) #Ensures we do not overcount reads   
    
    #Counts of loci
    documents = {}
    for _id in read_id_locus_map: 
        locus = read_id_locus_map[_id]["locus"] 
        if locus in documents:
            rna_dict = documents[locus]
            rna_dict["raw_total"] += 1
            sample_data = rna_dict["samples"]
            sample_data["raw_count"] += 1
            sample_data["cpm"] += cpm
        else:
            sample_data = {'sample type': sample, "tissue": tissue, "study": tcga_barcode, "abb":abb,
                            "raw_count": 1, "cpm": cpm}
            documents[locus] = {"locus":locus, "raw_total":1, "samples":sample_data}

    try:
        rnacol.insert_many(documents.values())
    except Exception as e:
        return f"Error: Insertion {e}. File path: {f}"
    myclient.close()
    return f"Complete: {tcga_barcode} Inserted: {len(documents)} documents"

def batch(ls, n):
    for i in range(0, len(ls), n):
        yield ls[i:i+n]

start = time.perf_counter()       
total_inputs = [f.path for f in os.scandir("data/TCGA") if f.name.endswith(".loci.bed")]


#Log to keep track of progress history
insert_output = open("log/insert_db.log", 'w')
error_output = open ("log/insert_db_error.log", "w")
num_batch = 1
print("Start insert of ", len(total_inputs), "files")

#Parallelize batches
for inputs in batch(total_inputs, 25):
    text = Parallel(n_jobs=25)(delayed(processInput)(i) for i in inputs)
    for result in text:
        if "Error" in result:
            insert_err += 1 #Count number of files that had errors.
            print(result, file=error_output)
        else:
            print(result, file=insert_output)
    batchtime = time.perf_counter()
    timediff = (batchtime - start)/60
    print("Complete batch:", num_batch, "Time Elasped:", timediff)
    start = batchtime
    num_batch += 1

insert_output.close()
error_output.close()
if insert_err > 0 or len(total_inputs) == 0:
    print(f"Unsuccessful run: {insert_err} insert errors and {len(total_inputs)} file inputs. Check log files.")
else:
    print(f"Complete {len(total_inputs)} bamfiles.")
