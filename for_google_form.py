from helpers import *
import time
import pandas as pd

while True:
    in_chr = "chr4"
    print("Enter the chromosome names required seperared by space [like: chrX]: chr4")
    chromosomes = []
    for i in in_chr.split():
        if i in ["chr" + str(i) for i in range(1,23)]:
            chromosomes.append(i)
        else:
            print(i, " is an invalid name.")
    if len(chromosomes) > 0:
        print("Taking in chromosomes: ", chromosomes)
        while True:
            asked = input("Choose any TF from:  ATAC  CTCF  REST  EP300   :     ")
            if asked in ["ATAC", "CTCF", "REST", "EP300"]:
                break
            else:
                print("Please choose correct TF")
        break

unbound = []
bound = []
all_seqs = []
for chrom in chromosomes:
    with open("data/"+chrom + ".fa") as f:
        chr1 = {i: line.strip() for i, line in enumerate(f)}
        
    chr1_task = pd.read_csv("ref/"+chrom + "_200bp_bins.tsv", sep="\t")

    for i,j,k in tqdm(chr1_task[["start","end",asked]].values, desc = "Reading chromosome:"+chrom, leave = False):
        string = ""
        for line_num in range(int(i/50), int(j/50)):
            string += chr1[line_num]
        all_seqs.append(string)
        if k == "U":
            unbound.append(string)
        else:
            bound.append(string)

for i in range(0,11):
    t1 = time.time()
    print("================         ORDER  ",i,"       =================")
    k_fold_cross_validation(bound, unbound, m = i+1, k = 5, plots = False)
    t2 = time.time()
    print("time taken :", t2-t1)