from helpers import *
import pandas as pd

while True:
    in_chr = input("Enter the chromosome names required seperared by space [like: chrX]: ")
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
        while True:   
            order = int(input("Enter the order of markov model desired:   "))
            if order >= 0:
                break
            else:
                print("Please select a positive integer")
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

m = order +1
unbound_counts = counts_all(unbound,m)
bound_counts = counts_all(bound,m)

model_unbound = {}
for x1 in get_all_kmers(m-1):
    denom = sum([unbound_counts.get(x1+x2,0) for x2 in "ACGT"]) + 1 # pseudo count added to avoid divsion by 0.
    for x2 in "ACGT":
        model_unbound[x1+x2] = unbound_counts.get(x1+x2,0)/denom
model_bound = {}
for x1 in get_all_kmers(m-1):
    denom = sum([bound_counts.get(x1+x2,0) for x2 in "ACGT"]) + 1 # pseudo count added to avoid divsion by 0.
    for x2 in "ACGT":
        model_bound[x1+x2] = bound_counts.get(x1+x2,0)/denom

for seq in all_seqs:
    print(compute_score(seq, model_bound, model_unbound, m))