from helpers import *
import time
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
        break
while True:   
    order = int(input("Enter the order of markov model desired:   "))
    if order >= 0:
        break
    else:
        print("Please select a positive integer")
m = order +1
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

ratio = float(input("Enter fraction to take for test set :  "))

plots_req = input("Plots required (y/n)  : ")
if plots_req == "y" or plots_req == "Y":
    plots = True
else:
    plots = False

t1 = time.time()

random.shuffle(bound)
random.shuffle(unbound)
unbound_train = unbound[int(ratio*len(unbound)):]
bound_train = bound[int(ratio*len(bound)):]
unbound_test = unbound[:int(ratio*len(unbound))]
bound_test = bound[:int(ratio*len(bound))]

auc, pr = do_everything(unbound_train, bound_train, unbound_test, bound_test, m,plots)
print(f"-------------- Area under ROC :    {auc} --------------")
print(f"-------------- Area under PR  :    {pr}  --------------")
t2 = time.time()
print("time taken :", t2-t1)