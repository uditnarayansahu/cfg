import itertools
from tqdm import tqdm
import random
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

def get_all_kmers(k):
    return ["".join(i) for i in itertools.product("ACGT", repeat=k)]

def counts_all(input_set,m):
    """
    Count all the m-mers present in the input set.
    """
    counts = {}
    for seq in tqdm(input_set, desc = "Building model", leave = False):
        for i in range(200-m):
            char = seq[i:i+m]
            counts[char] = counts.get(char, 0) + 1
    return counts

def compute_score(seq, model1,model2, m):
    """
    Compute score for a given sequence using models of order m.
    """
    score = 0.0
    for i in range(200-m):
        char = seq[i:i+m]
        score += np.log10((model1.get(char,0.01)+0.01)/(0.01+model2.get(char,0.01)))
    return score

def plot_roc_curve(y_true, scores, title="ROC and PR Curves", plots=False):
    """
    y_true : list or np.array of true labels (0/1)
    scores : list or np.array of Markov scores
    """
    fpr, tpr, _ = roc_curve(y_true, scores)
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y_true, scores)
    avg_precision = average_precision_score(y_true, scores)
    
    if plots :
        plt.figure(figsize=(14, 6))
        # ROC Curve
        plt.subplot(1, 2,1)
        plt.plot(fpr, tpr, linewidth=2, label=f"AUC = {roc_auc:.3f}")
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC Curve")
        plt.legend(loc="lower right")
        plt.grid(True)

        # Precision-Recall Curve
        plt.subplot(1,2,2)
        plt.plot(recall, precision, linewidth=2, label=f"AP = {avg_precision:.3f}")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title("Precision-Recall Curve")
        plt.legend(loc="lower left")
        plt.grid(True)
    
        plt.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout to prevent suptitle overlap
        plt.show()

    return roc_auc, avg_precision

def do_everything(unbound_train, bound_train, unbound_test, bound_test, m,plots=False):
    """
    Build model for training sets and then compute scores for testing sets, follwed by ROC curve and return AUC.
    """
    unbound_counts = counts_all(unbound_train,m)
    bound_counts = counts_all(bound_train,m)

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

    unbound_test_scores = [compute_score(seq, model_bound, model_unbound, m) for seq in tqdm(unbound_test, desc="Computing scores for unbound set",leave=False)]
    bound_test_scores = [compute_score(seq, model_bound, model_unbound, m) for seq in tqdm(bound_test, desc="Computing scores for unbound set",leave=False)]

    # Combine scores and labels
    y_true = [1]*len(bound_test_scores) + [0]*len(unbound_test_scores)
    scores = bound_test_scores + unbound_test_scores

    auc_value, pr_value = plot_roc_curve(
        y_true,
        scores,
        plots
    )
    return auc_value, pr_value

def k_fold_cross_validation(bound, unbound, m, k,plots):
    random.shuffle(bound)
    random.shuffle(unbound)
    aucc = 0.0
    prr = 0.0
    for i in range(k):
        print("...")
        unbound_train = unbound[:i*len(unbound)//k] + unbound[(i+1)*len(unbound)//k:]
        bound_train = bound[:i*len(bound)//k] + bound[(i+1)*len(bound)//k:]
        unbound_test = unbound[i*len(unbound)//k:(i+1)*len(unbound)//k]
        bound_test = bound[i*len(bound)//k:(i+1)*len(bound)//k]
        auc, pr = do_everything(unbound_train, bound_train, unbound_test, bound_test, m,plots)
        aucc += auc
        prr += pr
    print(f"-------------- Average AUC :    {aucc/k} --------------")
    print(f"-------------- Average PR  :    {prr/k}  --------------")