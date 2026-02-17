# Computational Functional Genomics : Project 1st half
## Markov Model for TF bound-unbound classification

### 1. Requirements
* Language: Python 3.6+
* Libraries: `numpy`, `pandas`,`scikit-learn`,`itertools`,`matplotlib`,`tqdm`

### 2. File Structure 
* `helpers.py` - Helper functions 
* `main.py` - Primary file to run model on a set of chromosome(s) to build a model. Does train:test split and runs a Markov model.
* `simplerVersion.py` - Simpler version of code to train a Markov model of order m on the full set of sequences in the given file and print out on the screen the log likelihood score of the same sequences using the trained model.
* `for_google_form.py` - Code with loops to get the outputs of avg AUC and avg area under PRC, after performing 5 fold cross-validation on "Chromosome4" file. 
* `data/` - Folder containing the `.fa` files
* `ref/` -Folder containing the annotaions for TFs.

### 3. How to Run

Run the  scriptS from the terminal using the following command to build and test a markov model.
    Example: python main.py
  then Give the inputs when prompted.
  
Input prompts and syntax asked for files are:
  
    1. Chromosome(s) name- "chrX", example: `chr1`, `chr2`
    
    2. TF name- Write one of `ATAC` `CTCF` `REST` `EP300`
    
    3. Order of markov model- Any positive integer less than 200
    
    4. k for k-fold cross validation- Any natural number (1,2,...)
