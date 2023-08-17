## Scripts to run models of coevolution of two discrete traits: 1) genomic features coded as low (below median = 0) or high (above median = 1) and pathogenicity (yes = 1, no = 0) / insect-association (yes = 1, no = 0).

This code is to run 7 models for 13 genomic features with 10 repetitions of subsampling.

### Specifications of the models

1) `command_Discrete_A.txt`  - model of independent evolution of two traits  
2) `command_Discrete_B.txt`  - model of dependent evolution of two traits  
3) `command_Discrete_C.txt`  - covarion model of evolution of two traits (variable across the tree)  
4) `command_DiscreteDependent_q12.txt` - model of dependent evolution of two traits with equal rates q12 and q21  
5) `command_DiscreteDependent_q13.txt` - model of dependent evolution of two traits with equal rates q13 and q31  
6) `command_DiscreteDependent_q24.txt` - model of dependent evolution of two traits with equal rates q24 and q42  
7) `command_DiscreteDependent_q34.txt` - model of dependent evolution of two traits with equal rates q34 and q43  

### Explanation of transition rates

* q12 = 00 -> 01 (gain)  
* q21 = 01 -> 00 (loss)  
* q13 = 00 -> 10 (gain)  
* q31 = 10 -> 00 (loss)  
* q24 = 01 -> 11 (gain)  
* q42 = 11 -> 01 (loss)  
* q34 = 10 -> 11 (gain)  
* q43 = 11 -> 10 (loss)  

**getSubset.ipynb**

Equal number of pathogenic and non-pathogenic species are subsampled for each genome size bin (every 10Mb between 20 and 70Mb). Species are sampled up to the abundance of the less frequent class.

**plotVariables.R**

Writing input files for bayesTraits, with genomic and lifestyle traites encoded as 0/1.

**convertTree.R**

Subsetting phylogenetic tree for each subsample.

**write_bayesTraits_models.sh**

Script for writing runner script to launch models with independent, dependent and covarion coevolution on computecanada.

**write_bayesTraits_transitions.sh**

Script for writing runner script to launch models with equal transition rates on computecanada.

**getBayesFactor.py**

Script to retrieve log likelihoods from different models and calculate log bayes factors for various model comparisons.

**getTransitions.py**

Script to retrieve a table with transition rates, plot transition rates, frequency of model strings, and counts of gains and losses.

**plotGainsLosses_Heatmap.R**

Script for plotting a pretty heatmap with gains and losses.
