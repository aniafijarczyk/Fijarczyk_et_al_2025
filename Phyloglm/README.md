## R script for running logistic regression model with phyloglm from phylolm package. The slope describes strength of association between pathogenicity/non-pathogenicity and the value of a trait.

This code runs logistic regression for 13 genomic features with 10 repetitions of subsampling. Subsampling is explained in "bayesTraits_discrete_models".


**run_phyloglm.R**

The script:  
1) Scales genomic features  
2) Removes species with missing values  
3) Prunes ultrametric tree and scales the distance from root to tip to 1.0  
4) Runs phyloglm using "logistic_MPLE" method  
5) Combines results for all the traits and adjusts p-value  
