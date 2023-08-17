import pandas as pd
import numpy as np
import glob
import os.path

traits_dict = {'assembly_length' : 'Genome',
 'assembly_noreps' : 'Genome w/o repeats',
 'effectors' : 'Effectors',
 'frac_repeats' : 'Repeats',
 'F_intron_genes' : 'Genes with introns',
 'genes' : 'Genes',
 'LEN_exons' : 'Exon length',
 'LEN_intergenic' : 'Intergenic length',
 'LEN_introns' : 'Intron length',
 'N_introns' : 'Introns',
 'propGC' : 'GC',
 'pseudo_tRNA' : 'Pseudo tRNA',
 'tRNA' :' tRNA'}



#inputs = glob.glob("bayestraits_*")
#traits = [ele.replace("bayestraits_","").replace(".txt","") for ele in inputs]
traits = list(traits_dict.keys())


def extractRuns(model, list_traits, list_runs, model_name):
    R = []
    for trait in list_traits:
        for run in list_runs:
            filename = 'model_' + trait + '_' + model + '/bayestraits_' + trait + '_' + run + '.txt.Stones.txt'
            if os.path.isfile(filename):
                fh = open(filename, 'r')
                linie = fh.readlines()
                k = float([ele.split() for ele in linie][-1][3])
                
                r = trait, run, k,
                R.append(r)
            else:
                print(model, trait, run, 'Run does not exist !!!')
    dR = pd.DataFrame(R, columns = ['trait', 'run', model_name])
    return(dR)


def calculateLogBF(df, complex_col, simple_col):
    df_ = df.fillna(0)[['trait','run',complex_col, simple_col]].reset_index(drop=True)
    df_['logBF'] = df_.apply(lambda x: 2 * (x[complex_col] - x[simple_col]), axis = 1)
    means = df_.groupby(['trait'])['logBF'].mean().reset_index().rename(columns = {'logBF' : 'mean_logBF'})
    df_ = pd.merge(df_, means, on = 'trait', how = 'left')
    return(df_[['logBF','mean_logBF']].reset_index(drop=True))



if __name__ == '__main__':

    # Extracting log likelihood
    df_independent = extractRuns('A', traits, ['run1', 'run2', 'run3'], 'independent')
    df_dependent = extractRuns('B', traits, ['run1', 'run2', 'run3'], 'dependent')
    df_covarion = extractRuns('C', traits, ['run1', 'run2', 'run3'], 'covarion')
    df_q12 = extractRuns('q12', traits, ['run1', 'run2', 'run3'], 'q12=q21')
    df_q13 = extractRuns('q13', traits, ['run1', 'run2', 'run3'], 'q13=q31')
    df_q24 = extractRuns('q24', traits, ['run1', 'run2', 'run3'], 'q24=q42')
    df_q34 = extractRuns('q34', traits, ['run1', 'run2', 'run3'], 'q34=q43')
    
    # Meging dataframes
    df1 = pd.merge(df_independent, df_dependent, on = ['trait', 'run'], how='left')
    df2 = pd.merge(df1, df_covarion, on = ['trait', 'run'], how='left')
    df3 = pd.merge(df2, df_q12, on = ['trait', 'run'], how='left')
    df4 = pd.merge(df3, df_q13, on = ['trait', 'run'], how='left')
    df5 = pd.merge(df4, df_q24, on = ['trait', 'run'], how='left')
    dfR = pd.merge(df5, df_q34, on = ['trait', 'run'], how='left')

    # Calculating log Bayes Factor
    dfR[['logBF_dependent','mean_logBF_dependent']] = calculateLogBF(dfR, 'dependent', 'independent')
    dfR[['logBF_covarion','mean_logBF_covarion']] = calculateLogBF(dfR, 'covarion', 'dependent')
    dfR[['logBF_q12=q21','mean_logBF_q12=q21']] = calculateLogBF(dfR, 'dependent', 'q12=q21')
    dfR[['logBF_q13=q31','mean_logBF_q13=q31']] = calculateLogBF(dfR, 'dependent', 'q13=q31')
    dfR[['logBF_q24=q42','mean_logBF_q24=q42']] = calculateLogBF(dfR, 'dependent', 'q24=q42')
    dfR[['logBF_q34=q43','mean_logBF_q34=q43']] = calculateLogBF(dfR, 'dependent', 'q34=q43')

    # Adding feature names
    dfR['feature'] = dfR['trait'].apply(lambda x: traits_dict[x])

    # Saving table
    dfR.to_csv('getBayesFactor.out', sep = '\t', header=True, index=False)
    
