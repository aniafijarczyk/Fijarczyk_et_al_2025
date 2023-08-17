import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import glob

#subset = '.'
subset = sys.argv[1]

if subset == '.':
    subset_id = '_'
else:
    subset_id = subset

straits = ["assembly_length","genes","effectors","assembly_noreps",
           "N_introns","tRNA","pseudo_tRNA","frac_repeats",
           "LEN_introns","propGC","LEN_intergenic","F_intron_genes",
           "LEN_exons"]


ltraits = ["Genome","Genes","Effectors","Genome w/o repeats",
           "Introns","tRNA","Pseudo tRNA","Repeats",
           "Intron length","GC","Intergenic length",
           "Genes with introns","Exon length"]


def readLog(filename):
    fh = open(filename, "r")
    linie = fh.readlines()
    k = [ele.split('\t') for ele in linie]

    check = 0
    n = []
    for ele in k:
        #print(check)
        if ele[0] == 'Iteration':
            check+=1
        if check == 1:
            n.append(ele)
    dN = pd.DataFrame(n[1:], columns = n[0])
    return(dN)


def analyzeLogs(all_traits_names, subset):

    output = {}

    B = []
    G = []
    L = []
    for strait in all_traits_names:
        print(strait)
        
        # Read file
        fname = subset + "/model_" + strait + "_B/bayestraits_" + strait + "_run1.txt.Log.txt"
        dN = readLog(fname)
        dN['trait'] = strait
        L.append(dN)
        
        n_iterations = dN.shape[0]
        print(n_iterations)
        
        # Get 10 most common model strings
        dS = dN[['Model string', 'Dep / InDep']].value_counts().reset_index()
        dS['Model frequency'] = dS['count'].apply(lambda x: x/n_iterations)
        dS = dS.iloc[0:10,:]
        dS['trait'] = strait
        B.append(dS)
        
        # Get losses and gains
        dN['00_01'] = dN.apply(lambda x: x['q12'] > x['q21'], axis = 1)
        dN['00=01'] = dN.apply(lambda x: x['q12'] == x['q21'], axis = 1)
        dN['10_11'] = dN.apply(lambda x: x['q34'] > x['q43'], axis = 1)
        dN['10=11'] = dN.apply(lambda x: x['q34'] == x['q43'], axis = 1)
        dN['01_11'] = dN.apply(lambda x: x['q24'] > x['q42'], axis = 1)
        dN['01=11'] = dN.apply(lambda x: x['q24'] == x['q42'], axis = 1)
        dN['00_10'] = dN.apply(lambda x: x['q13'] > x['q31'], axis = 1)
        dN['00=10'] = dN.apply(lambda x: x['q13'] == x['q31'], axis = 1)
        
        patho_gain_when_trait_small = dN['00_01'].sum() / n_iterations # q12 > q21
        nochange_when_trait_small = dN['00=01'].sum() / n_iterations # # q12 = q21
        patho_loss_when_trait_small = 1 - (patho_gain_when_trait_small + nochange_when_trait_small)
    
        patho_gain_when_trait_big = dN['10_11'].sum() / n_iterations # q34 > q43
        nochange_when_trait_big = dN['10=11'].sum() / n_iterations # q34 = q43
        patho_loss_when_trait_big = 1 - (patho_gain_when_trait_big + nochange_when_trait_big) # q34 < q43
    
        trait_gain_when_patho = dN['01_11'].sum() / n_iterations # q24 > q42
        nochange_when_patho = dN['01=11'].sum() / n_iterations # q24 = q42
        trait_loss_when_patho = 1 - (trait_gain_when_patho + nochange_when_patho) # q24 < q42
    
        trait_gain_when_nonpatho = dN['00_10'].sum() / n_iterations # q13 > q31
        nochange_when_nonpatho = dN['00=10'].sum() / n_iterations # q13 = q31
        trait_loss_when_nonpatho = 1 - (trait_gain_when_nonpatho + nochange_when_nonpatho) # q13 < q31
    
    
        tab = pd.DataFrame({"gains" : [patho_gain_when_trait_small,
                                 patho_gain_when_trait_big,
                                 trait_gain_when_patho,
                                 trait_gain_when_nonpatho],
                     "no_changes" : [nochange_when_trait_small,
                                    nochange_when_trait_big,
                                    nochange_when_patho,
                                    nochange_when_nonpatho],
                     'losses': [patho_loss_when_trait_small,
                               patho_loss_when_trait_big,
                               trait_loss_when_patho,
                               trait_loss_when_nonpatho]},
                     index = ['00_01','10_11','01_11','00_10'])              
        tab['trait'] = strait
        tab = tab.reset_index()
        G.append(tab)
        
        
    dB = pd.concat(B, ignore_index=True)
    dG = pd.concat(G, ignore_index=True)
    dL = pd.concat(L, ignore_index=True)
    
    output['rates'] = dL
    output['model_strings'] = dB
    output['changes'] = dG

    return(output)

def plotTransitionRates(dataframe, add_name):

    d_rates = dataframe[['trait','q12','q13','q21','q24','q31','q34','q42','q43']].reset_index(drop=True).set_index('trait').stack(
    ).reset_index().rename(columns = {'level_1':'type', 0:'rate'})
    d_rates['rate'] = d_rates['rate'].astype('float')

    g = sns.FacetGrid(d_rates, col="trait", col_wrap=3, sharey=False)
    g.map_dataframe(sns.boxplot, x = "type", y = "rate")
    g.savefig("getTransitions_transition_rates_" + add_name + ".png", dpi=150)


def plotModelStrings(dataframe, add_name):

    g = sns.FacetGrid(dataframe, col="trait", col_wrap=3, sharey=False)
    g.map_dataframe(sns.barplot, x = "Model frequency", y = "Model string", hue = 'Dep / InDep', dodge=False, palette = 'pastel')
    g.savefig("getTransitions_model_strings_" + add_name + ".png", dpi=150)


def plotLossesGains(dataframe, add_name):

    f, axs = plt.subplots(4,4, figsize=(20,20))
    axs = axs.ravel()
    for ax, strait in zip(axs, straits):
        df_ = dataframe.loc[dataframe['trait'] == strait,['index', 'gains', 'no_changes', 'losses']].set_index('index')
        sns.set(font_scale=0.7)
        sns.heatmap(df_, annot=True, linewidth=0.5, cmap="YlGnBu", ax=ax)
        ax.set_title(strait, fontsize=20)
    plt.savefig("getTransitions_losses_gains_" + add_name + ".png")




if __name__ == '__main__':

    # Generate dataframes
    datasets = analyzeLogs(straits, subset)

    # Save transitions
    datasets['changes'].to_csv("getTransitions_" + subset_id + ".tsv", sep='\t', header=True, index=False)

    # Plot transition rates
    plotTransitionRates(datasets['rates'], subset_id)

    # Plot model strings
    plotModelStrings(datasets['model_strings'], subset_id)

    # Plot losses and gains
    plotLossesGains(datasets['changes'], subset_id)


