import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import os
import re
import PIL
from matplotlib_venn import venn3, venn3_circles
import seaborn
from typing import Any, Dict


def combine_data(results_interaction: pd.DataFrame, results_relation: pd.DataFrame, parameters: Dict[str, Any]):

    gen_set = parameters['gen_set']
    
    res = results_interaction[results_interaction['Ligand'].isin(results_relation['relation']) | results_interaction['Receptor'].isin(results_relation['relation'])]

    df1 = pd.merge(res, results_relation, left_on = ['Ligand'], right_on= 'relation')
    df2 = pd.merge(res, results_relation, left_on = ['Receptor'], right_on= 'relation')

    df_summary = pd.concat([df1, df2])
    
    return df_summary
    

def venn(results_interaction: pd.DataFrame, results_relation: pd.DataFrame, parameters: Dict[str, Any]):

    gen_set = parameters['gen_set']
    
    for gen in gen_set:
        venn3(subsets = (set(results_interaction['Ligand']), set(results_interaction['Receptor']), set(results_relation['relation'][results_relation['gen'] == gen])) , set_labels = ('Ligand genes', 'Receptor genes', 'Relation genes'), alpha = 0.5)
        plt.title('Combine interaction / relation genes to' + str(gen))
        plt.savefig('data/plots/venn_plot_combine_'+ str(gen) + '.png', bbox_inches="tight")
        plt.clf()
        plt.close()
        
        
        
def summary_path(df_summary: pd.DataFrame, parameters: Dict[str, Any]):

    gen_set = parameters['gen_set']
    
    for gen in gen_set:
        tmp =  df_summary[df_summary['gen'] == gen]
        ab = []
        for i in tqdm(tmp['pathways'][tmp['pathways'] != str(['None'])]):
            ab = ab + eval(i)
        values, counts = np.unique(ab, return_counts=True)    
        count = pd.DataFrame({'names':values, 'n':counts})
        count = count[count['names'] != 'None']
        count = count.sort_values('n', ascending=False)
        count = count[0:20]
        fig, ax = plt.subplots()
        seaborn.barplot(count['names'], count['n'])
        plt.xticks(rotation=90)
        plt.xlabel('Pathways')
        plt.ylabel('Number of paths')
        plt.title('Summary of interaction / relation pathways to ' + str(gen))
        plt.savefig('data/plots/path_summary_to_'+ str(gen) + '.png', bbox_inches="tight")
        plt.clf()
        plt.close()
        