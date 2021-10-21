import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import networkx as nx
from pyvis.network import Network
import os
import re
import PIL
import scipy.stats as stats
import h2o
import psutil
from h2o.grid.grid_search import H2OGridSearch
from h2o.estimators.gbm import H2OGradientBoostingEstimator
from h2o.estimators.glm import H2OGeneralizedLinearEstimator
import matplotlib.pyplot as plt
import reactome2py
from reactome2py import analysis, content, utils
from collections import defaultdict
from itertools import chain
from operator import methodcaller
from typing import Any, Dict
import math as mt
import seaborn

def gene_relation(cell_data: pd.DataFrame, parameters: Dict[str, Any]):
    
    gen_set = parameters['gen_set']
    cell_data.index = [x.upper() for x in cell_data.index]
    cell_data = cell_data.transpose()

    mem = psutil.virtual_memory()
    mem = str(round(mem.free/1000000000*0.7)) + 'G'
    h2o.init(max_mem_size = mem, enable_assertions = False)

    results = pd.DataFrame()
    clusters_df = pd.DataFrame()


    for gen in gen_set:
        cell_data['expression'] = np.where(cell_data.loc[:,gen] >= np.mean(cell_data.loc[:,gen]), 'up', 'down')
        
        tmp =  h2o.H2OFrame(cell_data)
        tmp['expression'] = tmp['expression'].asfactor()
            
        t = tmp[0].runif()
            
        train = tmp[ t < 0.7 ]
        valid = tmp[ 0.7 <= t ] 

        #Supervisiored
        
        #Grid-share
        #Gradient Boosting Machine (GBM) 
            
        if (len(train)/2) < 30:
            ntrees_opt = list(range(1, mt.floor(len(train)/2), mt.floor(len(train)/4)))
            
        else:
            ntrees_opt = list(range(1, 15, 3))
            
        max_depth_opt = parameters['max_depth_opt']
        learn_rate_opt = parameters['learn_rate_opt']
        
        if (len(train)/2) < 20:
            min_rows = list(range(2, mt.floor(len(train)/2)-1, mt.floor(len(train)/4)))
            
        else:
            min_rows = list(range(2, 10, 3))


            
        hyper_parameters = {"ntrees": ntrees_opt, "max_depth":max_depth_opt, "learn_rate":learn_rate_opt, "min_rows": min_rows}
           
            
        gbm = H2OGridSearch(H2OGradientBoostingEstimator(distribution="multinomial"), hyper_params=hyper_parameters)
            
        gbm.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train, validation_frame= valid)
            
        gbm.get_grid('logloss', decreasing=False)
            
        mod_gbm = gbm.logloss(valid = True)
            
        mod_gbm = pd.DataFrame.from_dict(mod_gbm.items())
            
            
        #Generalized Linear Models (GLM)
            
        alpha = parameters['alpha']
        lambd = parameters['lambda']
        
        hyper_parameters = {'alpha': alpha, 'lambda': lambd}
            
        glm = H2OGridSearch(H2OGeneralizedLinearEstimator(family="binomial"), hyper_params=hyper_parameters)
            
        glm.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train, validation_frame= valid)
            
        glm.get_grid('logloss', decreasing=False)
            
        mod_glm = glm.logloss(valid = True)
            
        mod_glm = pd.DataFrame.from_dict(mod_glm.items())
              
            
        if min(mod_glm[1]) > min(mod_gbm[1]):   
            best_model = gbm.models[0]
            del(glm, mod_glm, gbm, mod_gbm)
            genes = best_model.varimp(use_pandas=True)
            best_model.train(x=tmp.names[1:len(tmp.names)-2], y=tmp.names[len(tmp.names)-1], training_frame=tmp)
        else:
            best_model = glm.models[0]
            del(glm, mod_glm, gbm, mod_gbm)
            genes = best_model.varimp(use_pandas=True)
            best_model.train(x=tmp.names[1:len(tmp.names)-2], y=tmp.names[len(tmp.names)-1], training_frame=tmp)
                
            
            
        genes = genes[genes['scaled_importance'] != 0]

           
            
            #Interaction selecting - p-val and FC
            
        tmp3 = pd.DataFrame(tmp.as_data_frame())
        test = []
        for i in tmp3.columns[tmp3.columns != 'expression']:
            f,p = stats.ttest_ind(tmp3[i][tmp3['expression'] == 'up'], tmp3[i][tmp3['expression'] == 'down'])
            test.append(p)
              
        
        
        tmp3 = tmp3.groupby(tmp3['expression']).mean()
        tmp3 = pd.DataFrame(tmp3.transpose())
        tmp3['p-val'] = test
        tmp3['-log10(p-val)'] = -np.log10(tmp3['p-val'])
        tmp3['FC'] = ((np.array(tmp3['up'])+0.01)/(np.array(tmp3['down'])+0.01))-1
        tmp3['log(FC)'] = np.log(tmp3['FC']+0.01)
        tmp3 = tmp3.sort_values(by = ['FC', 'p-val'], ascending = False)
        tmp3['top100'] = 'blue'
        tmp3['top100'][0:100] = 'red'
            
        df_tmp = tmp3[['FC', 'p-val']]
        df_tmp['relation'] = df_tmp.index
        df_tmp['gen'] = gen
        df_tmp = df_tmp[df_tmp['p-val'] < 0.05]
        df_tmp = df_tmp[0:1000]
            
            
            #Vulcano positive plot
        plt.style.use('seaborn-whitegrid')
        plt.scatter(x=tmp3['log(FC)'], y=tmp3['-log10(p-val)'], color = tmp3['top100'])
        plt.xlabel('log(FC+1)')
        plt.ylabel('-log10(p-val)')
        plt.title('Volcano plot - ' + str(gen))
        plt.savefig('data/plots/gen_' + str(gen) + '.png')
        plt.clf()
        plt.close()
            
                         
        results = pd.concat([results, df_tmp])
        
    h2o.shutdown() 
 
    return results



def path_drug(results: pd.DataFrame):
#####
    results['pathways'] = None

    for l in tqdm(range(0, len(results['relation']))):
        ids = ','.join([results['relation'][l]])
        result = analysis.identifiers(ids=ids)
        token = result['summary']['token']
        
        
        url = 'https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + token
        
        token_result = analysis.token(token, species = ['Homo sapiens', 'Mus musculus'], page_size='-1', page='-1', sort_by='ENTITIES_FDR', 
                                      order='ASC', resource='TOTAL', p_value='0.05', include_disease=True, 
                                      min_entities=None, max_entities=None)
        
       
        pathways = token_result['pathways']
        
        path = [p['name'] for p in pathways]
        if len(path) == 0:
            path = ['None']
        results['pathways'][l] = path
       
        
    #

    list_of_list = []
    for l in results['relation']:
        sub = l.split(', ')
        list_of_list.append(sub)


    result_drug_targets = [utils.genelist_drug_target(ids=",".join(d), source="drugcentral") for d in tqdm(list_of_list)]



    drugs = []
    for i in result_drug_targets:
        if i == None:
            drugs.append(i)
        elif i != None:
            c = []
            if len([ele for ele in i['drug'] if isinstance(ele, dict)]) == 0:
                drugs.append(i['drug']['drugName'])
            else:   
               for drug in i['drug']:
                   c.append(drug['drugName'])
            
               drugs.append(c)
            
    results['drug'] = drugs       

    return results
    
    
    
def plot(results: pd.DataFrame, parameters: Dict[str, Any]):

    gen_set = parameters['gen_set']
    for gen in tqdm(gen_set):
        tmp =  results[results['gen'] == gen]
        ab = []
        for i in tqdm(tmp['pathways'][tmp['pathways'] != str(['None'])]):
            ab = ab + eval(i)
        values, counts = np.unique(ab, return_counts=True)      
        count = pd.DataFrame({'names':values, 'n':counts})
        count = count[count['names'] != 'None']
        count = count.sort_values('n', ascending=False)
        count = count[0:20]
        seaborn.barplot(count['names'], count['n'])
        plt.xticks(rotation=90)
        plt.ylabel('Number of paths')
        plt.title('Pathways of genes relate to ' + str(gen))
        plt.savefig('data/plots/path_related_to_'+ str(gen) + '.png',  bbox_inches='tight')
        plt.clf()
        plt.close()