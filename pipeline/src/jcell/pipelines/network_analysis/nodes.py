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
import mca
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import cut_tree
from sklearn.cluster import KMeans
import h2o
import psutil
from h2o.estimators.kmeans import H2OKMeansEstimator
from h2o.grid.grid_search import H2OGridSearch
from h2o.estimators.gbm import H2OGradientBoostingEstimator
from h2o.estimators.glm import H2OGeneralizedLinearEstimator
import openpyxl
import itertools
from supervenn import supervenn
import logging
from typing import Any, Dict
from reactome2py import content, analysis
import pprint 
import webbrowser
import seaborn
import math as mt



def load_data(interaction_data: pd.DataFrame, cell_data: pd.DataFrame):

    # recepotrs ligands data
    d = interaction_data.iloc[:,[1,2]]
    d[d.iloc[:, [0]].columns[0]] = d[d.iloc[:, [0]].columns[0]].apply(lambda x: x.upper())
    d[d.iloc[:, [1]].columns[0]] = d[d.iloc[:, [1]].columns[0]].apply(lambda x: x.upper())

    # cells data

    cell_data.index = [x.upper() for x in cell_data.index]

    
    cell_data[cell_data > 0] = 1

    
    x = np.array(d.iloc[:,[0]])
    y = np.array(d.iloc[:,[1]])
    x = np.append(x,y)
    x = np.unique(x)
    x = [n.upper() for n in x]
    
    
    binary_data = cell_data.loc[set.intersection(set(cell_data.index), set(x)),:]

    return binary_data, d



#Ligand ~ receptor 

def inter_share (binary_data: pd.DataFrame, d: pd.DataFrame):
    df = {'Cell1':[],'Cell2':[],'Ligand':[],'Receptor':[]}
   

    for col1 in tqdm(binary_data.columns): 
        for col2 in binary_data.columns: 
            for rows in d.index:
                if (d.iloc[rows,int(0)] in binary_data.index and d.iloc[rows,int(1)] in binary_data.index):
                    if (binary_data.loc[d.iloc[rows,int(0)],col1] == int(1) and binary_data.loc[d.iloc[rows,int(1)],col2] == int(1)):
                        df['Cell1'].append(col1)
                        df['Cell2'].append(col2)
                        df['Ligand'].append(d.iloc[rows,int(0)])
                        df['Receptor'].append(d.iloc[rows,int(1)])
                       
                        
    df = pd.DataFrame.from_dict(df)
    print(df)
    df = df.dropna()                  
    df.loc[:,'direction'] = df.iloc[:,2] + ' -> ' + df.iloc[:,3]
    df['weight'] = df.groupby(['Cell1','Cell2'])['direction'].transform('nunique')
    df.loc[:,'connection'] = pd.factorize(df.loc[:,'direction'])[0]
    
    dfr = df[['Cell1', 'Cell2', 'Ligand', 'Receptor']]
    
    return df, dfr


#Analysis of interactions

def inter_analysis(df: pd.DataFrame):

    analysis_df_ligands = df[['Cell1', 'Cell2' ,'direction']]
    analysis_df_ligands = analysis_df_ligands[analysis_df_ligands['Cell1'] != analysis_df_ligands['Cell2']]
    analysis_df_ligands.loc[:,'cci'] = analysis_df_ligands.iloc[:,0] + ' -> ' + analysis_df_ligands.iloc[:,1]
    analysis_df_ligands['n'] = analysis_df_ligands.groupby('cci')['direction'].transform('nunique')
    analysis_df_ligands = analysis_df_ligands[['cci', 'direction', 'n']]




    analysis_dict = {'cci' : [], 'direction' : [], 'n' : []}
    directions = np.unique(analysis_df_ligands['direction'])
    cells = np.unique(analysis_df_ligands['cci'])

    for i in tqdm(cells):
        tmp = analysis_df_ligands[['cci','direction', 'n']].loc[analysis_df_ligands['cci'] == i]
        for d in directions:
            if not d in np.array(tmp['direction']):
                analysis_dict['direction'].append(d)
                analysis_dict['cci'].append(i)
                analysis_dict['n'].append(int(0))
            elif d in np.array(tmp['direction']):
                analysis_dict['direction'].append(d)
                analysis_dict['cci'].append(i)
                analysis_dict['n'].append(int(1))
                
                
    analysis_dict = pd.DataFrame.from_dict(analysis_dict)

    analysis_dict = analysis_dict.pivot(index='direction', columns='cci', values='n')

    analysis_dict = analysis_dict.transpose()
   
    return analysis_dict
 

#MCA instead of PCA

def mca_dim_red(analysis_dict: pd.DataFrame):
    mca_df = mca.MCA(analysis_dict)
    tmp = pd.DataFrame(mca_df.fs_r(1))
    tmp.index = analysis_dict.index
    label = analysis_dict.index
  
    return tmp, label


def plot_dendrogram(label, tmp, **kwargs):
    model = AgglomerativeClustering(distance_threshold=0, n_clusters=None)
    model = model.fit(tmp)
    labels = label
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,counts]).astype(float)

    # Plot the corresponding dendrogram
    den = dendrogram(linkage_matrix, labels = labels, **kwargs)
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel("Cell ~ Cell - interactome")
    plt.savefig('data/plots/tree.png', bbox_inches="tight")
    plt.clf()
    plt.close()
    return den


# analysis_dict['groups'] = cut_tree(linkage_matrix, n_clusters=len(np.unique(den['color_list'])), height=None)

#Marker conections
#


# function returns WSS score for k values from 1 to kmax
def calculate_k(tmp: pd.DataFrame, den):
    points = np.array(tmp)
    kmax = len(np.unique(den['color_list']))
    sse = []
    for k in range(1, kmax+10):
      kmeans = KMeans(n_clusters = k).fit(points)
      centroids = kmeans.cluster_centers_
      pred_clusters = kmeans.predict(points)
      curr_sse = 0
      
      # calculate square of Euclidean distance of each point from its cluster center and add to current WSS
      for i in range(len(points)):
        curr_center = centroids[pred_clusters[i]]
        curr_sse += (points[i, 0] - curr_center[0]) ** 2 + (points[i, 1] - curr_center[1]) ** 2
        
      sse.append(curr_sse)
    
    
    tmp3 = tmp
    print(tmp3)
    print(sse)
    return sse, tmp3
    
    
def k_estimator(sse):
    sse = pd.DataFrame(np.array(sse), columns = ['WSS'])
    sse['k'] = range(1,len(sse) + 1)
    for j in range(0,len(sse)):
        if j + 1 == len(sse) + 1 or j + 2 == len(sse) + 1 or j + 3 == len(sse) + 1:
            k = len(sse)
            break
        elif sse['WSS'][j] < sse['WSS'][j + 1]*1.3 and sse['WSS'][j]  < sse['WSS'][j + 2]*1.3 and sse['WSS'][j] < sse['WSS'][j + 3]*1.3:
            k = sse['k'][j] + 1
            break
    print(k)
    return k
    
    
# plot(sse['k'],sse['WSS'])

def preapre_data_java(tmp3: pd.DataFrame):
    mem = psutil.virtual_memory()
    mem = str(round(mem.free/1000000000*0.7)) + 'G'
    h2o.init(max_mem_size = mem, enable_assertions = False)

    df1 = h2o.H2OFrame(tmp3)

    r = df1[0].runif()
    train = df1[ r < 0.7 ]
    valid = df1[ 0.7 <= r ] 
    train = pd.DataFrame(train.as_data_frame())
    valid = pd.DataFrame(valid.as_data_frame())
    df1 = pd.DataFrame(df1.as_data_frame())
    
    return valid, train, df1


#Unsupervisored
#K-means


def k_mean(k, train: pd.DataFrame, valid: pd.DataFrame, df1: pd.DataFrame, analysis_dict: pd.DataFrame):

    train = h2o.H2OFrame(train)
    valid = h2o.H2OFrame(valid)
    df1 = h2o.H2OFrame(df1)
    
    cluster_estimator = H2OKMeansEstimator(k=int(k), estimate_k=True, standardize=False)

    cluster_estimator.train(x= train.names[1:], training_frame=train,  validation_frame= valid) 

    pred = cluster_estimator.predict(df1)

    #Conncet data
    #tmp4 = h2o.H2OFrame(analysis_dict).cbind(pred)
    #tmp4['predict'] = tmp4['predict'].asfactor()
    #tmp4 = pd.DataFrame(tmp4.as_data_frame())
    pred = pd.DataFrame(pred.as_data_frame())

    return pred



#Loop start

def marker_selection(analysis_dict: pd.DataFrame, pred, parameters: Dict[str, Any]):

    clusters = np.unique(pred)
    pred = h2o.H2OFrame(pred)
    tmp4 = h2o.H2OFrame(analysis_dict).cbind(pred)
    
    results = pd.DataFrame()
    clusters_df = pd.DataFrame()
    
    analysis_dict['cluster'] = np.array(pred.as_data_frame())
    
   
    for cluster in tqdm(clusters):  
        tmp5 = pd.DataFrame(tmp4.as_data_frame())
        tmp5['predict'][tmp5['predict'] == cluster] = 999
        tmp5['predict'][tmp5['predict'] != 999] = 0
        tmp5['predict'][tmp5['predict'] == 999] = 1
        tmp5 =  h2o.H2OFrame(tmp5)
        tmp5['predict']  = tmp5['predict'].asfactor()
        t = tmp5[0].runif()
        
        train = tmp5[ t < 0.7 ]
        valid = tmp5[ 0.7 <= t ] 
        
     
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
            best_model.train(x=tmp5.names[1:len(tmp5.names)-2], y=tmp5.names[len(tmp5.names)-1], training_frame=tmp5)
        else:
            best_model = glm.models[0]
            del(glm, mod_glm, gbm, mod_gbm)
            genes = best_model.varimp(use_pandas=True)
            best_model.train(x=tmp5.names[1:len(tmp5.names)-2], y=tmp5.names[len(tmp5.names)-1], training_frame=tmp5)
        
        
        
        genes = genes[genes['scaled_importance'] != 0]
        # best_model.varimp_plot(num_of_features = len(genes))
        
        #Interaction selecting - p-val and FC
        
        tmp6 = pd.DataFrame(tmp5.as_data_frame())
        test = []
        for i in tmp6.columns[tmp6.columns != 'predict']:
            f,p = stats.ttest_ind(tmp6[i][tmp6['predict'] == 0], tmp6[i][tmp6['predict'] == 1])
            test.append(p)
          
        tmp6 = tmp6.groupby(tmp6['predict']).mean()
        tmp6 = pd.DataFrame(tmp6.transpose())
        tmp6['p-val'] = test
        tmp6['-log10(p-val)'] = -np.log10(tmp6['p-val'])
        tmp6['FC'] = ((np.array(tmp6[1])+1)/(np.array(tmp6[0])+1))-1
        tmp6['log(FC)'] = np.log(tmp6['FC']+1)
        tmp6 = tmp6.sort_values(by = ['FC', 'p-val'], ascending = False)
        tmp6['top100'] = 'blue'
        tmp6['top100'][0:100] = 'red'
        
        df_tmp = tmp6[['FC', 'p-val']]
        df_tmp['interactors'] = df_tmp.index
        df_tmp['cluster'] = cluster
        df_tmp = df_tmp[df_tmp['p-val'] < 0.05]
        df_tmp = df_tmp[df_tmp['FC'] > 0]
        df_tmp = df_tmp[0:100]
        
        
        plt.style.use('seaborn-whitegrid')
        
        #Vulcano positive plot
        plt.scatter(x=tmp6['log(FC)'], y=tmp6['-log10(p-val)'], color = tmp6['top100'])
        plt.xlabel('log(FC+1)')
        plt.ylabel('-log10(p-val)')
        plt.title('Volcano plot')
        plt.savefig('data/plots/cluster' + str(cluster) + '.png')
        plt.clf()
        plt.close()
                     
        # clusters_df = pd.concat([clusters_df, cluster_tmp])
        results = pd.concat([results, df_tmp])
        names = np.array(results['interactors'])
        
        receptors = []
        for name in names:
            name = re.sub('.* ->', '',  name)
            receptors.append(re.sub('^ ', '',  name))
                      
        ligands = []
        for name in names:
            ligands.append(re.sub(' -> .*', '',  name))

        results['ligand'] = ligands
        results['receptor'] = receptors
        
    analysis_clusters = analysis_dict
    
    h2o.shutdown()
        
    return results, analysis_clusters, clusters
 





#Path analysis
def pathways_reactoma(results: pd.DataFrame, tmp3: pd.DataFrame):

    results['receptro_path'] = None
    results['ligand_path'] = None
    results['mutual_path'] = None


    for l in tqdm(range(0, len(results['ligand']))):
        idsr = ','.join([results['receptor'][l]])
        idsl = ','.join([results['ligand'][l]])
        resultr = analysis.identifiers(ids=idsr)
        resultl = analysis.identifiers(ids=idsl)
        tokenr = resultr['summary']['token']
        tokenl = resultl['summary']['token']
        
        
        url = 'https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + tokenr
        
        token_resultr = analysis.token(tokenr, species = ['Homo sapiens', 'Mus musculus'], page_size='-1', page='-1', sort_by='ENTITIES_FDR', 
                                      order='ASC', resource='TOTAL', p_value='0.05', include_disease=True, 
                                      min_entities=None, max_entities=None)
        
        url = 'https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + tokenl
        
        token_resultl = analysis.token(tokenl, species = ['Homo sapiens', 'Mus musculus'], page_size='-1', page='-1', sort_by='ENTITIES_FDR', 
                                      order='ASC', resource='TOTAL', p_value='0.05', include_disease=True, 
                                      min_entities=None, max_entities=None)
        
        pathwaysr = token_resultr['pathways']
        pathwaysl = token_resultl['pathways']
        
        pathways_receptor = [p['name'] for p in pathwaysr]
        if len(pathways_receptor) == 0:
            pathways_receptor = ['None']
        results['receptro_path'][l] = pathways_receptor
        pathways_ligand = [p['name'] for p in pathwaysl]
        if len(pathways_ligand) == 0:
            pathways_ligand = ['None']
        results['ligand_path'][l] = pathways_ligand
        
        
        
        path = []
        for i in pathways_receptor:
            if i in pathways_ligand:
                path.append(i)
        
        if len(path) == 0:
            path = ['None']
        results['mutual_path'][l] = path 
        path = None
        final_results = results
        
    return final_results


#Cell


def pathways_plot(analysis_clusters: pd.DataFrame, clusters, final_results):
    names = np.array(analysis_clusters.index)
    first = []
    for name in names:
        name = re.sub('.* ->', '',  name)
        first.append(re.sub('^ ', '',  name))
                  
    second = []
    for name in names:
        second.append(re.sub(' -> .*', '',  name))
         

     
    cluster_df = pd.DataFrame()
    cluster_df['efector_cells'] = first
    cluster_df['receptor_cells'] = second
    cluster_df['cluster'] = np.array(analysis_clusters['cluster'])


    #Barplot - interaction pathways

    for path in ['ligand_path','receptro_path', 'mutual_path']:
        for cluster in tqdm(clusters):
            ab = []
            for i in tqdm(final_results[path][final_results['cluster'] == cluster]):
                ab = ab + eval(i)
            values, counts = np.unique(ab, return_counts=True)    
            count = pd.DataFrame({'names':values, 'n':counts})
            count = count[count['names'] != 'None']
            count = count.sort_values('n', ascending=False)
            count = count[0:20]
            seaborn.barplot(count['names'], count['n'])
            plt.xticks(rotation=90)
            plt.ylabel('Number of paths')
            plt.title('Interactions ' + str(path) + ' cluster  ' + str(cluster))
            plt.savefig('data/plots/interactions_' + str(path) + '_'+ str(cluster) + '.png', bbox_inches="tight")
            plt.clf()
            plt.close()
    
    return cluster_df
            


#Venn graph
def venn_interactions(clusters, cluster_df):
    sets = str()
    nam = str()
    for clu in clusters:
        if clu == 0:
            sets = sets + "[set(cluster_df['efector_cells'][cluster_df['cluster'] == " + str(clu) + "]),"
            nam = nam + '"Efector_cells - Cluster ' + str(clu) + '",'
            sets = sets + "set(cluster_df['receptor_cells'][cluster_df['cluster'] == " + str(clu) + "]),"
            nam = nam + ' "Receptor cells - Cluster ' + str(clu) + '",'
        elif clu == len(clusters) - 1:
            sets = sets + " set(cluster_df['efector_cells'][cluster_df['cluster'] == " + str(clu) + "]),"
            nam = nam + ' "Efector_cells - Cluster ' + str(clu) + '",'
            sets = sets + " set(cluster_df['receptor_cells'][cluster_df['cluster'] == " + str(clu) + "])]"
            nam = nam + ' "Receptor cells - Cluster ' + str(clu) + '"'
        else:
            sets = sets + " set(cluster_df['efector_cells'][cluster_df['cluster'] == " + str(clu) + "]),"
            nam = nam + '"Efector_cells - Cluster ' + str(clu) + '",'
            sets = sets + " set(cluster_df['receptor_cells'][cluster_df['cluster'] == " + str(clu) + "]),"
            nam = nam + '"Receptor cells - Cluster ' + str(clu) + '",'
     
    nam = eval(nam)
    supervenn(eval(sets), nam, side_plots=False)
    plt.ylabel('Clusters')
    plt.title('Venn graph')
    plt.savefig('data/plots/venn.png', bbox_inches="tight")
    plt.clf()
    plt.close()




#AutoML - not neede dla selekcji markerów uzyjemy

# from h2o.automl import H2OAutoML
# aml = H2OAutoML(max_models=20, seed=1)
# aml.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train)

#Network
def network(df: pd.DataFrame, parameters: Dict[str, Any], cell_data: pd.DataFrame):
    
    gen_set = parameters['gen_set']
    cell_data.index = [x.upper() for x in cell_data.index]
    
    if gen_set == 0:
        icons = parameters['icons']
        nodes = df['Cell1'].drop_duplicates()
        nodes = np.array(nodes)
        names = np.array(nodes)

        n = 0
        for i in names:
            names[n] = re.sub(" .*", "", names[n])
            n+=1 

        df2 = pd.DataFrame()
        df2['cell'] = nodes
        df2['names'] = names


        images = {k: PIL.Image.open(fname) for k, fname in icons.items()}



                    
        df3 = df[['Cell1', 'Cell2', 'weight']].drop_duplicates()
        df3 = df3[df3['Cell1'] != df3['Cell2']]
        df3['weight'] = (df3['weight']-min(df3['weight']))/(max(df3['weight'])-min(df3['weight']))


                 
        G = nx.Graph() 

        for i in range(0,len(df2)):
            G.add_node(df2.iloc[i,0], image=images[df2.iloc[i,1]])
            
        for weight in range(0,len(df3)):
            G.add_edge(df3.iloc[weight, 0], df3.iloc[weight, 1], weight=df3.iloc[weight, 2])
            



        pos = nx.spring_layout(G, seed = 123)

        v = list(pos.values())
        w_range = []
        for i in v: w_range.append(i[0])
        h_range = []
        for i in v: h_range.append(i[1])

        fig, ax = plt.subplots(figsize=(30,30))


        nx.draw_networkx_edges(
            G,
            pos=pos,
            ax = ax,
            arrows=True,
            arrowstyle="-",
            min_source_margin=0,
            min_target_margin=0,
            alpha=0.1
        )


        pos_higher = {}

        for k, v in pos.items():
            if(v[0]>np.median(w_range) and v[1]>np.median(h_range)):
                pos_higher[k] = (v[0]+0.1, v[1]+0.1)
            elif(v[0]>np.median(w_range) and v[1]<np.median(h_range)):
                pos_higher[k] = (v[0]+0.1, v[1]-0.1)
            elif(v[0]<np.median(w_range) and v[1]>np.median(h_range)):
                pos_higher[k] = (v[0]-0.1, v[1]+0.1)
            else:
                pos_higher[k] = (v[0]-0.1, v[1]-0.1)
                
        nx.draw_networkx_labels(G,
            pos=pos_higher,
            font_color='red',
            font_size=10,
            ax = ax,
            font_family = 'Helvetica',
            font_weight="bold"

        )

        ax.set_yticks(np.arange(min(h_range)*1.25, max(h_range)*1.35, 0.1))
        ax.set_xticks(np.arange(min(w_range)*1.3, max(w_range)*1.45, 0.1))

        # Transform from data coordinates (scaled between xlim and ylim) to display coordinates
        tr_figure = ax.transData.transform
        # Transform from display to figure coordinates
        tr_axes = fig.transFigure.inverted().transform

        # Select the size of the image (o dorelative to the X axis)
        icon_size = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.02
        icon_center = icon_size / 2.0

        # Add the respective image to each node
        for n in G.nodes:
            xf, yf = tr_figure(pos[n])
            xa, ya = tr_axes((xf, yf))
            a = plt.axes([xa - icon_center, ya - icon_center, icon_size, icon_size])
            a.imshow(G.nodes[n]["image"])
            a.axis("off")
            
        plt.title('Interactome network')
        plt.savefig('data/plots/network.png', bbox_inches="tight")
        plt.clf()
        plt.close()
    
    else:
        for gen in gen_set:
        
            #gene compare network
            gene_expression = pd.DataFrame(cell_data.loc[gen,:])
            gene_expression['expression'] = np.where(cell_data.loc[gen,:] >= np.mean(cell_data.loc[gen,:]), 'red', 'blue')


            icons = parameters['icons']
            nodes = df['Cell1'].drop_duplicates()
            nodes = np.array(nodes)
            names = np.array(nodes)

            n = 0
            for i in names:
                names[n] = re.sub(" .*", "", names[n])
                n+=1 

            df2 = pd.DataFrame()
            df2['cell'] = nodes
            df2['names'] = names


            images = {k: PIL.Image.open(fname) for k, fname in icons.items()}



                        
            df3 = df[['Cell1', 'Cell2', 'weight']].drop_duplicates()
            df3 = df3[df3['Cell1'] != df3['Cell2']]
            df3['weight'] = (df3['weight']-min(df3['weight']))/(max(df3['weight'])-min(df3['weight']))


                     
            G = nx.Graph() 

            for i in range(0,len(df2)):
                G.add_node(df2.iloc[i,0], image=images[df2.iloc[i,1]])
                
            for weight in range(0,len(df3)):
                G.add_edge(df3.iloc[weight, 0], df3.iloc[weight, 1], weight=df3.iloc[weight, 2])
                



            pos = nx.spring_layout(G, seed = 123)

            v = list(pos.values())
            w_range = []
            for i in v: w_range.append(i[0])
            h_range = []
            for i in v: h_range.append(i[1])

            fig, ax = plt.subplots(figsize=(30,30))


            nx.draw_networkx_edges(
                G,
                pos=pos,
                ax = ax,
                arrows=True,
                arrowstyle="-",
                min_source_margin=0,
                min_target_margin=0,
                alpha=0.1
            )

            nx.draw_networkx_nodes(
                G,
                pos=pos,
                ax = ax,
                alpha=0.2,
                node_color = gene_expression.loc[:,'expression'],
                node_size = 3000

               
            )

            pos_higher = {}

            for k, v in pos.items():
                if(v[0]>np.median(w_range) and v[1]>np.median(h_range)):
                    pos_higher[k] = (v[0]+0.1, v[1]+0.1)
                elif(v[0]>np.median(w_range) and v[1]<np.median(h_range)):
                    pos_higher[k] = (v[0]+0.1, v[1]-0.1)
                elif(v[0]<np.median(w_range) and v[1]>np.median(h_range)):
                    pos_higher[k] = (v[0]-0.1, v[1]+0.1)
                else:
                    pos_higher[k] = (v[0]-0.1, v[1]-0.1)
                    
            nx.draw_networkx_labels(G,
                pos=pos_higher,
                font_color='red',
                font_size=10,
                ax = ax,
                font_family = 'Helvetica',
                font_weight="bold"

            )

            ax.set_yticks(np.arange(min(h_range)*1.25, max(h_range)*1.35, 0.1))
            ax.set_xticks(np.arange(min(w_range)*1.3, max(w_range)*1.45, 0.1))

            # Transform from data coordinates (scaled between xlim and ylim) to display coordinates
            tr_figure = ax.transData.transform
            # Transform from display to figure coordinates
            tr_axes = fig.transFigure.inverted().transform

            # Select the size of the image (o dorelative to the X axis)
            icon_size = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.02
            icon_center = icon_size / 2.0

            # Add the respective image to each node
            for n in G.nodes:
                xf, yf = tr_figure(pos[n])
                xa, ya = tr_axes((xf, yf))
                a = plt.axes([xa - icon_center, ya - icon_center, icon_size, icon_size])
                a.imshow(G.nodes[n]["image"])
                a.axis("off")
            plt.suptitle('Interactome network - ' + gen)
            plt.savefig('data/plots/network_' + gen + '.png', bbox_inches="tight")
            plt.clf()
            plt.close()
