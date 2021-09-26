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



os.chdir('C:\GIT_Projekty\JCELL')

# recepotrs ligands data
data_source = '/GIT_Projekty/JCELL/database/human_lr_pair.txt'
data = pd.read_csv(data_source, sep='\t')

d = data.iloc[:,[1,2]]


# cells data

matrix = pd.read_csv('/GIT_Projekty/JCELL/data/striatum.csv')
matrix.index = [x.upper() for x in matrix.index]


#gene compare expression
gene_expression = pd.DataFrame(matrix.loc['HTT',:])
gene_expression['expression'] = np.where(matrix.loc['HTT',:] >= np.mean(matrix.loc['HTT',:]), 'red', 'blue')


matrix[matrix > 0] = 1


x = np.array(d.iloc[:,[0]])
y = np.array(d.iloc[:,[1]])
x = np.append(x,y)
x = np.unique(x)


matrix = matrix.loc[set.intersection(set(matrix.index), set(x)),:]



#Ligand ~ receptor 

df = {'Cell1':[],'Cell2':[],'Ligand':[],'Receptor':[]}

for col1 in tqdm(matrix.columns): 
    # print(matrix[col1])
    for col2 in matrix.columns: 
        # print(matrix[col2])
        i = 0
        # print(col1 + ' -- ' + col2)
        for rows in d.index:
            if (d.iloc[rows,0] in matrix.index and d.iloc[rows,1] in matrix.index):
                # print('Next')
                # print(d.iloc[rows,0])
                if (matrix.loc[d.iloc[rows,0],col1] == 1 and matrix.loc[d.iloc[rows,1],col2] == 1):
                    i = i + 1
                   # dictionary[col1 + ' -> ' + col2].append([d.iloc[rows,0] + ' -> ' + d.iloc[rows,1]])
                    #df = df.append({'Cell1':col1,'Cell2':col2,'Ligand':d.iloc[rows,0],'Receptor':d.iloc[rows,1]}, ignore_index=True)
                    df['Cell1'].append(col1)
                    df['Cell2'].append(col2)
                    df['Ligand'].append(d.iloc[rows,0])
                    df['Receptor'].append(d.iloc[rows,1])
                   
                    
df = pd.DataFrame.from_dict(df)
df = df.dropna()                  
df.loc[:,'direction'] = df.iloc[:,2] + ' -> ' + df.iloc[:,3]
df['weight'] = df.groupby(['Cell1','Cell2'])['direction'].transform('nunique')
df.loc[:,'connection'] = pd.factorize(df.loc[:,'direction'])[0]


#Analysis of interactions


import scipy.stats as stats


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
 

#MCA instead of PCA

import mca
mca_df = mca.MCA(analysis_dict)
pd.DataFrame(mca_df.fs_r(1))

#Principal Components Analysis (PCA)

# import h2o
# import psutil
# mem = psutil.virtual_memory()
# mem = str(round(mem.free/1000000000*0.7)) + 'G'
# h2o.init(h2o.init(max_mem_size = mem))

# tmp = pd.DataFrame(analysis_dict)


# tmp = h2o.H2OFrame(tmp)

# from h2o.transforms.decomposition import H2OPCA

# pca_decomp = H2OPCA(k = 20, transform="NONE", pca_method="Power")

# pca_decomp.train(x=tmp.names[1:len(tmp.names)-1], training_frame=tmp)

# pred = pca_decomp.predict(tmp)

# # plt.scatter(pred.as_data_frame()['PC1'], pred.as_data_frame()['PC2'])

# pred.as_data_frame()

tmp = pd.DataFrame(pd.DataFrame(mca_df.fs_r(1)))
tmp.index = analysis_dict.index

import numpy as np

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

labels = tmp.index

def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
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

    

    
    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    den = dendrogram(linkage_matrix, labels = labels , **kwargs)
    return linkage_matrix, den

# setting distance_threshold=0 ensures we compute the full tree.
model = AgglomerativeClustering(distance_threshold=0, n_clusters=None)


model = model.fit(tmp)


plt.title('Hierarchical Clustering Dendrogram')

# plot the top three levels of the dendrogram

linkage_matrix, den = plot_dendrogram(model)
plt.xlabel("Cell ~ Cell - interactome")
plt.savefig('output/tree.png', bbox_inches="tight")
plt.clf()
plt.close()

from scipy.cluster.hierarchy import cut_tree

# analysis_dict['groups'] = cut_tree(linkage_matrix, n_clusters=len(np.unique(den['color_list'])), height=None)

#Marker conections
#

from sklearn.cluster import KMeans

# function returns WSS score for k values from 1 to kmax
def calculate_WSS(points, kmax):
  sse = []
  for k in range(1, kmax+1):
    kmeans = KMeans(n_clusters = k).fit(points)
    centroids = kmeans.cluster_centers_
    pred_clusters = kmeans.predict(points)
    curr_sse = 0
    
    # calculate square of Euclidean distance of each point from its cluster center and add to current WSS
    for i in range(len(points)):
      curr_center = centroids[pred_clusters[i]]
      curr_sse += (points[i, 0] - curr_center[0]) ** 2 + (points[i, 1] - curr_center[1]) ** 2
      
    sse.append(curr_sse)
  return sse

k_val = len(np.unique(den['color_list']))
sse = calculate_WSS(np.array(tmp), k_val + 10)
sse = pd.DataFrame(np.array(sse), columns = ['WSS'])
sse['k'] = range(1,len(sse) + 1)

for j in range(0,len(sse)):
    if j + 1 == len(sse) + 1 or j + 2 == len(sse) + 1 or j + 3 == len(sse) + 1:
        k = len(sse)
        break
    elif sse['WSS'][j] < sse['WSS'][j + 1]*1.3 and sse['WSS'][j]  < sse['WSS'][j + 2]*1.3 and sse['WSS'][j] < sse['WSS'][j + 3]*1.3:
        k = sse['k'][j] + 1
        break

# plot(sse['k'],sse['WSS'])

import h2o
import psutil
mem = psutil.virtual_memory()
mem = str(round(mem.free/1000000000*0.7)) + 'G'
h2o.init(h2o.init(max_mem_size = mem))

df1 = h2o.H2OFrame(tmp)

r = df1[0].runif()

train = df1[ r < 0.7 ]
valid = df1[ 0.7 <= r ] 

#Unsupervisored


#K-means


###



from h2o.estimators.kmeans import H2OKMeansEstimator


cluster_estimator = H2OKMeansEstimator(k=int(k), estimate_k=True, standardize=False)

cluster_estimator.train(x= train.names[1:], training_frame=train,  validation_frame= valid) 

pred = cluster_estimator.predict(df1)



#Conncet data

tmp = h2o.H2OFrame(analysis_dict).cbind(pred)
tmp['predict'] = tmp['predict'].asfactor()

#Mca 2d plot

mca = pd.DataFrame(tmp.as_data_frame())

#Loop start
results = pd.DataFrame()
clusters_df = pd.DataFrame()
clusters = np.unique(tmp['predict'].as_data_frame())
analysis_dict['cluster'] = np.array(pred.as_data_frame())


for cluster in tqdm(clusters):
    print(cluster)   
    tmp2 = pd.DataFrame(tmp.as_data_frame())
    tmp2['predict'][tmp2['predict'] == cluster] = 999
    tmp2['predict'][tmp2['predict'] != 999] = 0
    tmp2['predict'][tmp2['predict'] == 999] = 1
    tmp2 =  h2o.H2OFrame(tmp2)
    tmp2['predict']  = tmp2['predict'].asfactor()
    
    t = tmp2[0].runif()
    
    train = tmp2[ t < 0.7 ]
    valid = tmp2[ 0.7 <= t ] 
    
    
    
    
    
    #Supervisiored
    
    #Grid-share
    #Gradient Boosting Machine (GBM) 
    
    
    ntrees_opt = [2, 5, 10]
    
    max_depth_opt = [2, 3, 4]
    
    learn_rate_opt = [0.1, 0.2]
    
    hyper_parameters = {"ntrees": ntrees_opt, "max_depth":max_depth_opt, "learn_rate":learn_rate_opt}
    
    from h2o.grid.grid_search import H2OGridSearch
    
    from h2o.estimators.gbm import H2OGradientBoostingEstimator
    
    
    gbm = H2OGridSearch(H2OGradientBoostingEstimator(distribution="multinomial"), hyper_params=hyper_parameters)
    
    gbm.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train, validation_frame= valid, nfolds=10)
    
    gbm.sort_by('logloss', increasing=True)
    
    mod_gbm = pd.DataFrame(gbm.sort_by('logloss', increasing=True).as_data_frame())
    
    
    
    
    #Generalized Linear Models (GLM)
    
    from h2o.estimators.glm import H2OGeneralizedLinearEstimator
    
    
    hyper_parameters = {'alpha': [0.01],
                    'lambda': [1e-5, 1e-6]}
    
    
    glm = H2OGridSearch(H2OGeneralizedLinearEstimator(family="binomial"), hyper_params=hyper_parameters)
    
    glm.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train, validation_frame= valid, nfolds=10)
    
    glm.sort_by('logloss', increasing=True)
    
    mod_glm = pd.DataFrame(glm.sort_by('logloss', increasing=True).as_data_frame())
    
    
    
    
    
    
    
    if min(mod_glm['logloss']) > min(mod_gbm['logloss']):   
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
    # best_model.varimp_plot(num_of_features = len(genes))
    
    
    #Interaction selecting - p-val and FC
    
    tmp3 = pd.DataFrame(tmp2.as_data_frame())
    test = []
    for i in tmp3.columns[tmp3.columns != 'predict']:
        f,p = stats.ttest_ind(tmp3[i][tmp3['predict'] == 0], tmp3[i][tmp3['predict'] == 1])
        test.append(p)
      
    tmp3 = tmp3.groupby(tmp3['predict']).mean()
    tmp3 = pd.DataFrame(tmp3.transpose())
    tmp3['p-val'] = test
    tmp3['-log10(p-val)'] = -np.log10(tmp3['p-val'])
    tmp3['FC'] = ((np.array(tmp3[1])+1)/(np.array(tmp3[0])+1))-1
    tmp3['log(FC)'] = np.log(tmp3['FC']+1)
    tmp3 = tmp3.sort_values(by = ['FC', 'p-val'], ascending = False)
    tmp3['top100'] = 'blue'
    tmp3['top100'][0:99] = 'red'
    
    df_tmp = tmp3[['FC', 'p-val']]
    df_tmp['interactors'] = df_tmp.index
    df_tmp['cluster'] = cluster
    df_tmp = df_tmp[df_tmp['p-val'] < 0.05]
    df_tmp = df_tmp[df_tmp['FC'] > 0]
    df_tmp = df_tmp[0:99]
    
    
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-whitegrid')
    
    #Vulcano positive plot
    plt.scatter(x=tmp3['log(FC)'], y=tmp3['-log10(p-val)'], color = tmp3['top100'])
    plt.xlabel('log(FC+1)')
    plt.ylabel('-log10(p-val)')
    plt.title('Volcano plot')
    plt.savefig('output/cluster' + str(cluster) + '.png')
    plt.clf()
    plt.close()
    
                 
    # clusters_df = pd.concat([clusters_df, cluster_tmp])
    results = pd.concat([results, df_tmp])
 


import openpyxl

import re
  
names = np.array(results['interactors'])
receptors = []
for name in names:
    name = re.sub('.* ->', '',  name)
    receptors.append(re.sub('^\s', '',  name))
              
ligands = []
for name in names:
    ligands.append(re.sub(' -> .*', '',  name))

results['ligand'] = ligands
results['receptor'] = receptors

#Path analysis

results['receptro_path'] = None
results['ligand_path'] = None
results['mutual_path'] = None

from reactome2py import content, analysis
import pprint 
import webbrowser
import itertools
import os

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


#Cell
import re
  
names = np.array(analysis_dict.index)
first = []
for name in names:
    name = re.sub('.* ->', '',  name)
    first.append(re.sub('^\s', '',  name))
              
second = []
for name in names:
    second.append(re.sub(' -> .*', '',  name))
     

 
cluster_df = pd.DataFrame()
cluster_df['efector_cells'] = first
cluster_df['receptor_cells'] = second
cluster_df['cluster'] = np.array(analysis_dict['cluster'])



results.to_excel('output/clusters_interactors.xlsx')
# clusters_df.to_excel('output/clusters_cells.xlsx')


#Barplot - interaction pathways

import itertools

clusters = np.unique(tmp['predict'].as_data_frame())

for path in ['ligand_path','receptro_path', 'mutual_path']:
    for cluster in tqdm(clusters):
        ab = list(itertools.chain.from_iterable(results[path][results['cluster'] == cluster]))
        values, counts = np.unique(ab, return_counts=True)    
        count = pd.DataFrame({'names':values, 'n':counts})
        count = count[count['names'] != 'None']
        count = count.sort_values('n', ascending=False)
        count = count[0:20]
        import seaborn
        seaborn.barplot(count['names'], count['n'])
        plt.xticks(rotation=90)
        # plt.tight_layout()
        plt.xlabel('Pathways')
        plt.ylabel('Number of paths')
        plt.title('Interactions ' + str(path) + ' ' + str(cluster))
        plt.savefig('output/interactions_' + str(path) + '_'+ str(cluster) + '.png', bbox_inches="tight")
        plt.clf()
        plt.close()


#Venn graph - potrzeba lepiej zintegrowac go dla ligand - receptor
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
from supervenn import supervenn
supervenn(eval(sets), nam, side_plots=False)
plt.ylabel('Clusters')
plt.title('Venn graph')
plt.savefig('output/venn.png', bbox_inches="tight")
plt.clf()
plt.close()




#AutoML - not neede dla selekcji markerów uzyjemy

# from h2o.automl import H2OAutoML
# aml = H2OAutoML(max_models=20, seed=1)
# aml.train(x=train.names[1:len(train.names)-2], y=train.names[len(train.names)-1], training_frame=train)

#Network

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


icons = {
    "Astrocytes": "img/astrocyte.PNG",
    "Oligodendrocytes": "img/oligodendrocyte.PNG",
    "Glutamatergic": "img/neurons.PNG",
    "Cholinergic": "img/neurons.PNG",
    "GABAergic": "img/neurons.PNG",
    "Dopaminergic" : "img/neurons.PNG",
    "OPC" : "img/opc.PNG",
    "Macrophages" : "img/immuno.PNG",   
    "Microglial" : "img/immuno.PNG",
    "Purkinje" : "img/purkinje.PNG",
    "Granule" : "img/granule.PNG",
    "Golgi" : "img/golgi.PNG",
    "GPCs" : "img/gpc.PNG",
    "NSCs" : "img/nsc.PNG",
    "NPCs" : "img/npc.PNG",
    "Pericytes" : "img/pericyte.PNG",
    "Endothelial" : "img/endothelial.PNG",
    "Tanacytes" : "img/endothelial.PNG",
    "Fibroblasts" : "img/endothelial.PNG"




} 

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

plt.savefig('output/network.png', bbox_inches="tight")
plt.clf()
plt.close()
