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

matrix = pd.read_csv('/GIT_Projekty/JCELL/data/hippo.csv')
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

#dictionary = defaultdict(list)



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



# net = Network(notebook = False)
# net.from_nx(network)
# net.show_buttons(filter_=['physics'])
# net.show('example.html')

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

# df2['img'] = np.nan

icons = {
    "Astrocytes": "img/astrocyte.PNG",
    "Oligodendrocytes": "img/oligodendrocyte.PNG",
    "Glutamatergic": "img/neurons.PNG",
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

fig, ax = plt.subplots(figsize=(15,15))


# Note: the min_source/target_margin kwargs only work with FancyArrowPatch objects.
# Force the use of FancyArrowPatch for edge drawing by setting `arrows=True`,
# but suppress arrowheads with `arrowstyle="-"`
nx.draw_networkx_edges(
    G,
    pos=pos,
    ax = ax,
    arrows=True,
    arrowstyle="-",
    min_source_margin=0,
    min_target_margin=0,
    alpha=0.2
)

nx.draw_networkx_nodes(
    G,
    pos=pos,
    ax = ax,
    alpha=0.3,
    node_color = gene_expression.loc[:,'expression'],
    node_size = 500

   
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
    font_size=6,
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
icon_size = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.025
icon_center = icon_size / 2.0

# Add the respective image to each node
for n in G.nodes:
    xf, yf = tr_figure(pos[n])
    xa, ya = tr_axes((xf, yf))
    a = plt.axes([xa - icon_center, ya - icon_center, icon_size, icon_size])
    a.imshow(G.nodes[n]["image"])
    a.axis("off")

plt.show()