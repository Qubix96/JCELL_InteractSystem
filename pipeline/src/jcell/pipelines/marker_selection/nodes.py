import pandas as pd
import numpy as np
from tqdm import tqdm
from typing import Any, Dict




def marker_selector(cell_data: pd.DataFrame):

    cell_data.index = [x.upper() for x in cell_data.index]

    cell_data = cell_data.transpose()
    cell_data['subtype'] = cell_data.index


    #Loop start
    results = pd.DataFrame()
    clusters_df = pd.DataFrame()
    clusters = np.array(cell_data.index)


    for cluster in tqdm(clusters): 
        tmp = pd.DataFrame(cell_data)
        tmp.index = cell_data.index
        tmp['subtype'] = None
        tmp['subtype'][tmp.index == cluster] = 1
        tmp['subtype'][tmp['subtype'] != 1] = 0
        tmp = tmp.groupby(tmp['subtype']).mean()
        tmp['subtype'] = tmp.index 
        tmp2 = pd.DataFrame(tmp.transpose())
        tmp2['MCP'] = (np.array(tmp2[1]))-(np.array(tmp2[0]))
        tmp2 = tmp2.sort_values(by = ['MCP'], ascending = False)
        tmp2 = tmp2[tmp2['MCP'] > 0]
        tmp2['subtype'] = str(cluster)
        tmp2['gene'] = tmp2.index
        
        df_tmp = pd.DataFrame(tmp2[['gene','MCP','subtype']])
       
        
        df_tmp = df_tmp[0:99]
        
       
        results = pd.concat([results, df_tmp])
        
    return results

