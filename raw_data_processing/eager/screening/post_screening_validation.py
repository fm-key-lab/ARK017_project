#%% 
import pandas as pd
import numpy as np
import os
import glob

os.chdir('')
ARK_paths_edit_dists=glob.glob('/Users/ad_loris/Documents/key_lab/outputs/eagar/screening/dairy_screening/dairy_screening_after_jan2023/default/editDistance/ARK*')
ARK_paths_ancient_edit_dists=glob.glob('/Users/ad_loris/Documents/key_lab/outputs/eagar/screening/dairy_screening/dairy_screening_after_jan2023/ancient/editDistance/ARK*')

# %% edit distances
edit_dist_outputs=pd.DataFrame()
for edit_dist_path in ARK_paths_edit_dists:
    sample=edit_dist_path.split('/')[-1].split('.unmapped')[0]
    edit_dist_outputs=pd.concat([edit_dist_outputs,pd.DataFrame(pd.read_csv(edit_dist_path,sep='\t',index_col=0).loc['Yersinia_pestis'][[x for x in range (6)]].rename(sample)).transpose()])
    
# %% ancient edit distances
edit_dist_outputs=pd.DataFrame()
for edit_dist_path in ARK_paths_ancient_edit_dists:
    sample=edit_dist_path.split('/')[-1].split('.unmapped')[0]
    edit_dist_outputs=pd.concat([edit_dist_outputs,pd.DataFrame(pd.read_csv(edit_dist_path,sep='\t',index_col=0).loc['Yersinia_pestis'][[x for x in range (6)]].rename(sample)).transpose()])
