# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
# %% 
# Move to relevant directory
pca_ouptput_dir='/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/figures'
os.chdir(pca_ouptput_dir)

# %% load eigen vectors and values
pca_coords=pd.read_csv('/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/metadata_tables/pca_tables/host_validation_pca_with_corrected_labels.csv',header=0,sep=';', decimal=',')
pca_values=pd.read_csv('/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/metadata_tables/pca_tables/PCA_240607_3.eval.txt',header=None)

# %% 
# standardize data for plotting labels
pca_to_plot=pca_coords.query('Plotted == True')

plotting_labels=[]
for i,entry in pca_to_plot.iterrows():
    if entry['Sample label'] == 'ARK017':
        plotting_labels.append(entry[0])
    else:
        plotting_labels.append(entry['Corrected Species'])
pca_to_plot['plotting_labels_pc1_pc2']=plotting_labels
# %%
# plotting pc 1 and pc 2
#plt.rcParams["text.usetex"] = True
font_size=12
c_conversions={'ARK017':'#8e0152', 
                'Ovis aries':      '#e5f5e0', 
                'Ovis orientalis': '#a1d99b',
                'Ovis canadensis':'#31a354',
                'Ovis dalli':'#31a354'}
l_conversions={'Ovis aries':      '$\it{O. aries}$\n(domesticated sheep)', 
                'Ovis orientalis':'$\it{O. orientalis}$',
                'Ovis canadensis':'$\it{O. canadensis}$',
                'Ovis dalli':'$\it{O. dalli}$'}
fig, ax = plt.subplots(facecolor='white', figsize=(4, 4))
for species in ['Ovis canadensis','Ovis dalli', 'Ovis orientalis','Ovis aries', 'ARK017']:
    x_coords=pca_to_plot.query(f'plotting_labels_pc1_pc2 == "{species}"').pc1
    y_coords=pca_to_plot.query(f'plotting_labels_pc1_pc2 == "{species}"').pc2
    if species == 'ARK017':
        plt.scatter(x=x_coords,y=y_coords,c=c_conversions[species],label=species,s=150,marker='x')
        plt.text(np.mean(x_coords), np.mean(y_coords)+0.025 , species, fontsize=font_size,horizontalalignment='left',c=c_conversions[species])
    else:        
        plt.scatter(x=x_coords,y=y_coords,c=c_conversions[species],label=species,edgecolors='black',linewidths=0.5,s=100,alpha=0.7)
        if species=='Ovis canadensis':
            plt.text(np.mean(x_coords) * (1 - 0.05), np.mean(y_coords)+0.01 , l_conversions[species], fontsize=font_size, horizontalalignment='right',verticalalignment='center')
        elif species == 'Ovis dalli':
            plt.text(np.mean(x_coords) * (1 - 0.05), np.mean(y_coords)-0.01 , l_conversions[species], fontsize=font_size, horizontalalignment='right',verticalalignment='center')
        else:
            plt.text(np.max(x_coords)+0.02, np.mean(y_coords)-0.02 ,l_conversions[species], fontsize=font_size, horizontalalignment='left',verticalalignment='center')

pc1_expl=round(pca_values.iloc[0][0]/np.sum(pca_values)[0]*100,2)
pc2_expl=round(pca_values.iloc[1][0]/np.sum(pca_values)[0]*100,2)
plt.xlabel(f'PC 1 ({pc1_expl}%)')
plt.ylabel(f'PC 2 ({pc2_expl}%)')
plt.savefig('ovis_pca_projection_pc1_pc2.svg',bbox_inches='tight')

#c7eae5
#80cdc1
#35978f

# greens
#e5f5e0
#a1d99b
#31a354
# %%
# 3blues
#ece7f2
#a6bddb
#2b8cbe

#species_name_conversions={'OvisOrientalis': 'O. orientalis','AmericanWild':'O. canadensis','Arkaim':'ARK017'}
plotting_labels=[]
for i,entry in pca_to_plot.iterrows():
    if entry['Sample label'] == 'ARK017':
        plotting_labels.append(entry[0])
    elif entry['Corrected Species'] in ['Ovis orientalis', 'Ovis canadensis','Ovis dalli']:
        plotting_labels.append(entry['Corrected Species'])
    else:
        plotting_labels.append(entry['Corrected Region'])
pca_to_plot['plotting_labels_pc2_pc3']=plotting_labels

fig, ax = plt.subplots(facecolor='white', figsize=(4,4))
c_conversions={'Ovis dalli':'#31a354',
               'Ovis canadensis':'#31a354',
                    'Ovis orientalis':'#a1d99b', 
                    'Africa':'#0570b0',
                    'C Europe':'#74a9cf',
                    'SW Europe':'#74a9cf',
                    'N Europe':'#74a9cf',
                    'America':'#bdc9e1',
                    'S Asia':'#f1eef6',
                    'SW Asia':'#f1eef6',
                    'Tibet':'#f1eef6',
                    'Indonesia':'#f1eef6',
                    'ARK017':'#8e0152' }
m_conversion={'Ovis dalli':'o',
              'Ovis canadensis':'o',
              'Ovis orientalis':'o', 
              'Africa':'o',
              'C Europe':'o',
              'SW Europe':'s',
              'N Europe':'d',
              'America':'o',
              'S Asia':'o',
              'SW Asia':'s',
              'Tibet':'d',
              'Indonesia':'v'}
l_conversions={'Ovis dalli':'$\it{O. dalli}$',
               'Ovis canadensis':'$\it{O. canadensis}$',
                    'Ovis orientalis':'$\it{O. orientalis}$', 
                    'Africa':'$\it{O. aries}$ Africa',
                    'C Europe':'$\it{O. aries}$ C Europe',
                    'Indonesia':'$\it{O. aries}$ Indonesia',
                    'N Europe':'$\it{O. aries}$ N Europe',
                    'America':'$\it{O. aries}$ Americas',
                    'S Asia':'$\it{O. aries}$ S Asia',
                    'SW Asia':'$\it{O. aries}$ SW Asia',
                    'SW Europe':'$\it{O. aries}$ SW Europe',
                    'Tibet':'$\it{O. aries}$ Tibet',
                    'ARK017':'ARK017'}
for group in c_conversions:
    x_coords=pca_to_plot.query(f'plotting_labels_pc2_pc3 == "{group}"').pc1
    y_coords=pca_to_plot.query(f'plotting_labels_pc2_pc3 == "{group}"').pc3
    if group == 'ARK017':
        plt.scatter(x=x_coords,y=y_coords,c=c_conversions[group],label=l_conversions[group],s=150,marker='x')
        plt.text(np.mean(x_coords)+0.025, np.mean(y_coords) , species, fontsize=font_size,horizontalalignment='left',c=c_conversions[group])
    else:
        plt.scatter(x=x_coords,y=y_coords,c=c_conversions[group],label=l_conversions[group],edgecolors='black',linewidths=0.5,s=30,marker=m_conversion[group])
        if group=='Ovis canadensis':
            plt.text(np.mean(x_coords) * (1 - 0.025), np.mean(y_coords) , l_conversions[group], fontsize=font_size, horizontalalignment='right',verticalalignment='center')
        elif group=='Ovis dalli':
            plt.text(np.mean(x_coords) * (1 - 0.025), np.mean(y_coords)-0.025 , l_conversions[group], fontsize=font_size, horizontalalignment='right',verticalalignment='center')
        elif group=='Ovis orientalis':
            plt.text(np.mean(x_coords) * (1 + 0.05), np.mean(y_coords)+0.025 , l_conversions[group], fontsize=font_size, horizontalalignment='left',verticalalignment='center')
        #else:
         #   plt.text(np.mean(x_coords)+0.02, np.max(y_coords)-0.02 , l_conversions[group], fontsize=font_size, horizontalalignment='left',verticalalignment='center',label=l_conversions[group])

plt.legend(bbox_to_anchor=(1, 1.025))
pc1_expl=round(pca_values.iloc[0][0]/np.sum(pca_values)[0]*100,2)
pc3_expl=round(pca_values.iloc[2][0]/np.sum(pca_values)[0]*100,2)
plt.xlabel(f'PC 1 ({pc1_expl}%)')
plt.ylabel(f'PC 3 ({pc3_expl}%)')
plt.savefig('ovis_pca_projection_pc1_pc3.svg',bbox_inches='tight')
