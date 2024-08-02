# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import glob
import gzip
import sys
import subprocess
from sklearn.decomposition import PCA
import re
import scipy.stats as sts
from statsmodels.stats.multitest import multipletests

# %% 
# Functions for data loading and standardization
def parse_breadth_bed_files(breadth_bed_files):
    collector=[]
    for breadth in breadth_bed_files:      
        with gzip.open(breadth) as breadth_data: # define opening procedure gzip.open or open
            current_sample_name=breadth.split('.breadth.gz')[0].split('/')[-1]
            current_data = pd.read_csv(
                breadth_data,sep='\t',names=['chr','start','end','reads','bases_covered','bases_total',current_sample_name],
                dtype={'chr':str,'start':str,'end':str,'reads':float,'bases_covered':float,'bases_total':float,current_sample_name:float})
            current_data['chr_start_end']=current_data['chr']+'_'+current_data['start']+'_'+current_data['end']
            subset_current_data=current_data[['chr_start_end',current_sample_name]].copy()
            subset_current_data.set_index('chr_start_end',inplace=True)
            collector.append(subset_current_data)
    breadth_concat=pd.concat(collector,axis=1)
    return breadth_concat

def parse_depth_bed_files(depth_bed_files):
    collector=[]
    for depth in depth_bed_files:      
        with gzip.open(depth) as depth_data: # define opening procedure gzip.open or open
            current_sample_name=depth.split('.depth.gz')[0].split('/')[-1]
            current_data = pd.read_csv(
                depth_data,sep='\t',names=['chr','start','end',current_sample_name],
                dtype={'chr':str,'start':str,'end':str,'depth':float,current_sample_name:float})
            current_data['chr_start_end']=current_data['chr']+'_'+current_data['start']+'_'+current_data['end']
            subset_current_data=current_data[['chr_start_end',current_sample_name]].copy()
            subset_current_data.set_index('chr_start_end',inplace=True)
            collector.append(subset_current_data)
    depth_concat=pd.concat(collector,axis=1)
    return depth_concat


def get_chr_start_stop_locus_tag_conversions(gff_path):
    conversion_locus_tag={}
    with open(gff_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if not line.startswith('#'):
                chr_name_start_stop_this_line=line.strip().split('\t')[:5]
                parsed=line.strip().split('\t')[8].split(';')
                parsed_dict={x.split('=')[0]:x.split('=')[1] for x in parsed}
                chr_name_start_stop=chr_name_start_stop_this_line[0]+'_'+chr_name_start_stop_this_line[3]+'_'+chr_name_start_stop_this_line[4]
                if parsed_dict['ID'] not in conversion_locus_tag and 'old_locus_tag' in parsed_dict:
                    conversion_locus_tag[chr_name_start_stop]={'ID':parsed_dict['ID'],'old_locus_tag':parsed_dict['old_locus_tag']}
                elif 'Parent' in parsed_dict:
                    if chr_name_start_stop in conversion_locus_tag:
                        conversion_locus_tag[chr_name_start_stop]['accession_id']=parsed_dict['ID'].strip('cds-')
                        if 'gene' in parsed_dict:
                            conversion_locus_tag[chr_name_start_stop]['gene']=parsed_dict['gene']
                        else: 
                            conversion_locus_tag[chr_name_start_stop]['gene']=''
                        if 'product' in parsed_dict:
                            conversion_locus_tag[chr_name_start_stop]['product']=parsed_dict['product']
                        else: 
                            conversion_locus_tag[chr_name_start_stop]['product']=''                
    return conversion_locus_tag

# functions for statistics
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), sts.sem(a)
    h = se * sts.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

ref_genome_folder='/Users/ad_loris/Nextcloud/keylab/reference_genomes/Ypseudotuberculosis_GCF000047365'
gff_file='genome.gff'
chr_start_stop_locus_tag_conversions={}
chr_start_stop_locus_tag_conversions.update(get_chr_start_stop_locus_tag_conversions(f'{ref_genome_folder}/{gff_file}'))
index_tag_conversions=dict()
for x in chr_start_stop_locus_tag_conversions:
    if 'gene' in chr_start_stop_locus_tag_conversions[x]:
        if chr_start_stop_locus_tag_conversions[x]['gene'] != '':
            index_tag_conversions[x] = chr_start_stop_locus_tag_conversions[x]['gene']
        else:
            index_tag_conversions[x] = chr_start_stop_locus_tag_conversions[x]['old_locus_tag']
    else:
        index_tag_conversions[x] = chr_start_stop_locus_tag_conversions[x]['old_locus_tag']


index_tag_conversions_old_locus_tag=dict()
for x in chr_start_stop_locus_tag_conversions:
    index_tag_conversions_old_locus_tag[x] = chr_start_stop_locus_tag_conversions[x]['old_locus_tag']
    if 'gene' in chr_start_stop_locus_tag_conversions[x]:
        if len(chr_start_stop_locus_tag_conversions[x]['gene']) > 0:
            gene_name=chr_start_stop_locus_tag_conversions[x]['gene']
            index_tag_conversions_old_locus_tag[x] = f'({gene_name}) '+  chr_start_stop_locus_tag_conversions[x]['old_locus_tag']       
index_tag_conversions_old_locus_tag_no_additional_info=dict()
for x in chr_start_stop_locus_tag_conversions:
    index_tag_conversions_old_locus_tag_no_additional_info[x] = chr_start_stop_locus_tag_conversions[x]['old_locus_tag']


seen_gene_names={}
for x in index_tag_conversions:
    if index_tag_conversions[x] in seen_gene_names:
        seen_gene_names[index_tag_conversions[x]]+=1
        index_tag_conversions[x]=index_tag_conversions[x]+'_'+str(seen_gene_names[index_tag_conversions[x]])
    else:
        seen_gene_names[index_tag_conversions[x]]=0


# %%
# Create analysis directories
pseudotb_analysis_dir='/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/pangeome/pseudotb_genome_analysis'
os.makedirs(pseudotb_analysis_dir, exist_ok=True)
os.chdir(pseudotb_analysis_dir)

# load order of samples for plotting:
tree_order=np.loadtxt('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/production_run/tree_order.txt', dtype=str)
lnba_sample_names=np.loadtxt('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/production_run/lnba_samples_base_name.txt',dtype=str)
lnba_sample_names=np.array(list(lnba_sample_names)+['HGC068'])

# load all pseudotb samples for heatmap organization
modern_metadata=pd.read_csv('/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/raw_data_processing/eager/modern_accessions.tsv',sep='\t',header=0)
pseudotb_species_samples=modern_metadata.query('species == "Y. pseudotuberculosis"').sample_name

#   LNBA and pre-LNBA metadata
lnba_metadata_path='/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/metadata_tables/lnba_prelnba_metadata.tsv'
lnba_metadata=pd.read_csv(lnba_metadata_path,sep='\t')

# load breadth files
breadth_bed_files=glob.glob('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/pangeome/pseudotb_genome_analysis/bed_files/*breadth*')
depth_bed_files=glob.glob('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/pangeome/pseudotb_genome_analysis/bed_files/*depth*')

breadth=parse_breadth_bed_files(breadth_bed_files)
#depth=parse_depth_bed_files(depth_bed_files)
#tree_order = np.array(list(tree_order) + [x for x in breadth.columns if x not in tree_order])
tree_order = ['pestis_probes'] + [x for x in tree_order if x in breadth.columns]
breadth=breadth[tree_order]
#depth=depth[tree_order]
# rename YAC --> IP32953
breadth=breadth.rename({'YAC': 'IP32953'},axis=1)
breadth=breadth.rename({'pestis_probes': 'Probeset'},axis=1)

non_capture_samples_lnba=np.array(['1343UnTal85', 'RISE509', 'RISE505', 'KunilaII', 'Gyvakarai1','6Post'])
non_capture_samples_all=np.concatenate([non_capture_samples_lnba,['RV2039','Gok2']])
pattern=re.compile('^[0-4]{1}\.')
ancient_samples=breadth.columns[(~np.in1d(breadth.columns,list(filter(pattern.search, list(breadth.columns)))) & ~np.isin(breadth.columns,['IP32953','Probeset']))]
ancient_samples_non_lnba_samples=np.setdiff1d(ancient_samples,lnba_sample_names)


# %% 
# finding which pseudotb ref genes are corresponding to core-genes in the pseudo pangenom
corresponding_line_gff='corresponding_line.gff'
corresponding_line_gff_parsed=pd.read_csv(f'{ref_genome_folder}/{corresponding_line_gff}',sep=' ',names=['chr','ref','gene','start','stop','.','direction','.2','info'])
corresponding_line_gff_parsed['chr_name_start_stop']=[x.chr +'_'+ str(x.start) +'_'+ str(x.stop) for index,x in corresponding_line_gff_parsed.iterrows()]

gff_file_pangenome='test.gff'
gff_file_pangenome_parsed=pd.read_csv(f'{ref_genome_folder}/{gff_file_pangenome}',sep='\t',names=['chr','ref','gene','start','stop','.','direction','.2','info'])
pangenome_id=[x.split(';')[0].split('=')[1] for x in gff_file_pangenome_parsed['info']]
corresponding_line_gff_parsed['pangenome_id']=pangenome_id

pangenome_metadata_dir='/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/pangeome/roary_output_pestis_pseudotb_95'
if not os.path.exists(f'{pangenome_metadata_dir}/header_to_gene_name.txt'):
    subprocess.run([f'grep ">" {pangenome_metadata_dir}/pan_genome_reference.fa | tr -d ">" > {pangenome_metadata_dir}/header_to_gene_name.txt'],shell=True)
name_conversions=pd.read_csv(f'{pangenome_metadata_dir}/header_to_gene_name.txt',sep=' ',header=None)
genome_metadata=pd.read_csv(f'{pangenome_metadata_dir}/pestis_pseudotb_pangenome_genomes_with_header.tsv', sep='\t', header=0)
gene_presence_absence=pd.read_csv(f'{pangenome_metadata_dir}/gene_presence_absence.Rtab', sep='\t',header=0,index_col=0)

prokka_name_to_gene_name={}
with open(f'{pangenome_metadata_dir}/gene_presence_absence.csv') as f:
    for l in f:
        line=l.replace('"','').strip().split(',')
        gene_name=line[0]
        samplesnames_for_gene=np.array(line[14:])
        for s in np.unique(samplesnames_for_gene):
            if 'GCF_000047365' in s:
                prokka_name_to_gene_name[gene_name]=s


"""
chrom_to_gene_name=pd.read_csv(f'{pangenome_metadata_dir}/header_to_gene_name.txt',sep=' ',header=None)
chrom_to_gene_name_conversion={str(x):str(y) for x, y in zip(chrom_to_gene_name[0],chrom_to_gene_name[1])}
"""
ypestis_genomes=genome_metadata.query('species_taxid == 632')
ypseudotb_genomes=genome_metadata.query('species_taxid == 633')

ypestis_genome_gcf = [x.split('.')[0] for x in ypestis_genomes['#assembly_accession']]
ypseudotb_genome_gcf = [x.split('.')[0] for x in ypseudotb_genomes['#assembly_accession']]

# standardize order of reference genomes in presence_absence
gene_presence_absence = gene_presence_absence[ypseudotb_genome_gcf+ypestis_genome_gcf]

# generate analysis sets
core_cutoff = 0.9

present_across_all=gene_presence_absence.index[gene_presence_absence.sum(axis=1)==gene_presence_absence.shape[1]]

# generate present across only ypseudo
present_ypseudo=gene_presence_absence[ypseudotb_genome_gcf].sum(axis=1)
present_ypestis=gene_presence_absence[ypestis_genome_gcf].sum(axis=1)

core_pseudo=(present_ypseudo>len(ypseudotb_genome_gcf)*core_cutoff)
pangenome_pseudo=(present_ypseudo>0)
true_core_pseudo=(present_ypseudo==len(ypseudotb_genome_gcf))
true_core_pestis=(present_ypestis==len(ypestis_genome_gcf))

true_core_pestis_core_pseudo=(true_core_pestis & true_core_pseudo)

def generate_indices_subset_from_pangenome(pangenome_set,corresponding_line_gff_parsed,index_tag_conversions,prokka_name_to_gene_name):
    core_genes_prokka_name=[]
    for x in pangenome_set.index[pangenome_set]:
        if x in prokka_name_to_gene_name:
            core_genes_prokka_name.append(prokka_name_to_gene_name[x])
    breadth_indices_to_subset=corresponding_line_gff_parsed.chr_name_start_stop[corresponding_line_gff_parsed.pangenome_id.isin(core_genes_prokka_name)]
    ypseudo_breadth_indices_to_subset=[index_tag_conversions[x] for x in breadth_indices_to_subset if x in index_tag_conversions]
    return breadth_indices_to_subset,ypseudo_breadth_indices_to_subset

pseudo_core_90,pseudo_core_90_converted_gene_names=generate_indices_subset_from_pangenome(core_pseudo,corresponding_line_gff_parsed,index_tag_conversions,prokka_name_to_gene_name)
pseudo_core_100,pseudo_core_100_converted_gene_names=generate_indices_subset_from_pangenome(true_core_pseudo,corresponding_line_gff_parsed,index_tag_conversions,prokka_name_to_gene_name)

true_core_pestis_tru_core_pseudo_gene_indices,true_core_pestis_tru_core_pseudo_gene_indices_converted_gene_names=generate_indices_subset_from_pangenome(true_core_pestis_core_pseudo,corresponding_line_gff_parsed,index_tag_conversions,prokka_name_to_gene_name)



# %%
# Run QC filtering
#

# Rescale breadth to mean of breadht per sample
rescaled_breadth = breadth / np.mean(breadth,axis=0) 
rescaled_breadth[rescaled_breadth>1] = 1
rescaled_breadth=rescaled_breadth.transpose().iloc[::-1].transpose()

# %% generate collapsed breadth data:
def generate_collapsed_heatmap_data(breadth,include_probes=True):
    branches_1_to_4=(re.compile(r'^[1-4]{1}\.'),'Branches 1-4')
    pe2=(re.compile(r'^0{1}\.PE2'),'0.PE2')
    pe4=(re.compile(r'^0{1}\.PE4'),'0.PE4')
    pe5=(re.compile(r'^0{1}\.PE5'),'0.PE5')
    pe7=(re.compile(r'^0{1}\.PE7'),'0.PE7')
    ant0=(re.compile(r'^0{1}\.ANT'),'0.ANT')
    collapsed_indices=np.full(breadth.columns.shape,False)
    data_to_collapse=[]
    for collapse_pattern,pattern_name in [branches_1_to_4,ant0,pe2,pe4,pe5,pe7]:
        set_to_analyze = np.in1d(breadth.columns, np.array(list(filter(collapse_pattern.search, list(breadth.columns)))))
        data_to_collapse.append(pd.DataFrame({pattern_name:np.max(breadth[breadth.columns[set_to_analyze]],axis=1)}))
        collapsed_indices[set_to_analyze] = True

    ancient_indices_post_lnba = np.where(np.isin(breadth.columns,ancient_samples_non_lnba_samples[~np.isin(ancient_samples_non_lnba_samples,['RV2039','Gok2','Ypseudo_FDA_ARGOS_665','Probeset','RT5','RT6','I2470'])]))[0]
    first_pandemic_genomes_indices = np.where(np.isin(breadth.columns,breadth.columns[ancient_indices_post_lnba][np.isin(breadth.columns[ancient_indices_post_lnba],['EDI001.A','AE1175'])]))
    second_pandemic_genomes_indices = np.where(np.isin(breadth.columns,breadth.columns[ancient_indices_post_lnba][~np.isin(breadth.columns[ancient_indices_post_lnba],['EDI001.A','AE1175'])]))
    ancient_indices_contemp_lnba = np.where(np.isin(breadth.columns,ancient_samples_non_lnba_samples[np.isin(ancient_samples_non_lnba_samples,['RT5','RT6','I2470'])]))[0]
    Probeset_index = np.where(np.isin(breadth.columns,['Probeset']))
    ip32953_index = np.where(np.isin(breadth.columns,['IP32953']))

    # mark as collapsed
    collapsed_indices[Probeset_index] = True
    collapsed_indices[ip32953_index] = True
    collapsed_indices[ancient_indices_post_lnba] = True
    collapsed_indices[ancient_indices_contemp_lnba] = True

    data_to_collapse = data_to_collapse[:1]  + [pd.DataFrame({'2nd Pandemic':np.max(breadth[breadth.columns[second_pandemic_genomes_indices]],axis=1)})]+ data_to_collapse[1:2] + [pd.DataFrame({'1st Pandemic':np.max(breadth[breadth.columns[first_pandemic_genomes_indices]],axis=1)})]+ data_to_collapse[2:-1] + [pd.DataFrame({'RT5/RT6/I2470':np.max(breadth[breadth.columns[ancient_indices_contemp_lnba]],axis=1)})] +  data_to_collapse[-1:]
    data_to_collapse.append(breadth[breadth.columns[~collapsed_indices]])
    data_to_collapse.append(pd.DataFrame({'IP32953':np.max(breadth[breadth.columns[ip32953_index]],axis=1)}))
    if include_probes:
        data_to_collapse.append(pd.DataFrame({'Probeset':np.max(breadth[breadth.columns[Probeset_index]],axis=1)}))
    concatenated_collapsed_breadth=pd.concat(data_to_collapse,axis=1)
    return concatenated_collapsed_breadth

# %% 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# PRODUCTION CODE: FOR GENERATING FIGURES

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# %% PCA
# PCA
# PCA

def generate_pca_big_bang_collapsed(subset_breadth,cutoff,title,filename=''):
    pca=PCA()
    if cutoff != None:
        subset_breadth = subset_breadth > cutoff
    transformed_data=pca.fit(subset_breadth)

    index_rt5_rt6_i2470 = np.where(np.isin(subset_breadth.columns, ['RT5','RT6','I2470']))[0]
    indices_pre_lnba = np.where(np.isin(subset_breadth.columns, ['RV2039','Gok2']))[0]
    indices_lnba= np.where(np.isin(subset_breadth.columns, lnba_sample_names))[0]
    indices_ark= np.where(np.isin(subset_breadth.columns, ['ARK017']))[0]

    indices_non_lnba_ancient = np.where(np.isin(subset_breadth.columns, ancient_samples_non_lnba_samples))[0]
    indices_1st_2nd=np.setdiff1d(indices_non_lnba_ancient,np.concatenate([indices_pre_lnba,index_rt5_rt6_i2470]))
    index_yac = np.where(np.isin(subset_breadth.columns, ['IP32953','Ypseudo_FDA_ARGOS_665']))[0]
    
    # Plot the scatter plot
    fig, ax = plt.subplots(facecolor='white',)
    selected_colors = ['#f1eef6','#d0d1e6','#a6bddb','#74a9cf','#2b8cbe','#045a8d']
    ark_color='#c90076'
    #['#ffffe5','#c2e699','#78c679','#238443']
    #['#ffffe5','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32']
    selected_colors.reverse()
#'#e66101','#fdb863','#b2abd2','#5e3c99']
    plt.rcParams["text.usetex"] = True
    ax.scatter(transformed_data.components_[0][index_yac], transformed_data.components_[1][index_yac], color='Black',edgecolors='Black',linewidths=0.75,  label='$\it{Y. pseudotuberculosis}$',marker='X',s=50)
    plt.rcParams["text.usetex"] = False
    ax.scatter(transformed_data.components_[0][indices_pre_lnba], transformed_data.components_[1][indices_pre_lnba], color=selected_colors[0],edgecolors='Black', linewidths=0.75,label='Pre-LNBA',marker='o')
    ax.scatter(transformed_data.components_[0][indices_lnba], transformed_data.components_[1][indices_lnba], color=selected_colors[1], label='LNBA',edgecolors='Black',linewidths=0.75,marker='s')
    ax.scatter(transformed_data.components_[0][indices_ark], transformed_data.components_[1][indices_ark], color=selected_colors[1],edgecolors=ark_color,linewidths=3,label='ARK017',marker='D')
    ax.scatter(transformed_data.components_[0][index_rt5_rt6_i2470], transformed_data.components_[1][index_rt5_rt6_i2470], color=selected_colors[2],edgecolors='Black',linewidths=0.75, label='RT5/RT6/I2470',marker='o')
    ax.scatter(transformed_data.components_[0][indices_1st_2nd], transformed_data.components_[1][indices_1st_2nd], color=selected_colors[3], edgecolors='Black',linewidths=0.75, label='1st/2nd Pandemic',marker='s')
    order_for_entries=[('0.PE','0.PE',4,'Black','o'),('0.ANT','0.ANT',4,'Black','s'),(r'^[1-4]{1}\.','Branches 1-4',4,'Black','D')]
    for pattern_to_compile,basename,color_index,edge_color,marker_style in order_for_entries:
        pattern=re.compile(pattern_to_compile)
        set_to_analyze = np.in1d(subset_breadth.columns, np.array(list(filter(pattern.search, list(subset_breadth.columns)))))
        ax.scatter(transformed_data.components_[0][set_to_analyze], transformed_data.components_[1][set_to_analyze],color=selected_colors[color_index],edgecolors=edge_color,marker=marker_style,linewidths=0.75, label=f'{basename}')
    # Add labels and legend
    ax.set_ylim(-0.35,0.15)
    ax.set_xlabel('PC 1 ('+str(round(pca.explained_variance_[0],2))+'%)')
    ax.set_ylabel('PC 2 ('+str(round(pca.explained_variance_[1],2))+'%)')
    ax.set_title(f'{title}')
    plt.legend(loc='lower right',borderaxespad=0.1)
    plt.savefig(f'{pseudotb_analysis_dir}/{filename}.svg',bbox_inches='tight')

generate_pca_big_bang_collapsed(rescaled_breadth[rescaled_breadth.index.isin(pseudo_core_90)],None,'',filename='')
generate_pca_big_bang_collapsed(rescaled_breadth[rescaled_breadth.index.isin(pseudo_core_90)],None,'',filename='production_pca_ypseudo_core_gene_content_rescaled_breadth')

######################
# %% BARPLOT
# BARPLOT
# BARPLOT
# %%
# additional QC for polarization of samples, analyses sensitive to within-set variability of samples
# : Remove outliers from matrices
true_core_breadth_rescaled=rescaled_breadth[rescaled_breadth.index.isin(true_core_pestis_tru_core_pseudo_gene_indices)]
samples_passing_breadth_outlier_test=true_core_breadth_rescaled.columns[~(np.sum(true_core_breadth_rescaled>0.9) < np.mean(np.sum(true_core_breadth_rescaled>0.9))-np.std(np.sum(true_core_breadth_rescaled>0.9))*2)]

core_genome_outliers_removed_breadth=rescaled_breadth[samples_passing_breadth_outlier_test]


def generate_paired_barplot_with_95ci_production(breadth_to_use,cutoff,title,filename=''):
    # Assuming you have the following data:
    indices_lnba = np.where(np.isin(breadth_to_use.columns, lnba_sample_names))[0]
    indices_modern = np.where(~np.isin(breadth_to_use.columns, np.concatenate((ancient_samples, np.array(['RV2039','Gok2', 'IP32953','Ypseudo_FDA_ARGOS_665','Probeset'])))))[0]
    indices_basal_contemp_lnba = np.where(np.isin(breadth_to_use.columns, np.array(['RT5','RT6','I2470', 'RV2039','Gok2'])))[0]
    indices_1st_2nd_pandemic=np.where(np.isin(breadth_to_use.columns,np.setdiff1d(ancient_samples,np.concatenate([lnba_sample_names,['RT5','RT6','I2470', 'RV2039','Gok2']]))))

    # Calculate the mean values for each group
    modern_core_pseudo_mean = np.sum(breadth_to_use[breadth_to_use.columns[indices_modern]] > cutoff)
    lnba_core_pseudo_mean = np.sum(breadth_to_use[breadth_to_use.columns[indices_lnba]] > cutoff)
    non_lnba_ancient_core_pseudo_mean =  np.sum(breadth_to_use[breadth_to_use.columns[indices_basal_contemp_lnba]] > cutoff)
    first_second_core_pseudo_mean= np.sum(breadth_to_use[breadth_to_use.columns[indices_1st_2nd_pandemic]] > cutoff)

    # min error,max error
    modern_core_pseudo_mean,modern_core_pseudo_mean_low,modern_core_pseudo_mean_high=mean_confidence_interval(modern_core_pseudo_mean)
    lnba_core_pseudo_mean,lnba_core_pseudo_mean_low,lnba_core_pseudo_mean_high=mean_confidence_interval(lnba_core_pseudo_mean)
    non_lnba_ancient_core_pseudo_mean,non_lnba_ancient_core_pseudo_mean_low,non_lnba_ancient_core_pseudo_mean_high=mean_confidence_interval(non_lnba_ancient_core_pseudo_mean)
    first_second_core_pseudo_mean,first_second_core_pseudo_mean_low,first_second_core_pseudo_mean_high=mean_confidence_interval(first_second_core_pseudo_mean)

    lower_bounds=np.min([modern_core_pseudo_mean_low,lnba_core_pseudo_mean_low,non_lnba_ancient_core_pseudo_mean_low,first_second_core_pseudo_mean_low])
    upper_bounds=np.max([modern_core_pseudo_mean_high,lnba_core_pseudo_mean_high,non_lnba_ancient_core_pseudo_mean_high,first_second_core_pseudo_mean_high])

    # run statistics
    # lnba vs modern
    p_val_lnba_vs_modern=sts.mannwhitneyu(np.sum(breadth_to_use[breadth_to_use.columns[indices_lnba]] > cutoff),np.sum(breadth_to_use[breadth_to_use.columns[indices_modern]]> cutoff))
    print('lnba vs modern',p_val_lnba_vs_modern)
    # lnba vs early
    p_val_lnba_vs_early=sts.mannwhitneyu(np.sum(breadth_to_use[breadth_to_use.columns[indices_lnba]] > cutoff),np.sum(breadth_to_use[breadth_to_use.columns[indices_basal_contemp_lnba]] > cutoff))
    print('lnba vs early',p_val_lnba_vs_early)
    # lnba vs late
    p_val_lnba_vs_late=sts.mannwhitneyu(np.sum(breadth_to_use[breadth_to_use.columns[indices_lnba]] > cutoff),np.sum(breadth_to_use[breadth_to_use.columns[indices_1st_2nd_pandemic]] > cutoff))
    print('lnba vs late',p_val_lnba_vs_late)

    # Set the x positions for the bars
    labels=['Branch 0-4','LNBA', 'Pre-Histoic', 'Historic']
    x = np.arange(len(labels))

    # Set the width of the bars

    # Create the figure and axes
    fig, ax = plt.subplots(facecolor='white',figsize=(5,4))

    # Plot the grouped bars
    ax.bar(x, [modern_core_pseudo_mean,lnba_core_pseudo_mean,non_lnba_ancient_core_pseudo_mean,first_second_core_pseudo_mean], 
           color=['#d0d1e6','#2b8cbe', '#045a8d','#a6bddb'],
           yerr=[modern_core_pseudo_mean_high-modern_core_pseudo_mean,
                lnba_core_pseudo_mean_high-lnba_core_pseudo_mean,
                non_lnba_ancient_core_pseudo_mean_high-non_lnba_ancient_core_pseudo_mean,
                first_second_core_pseudo_mean_high-first_second_core_pseudo_mean],
            capsize=10)

    # fdr correct and covert pvalues
    lnba_comparisons_pvalues=[p_val_lnba_vs_modern.pvalue,p_val_lnba_vs_early.pvalue,p_val_lnba_vs_late.pvalue]
    lnba_comparisons_pvalues_corrected=multipletests(lnba_comparisons_pvalues,0.05,'fdr_by')[1]
    print(lnba_comparisons_pvalues)
    def convert_pval(pval):
        if pval < 0.0005:
            converted_pval='***'
        elif pval < 0.005:
            converted_pval='**'
        elif pval < 0.05:
            converted_pval='*'
        else:
            converted_pval= 'ns'
        return converted_pval
            
    p_val_a_to_b=convert_pval(lnba_comparisons_pvalues[0])
    p_val_b_to_c=convert_pval(lnba_comparisons_pvalues[1])
    p_val_b_to_d=convert_pval(lnba_comparisons_pvalues[2])
    #p_val_a_to_b=np.format_float_scientific(p_val_lnba_vs_modern.pvalue,3)
    #p_val_b_to_c=np.format_float_scientific(p_val_lnba_vs_early.pvalue,3)
    #p_val_b_to_d=np.format_float_scientific(p_val_lnba_vs_late.pvalue,3)

    ax.set_ylim(lower_bounds-.01*lower_bounds,upper_bounds)
    bounds=upper_bounds-lower_bounds-.01*lower_bounds
    print(bounds)
    x1,x2,x3,x4=0,1,2,3
    y=lower_bounds+bounds*1.75
    h=bounds*.1
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.25, c='black')
    ax.text((x1+x2)*0.5,y+h, f'{p_val_a_to_b}', ha='center', va='bottom', color='black')
    ax.plot([ x2, x2, x3,x3], [y+bounds*.2, y+bounds*.2+h, y+bounds*.2+h, y+bounds*.2], lw=1.25, c='black')
    ax.text((x3+x2)*0.5,y+bounds*.2+h, f'{p_val_b_to_c}', ha='center', va='bottom', color='black')
    ax.plot([x2, x2, x4,x4], [y+bounds*.45, y+bounds*.45+h, y+bounds*.45+h, y+bounds*.45], lw=1.25, c='black')
    ax.text((x4+x2)*0.5,y+bounds*.45+h, f'{p_val_b_to_d}', ha='center', va='bottom', color='black')



    # Add labels, title, and legend
    ax.set_xlabel('')
    plt.rcParams["text.usetex"] = True
    ax.set_ylabel('$\it{Y. pseudotuberculosis}$\ncore genes retained in $\it{Y. pestis}$')
    plt.rcParams["text.usetex"] = False
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(max(0,abs(lower_bounds-10)),upper_bounds+50)
    # Display the plot
    if len(filename)>0:
        plt.savefig(f'{filename}.svg',bbox_inches='tight')
    else:
        plt.show()
generate_paired_barplot_with_95ci_production(core_genome_outliers_removed_breadth[core_genome_outliers_removed_breadth.index.isin(pseudo_core_90)],0.9,'',filename='production_barplot_ypseudo_core_gene_cutoff90_modern_lnba_early_late')


######################
# %% HEATMAPS
# HEATMAPS
# HEATMAPS

# formatting names of genic regions
ypseudo_genes_breadth_old_tag_names=rescaled_breadth.rename(index=index_tag_conversions_old_locus_tag_no_additional_info)
ypseudo_genes_breadth_old_tag_names_with_extension=rescaled_breadth.rename(index=index_tag_conversions_old_locus_tag)

# getting subsets of samples for IDing genes fitting patterns of interest
ypseudo_genes_breadth=rescaled_breadth.rename(index=index_tag_conversions)
indices_lnba = np.where(np.isin(ypseudo_genes_breadth.columns, lnba_sample_names))[0]
indices_post_lnba = np.where(~np.isin(ypseudo_genes_breadth.columns, np.concatenate((lnba_sample_names, np.array(['RV2039','Gok2','Probeset', 'IP32953','Ypseudo_FDA_ARGOS_665'])))))[0]


indices_lnba_and_pre= np.where(np.isin(rescaled_breadth.columns, list(lnba_sample_names)+['RV2039','Gok2']))[0]

genes_present_few_lnba_few_present_non=(np.sum(rescaled_breadth[rescaled_breadth.columns[indices_lnba_and_pre]]>0.99,axis=1) > 0) & (np.sum(rescaled_breadth[rescaled_breadth.columns[indices_post_lnba]]<0.9,axis=1) > len(indices_post_lnba)*0.75)

not_in_core=[]
in_core=[]
for region in genes_present_few_lnba_few_present_non[genes_present_few_lnba_few_present_non].index:
    if region not in list(pseudo_core_90):
        not_in_core.append(region)
    else: 
        in_core.append(region)

heatmap_genes_present_few_lnba_few_present_non=rescaled_breadth.loc[rescaled_breadth.index[rescaled_breadth.index.isin(in_core)]].rename(index=index_tag_conversions_old_locus_tag)
heatmap_concatenated_collapsed_breadth=generate_collapsed_heatmap_data(heatmap_genes_present_few_lnba_few_present_non,include_probes=False)


fig, ax = plt.subplots(facecolor='white', figsize=(10, 10))
cbar_ax = sns.heatmap(heatmap_concatenated_collapsed_breadth.transpose(), cmap=sns.dark_palette((10, 60, 100), input="husl", as_cmap=True), linewidths=1, linecolor='lightgray', clip_on=False, ax=ax, cbar=True)
#plt.title('Pseudotuberculosis gene differentiation between LNBA and modern lineages',fontsize=25)
plt.xlabel('')
plt.xticks(rotation=90,fontsize=15)
plt.yticks(fontsize=12)

# Adjust the color bar size
cbar = ax.collections[0].colorbar
cbar.ax.set_position([0.75, 0.2, 0.01, 0.5])
cbar.set_ticks([0, 0.5, 1])
cbar.set_ticklabels(['0', '0.5', '1'])
cbar.ax.tick_params(labelsize=15)

cbar.ax.spines['top'].set_visible(True)
cbar.ax.spines['right'].set_visible(True)
cbar.ax.spines['bottom'].set_visible(True)
cbar.ax.spines['left'].set_visible(True)
cbar.ax.spines['top'].set_color('lightgray')
cbar.ax.spines['right'].set_color('lightgray')
cbar.ax.spines['bottom'].set_color('lightgray')
cbar.ax.spines['left'].set_color('lightgray')
cbar.ax.spines['top'].set_linewidth(1)
cbar.ax.spines['right'].set_linewidth(1)
cbar.ax.spines['bottom'].set_linewidth(1)
cbar.ax.spines['left'].set_linewidth(1)

plt.savefig(f'{pseudotb_analysis_dir}/heatmap_present_lnba_pre_collapsed_production.svg',bbox_inches='tight')


# uncollapsed 
indices_lnba_and_pre= np.where(np.isin(rescaled_breadth.columns, list(lnba_sample_names)+['RV2039','Gok2']))[0]

genes_present_few_lnba_few_present_non=(np.sum(rescaled_breadth[rescaled_breadth.columns[indices_lnba_and_pre]]>0.99,axis=1) > 0) & (np.sum(rescaled_breadth[rescaled_breadth.columns[indices_post_lnba]]<0.9,axis=1) > len(indices_post_lnba)*0.75)

not_in_core=[]
in_core=[]
for region in genes_present_few_lnba_few_present_non[genes_present_few_lnba_few_present_non].index:
    if region not in list(pseudo_core_90):
        not_in_core.append(region)
    else: 
        in_core.append(region)

heatmap_genes_present_few_lnba_few_present_non=rescaled_breadth.loc[rescaled_breadth.index[rescaled_breadth.index.isin(in_core)]].rename(index=index_tag_conversions_old_locus_tag)

fig, ax = plt.subplots(facecolor='white', figsize=(10, 40))
cbar_ax = sns.heatmap(heatmap_genes_present_few_lnba_few_present_non.transpose(), cmap=sns.dark_palette((10, 60, 100), input="husl", as_cmap=True), linewidths=1, linecolor='lightgray', clip_on=False, ax=ax, cbar=True)
#plt.title('Pseudotuberculosis gene differentiation between LNBA and modern lineages',fontsize=25)
plt.xlabel('')
plt.xticks(rotation=90,fontsize=12)
plt.yticks(fontsize=10)

# Adjust the color bar size
cbar = ax.collections[0].colorbar
cbar.ax.set_position([0.75, 0.2, 0.01, 0.5])
cbar.set_ticks([0, 0.5, 1])
cbar.set_ticklabels(['0', '0.5', '1'])
cbar.ax.tick_params(labelsize=15)

cbar.ax.spines['top'].set_visible(True)
cbar.ax.spines['right'].set_visible(True)
cbar.ax.spines['bottom'].set_visible(True)
cbar.ax.spines['left'].set_visible(True)
cbar.ax.spines['top'].set_color('lightgray')
cbar.ax.spines['right'].set_color('lightgray')
cbar.ax.spines['bottom'].set_color('lightgray')
cbar.ax.spines['left'].set_color('lightgray')
cbar.ax.spines['top'].set_linewidth(1)
cbar.ax.spines['right'].set_linewidth(1)
cbar.ax.spines['bottom'].set_linewidth(1)
cbar.ax.spines['left'].set_linewidth(1)

plt.savefig(f'{pseudotb_analysis_dir}/heatmap_present_lnba_pre_production.svg',bbox_inches='tight')

# also just outputting as a tsv
heatmap_genes_present_few_lnba_few_present_non.to_csv(f'{pseudotb_analysis_dir}/normalized_breadth_present_lnba_pre_production.tsv',sep='\t')

# 
# califf regions
all_califf_regions=np.where(np.isin(ypseudo_genes_breadth_old_tag_names.index,np.array(['YPTB2793','YPTB1495','YPTB1058','YPTB3368','YPTB0872','YPTB0873','YPTB0874','YPTB0875','YPTB0876','YPTB0877','YPTB0878','YPTB2180','YPTB2181','YPTB2182','YPTB2183','YPTB2193','YPTB2194','YPTB2195','YPTB2196','YPTB2197','YPTB2198','YPTB2199','YPTB2200','YPTB2201','YPTB2205','YPTB2206','YPTB2207','YPTB2490','YPTB2491','YPTB2492','YPTB2493','YPTB2494','YPTB2495','YPTB2496','YPTB2497','YPTB3450','YPTB3451','YPTB3452','YPTB3453','YPTB3454','YPTB3455','YPTB3456','YPTB3457','YPTB3458','YPTB3459','YPTB2539'])))[0]

heatmap_data=ypseudo_genes_breadth_old_tag_names_with_extension.iloc[all_califf_regions]
heatmap_data_concatenated_collapsed_breadth=generate_collapsed_heatmap_data(heatmap_data,False)

fig, ax = plt.subplots(facecolor='white', figsize=(10, 10))
cbar_ax = sns.heatmap(heatmap_data_concatenated_collapsed_breadth.transpose(), cmap=sns.dark_palette((10, 60, 100), input="husl", as_cmap=True), linewidths=1, linecolor='lightgray', clip_on=False, ax=ax, cbar=True)
#plt.title('Pseudotuberculosis gene differentiation between LNBA and modern lineages',fontsize=25)
plt.xlabel('')
plt.xticks(rotation=90,fontsize=12)
plt.yticks(fontsize=12)

# Adjust the color bar size
cbar = ax.collections[0].colorbar
cbar.ax.set_position([0.75, 0.2, 0.01, 0.5])
cbar.set_ticks([0, 0.5, 1])
cbar.set_ticklabels(['0', '0.5', '1'])
cbar.ax.tick_params(labelsize=15)

cbar.ax.spines['top'].set_visible(True)
cbar.ax.spines['right'].set_visible(True)
cbar.ax.spines['bottom'].set_visible(True)
cbar.ax.spines['left'].set_visible(True)
cbar.ax.spines['top'].set_color('lightgray')
cbar.ax.spines['right'].set_color('lightgray')
cbar.ax.spines['bottom'].set_color('lightgray')
cbar.ax.spines['left'].set_color('lightgray')
cbar.ax.spines['top'].set_linewidth(1)
cbar.ax.spines['right'].set_linewidth(1)
cbar.ax.spines['bottom'].set_linewidth(1)
cbar.ax.spines['left'].set_linewidth(1)

plt.savefig(f'{pseudotb_analysis_dir}/heatmap_califf_branches_collapsed_production.svg',bbox_inches='tight')



# califf regions uncollapsed
heatmap_data=ypseudo_genes_breadth_old_tag_names_with_extension.iloc[all_califf_regions]
fig, ax = plt.subplots(facecolor='white', figsize=(10, 40))
cbar_ax = sns.heatmap(heatmap_data.transpose(), cmap=sns.dark_palette((10, 60, 100), input="husl", as_cmap=True), linewidths=1, linecolor='lightgray', clip_on=False, ax=ax, cbar=True)
#plt.title('Pseudotuberculosis gene differentiation between LNBA and modern lineages',fontsize=25)
plt.xlabel('')
plt.xticks(rotation=90,fontsize=12)
plt.yticks(fontsize=12)

# Adjust the color bar size
cbar = ax.collections[0].colorbar
cbar.ax.set_position([0.75, 0.2, 0.01, 0.5])
cbar.set_ticks([0, 0.5, 1])
cbar.set_ticklabels(['0', '0.5', '1'])
cbar.ax.tick_params(labelsize=15)

cbar.ax.spines['top'].set_visible(True)
cbar.ax.spines['right'].set_visible(True)
cbar.ax.spines['bottom'].set_visible(True)
cbar.ax.spines['left'].set_visible(True)
cbar.ax.spines['top'].set_color('lightgray')
cbar.ax.spines['right'].set_color('lightgray')
cbar.ax.spines['bottom'].set_color('lightgray')
cbar.ax.spines['left'].set_color('lightgray')
cbar.ax.spines['top'].set_linewidth(1)
cbar.ax.spines['right'].set_linewidth(1)
cbar.ax.spines['bottom'].set_linewidth(1)
cbar.ax.spines['left'].set_linewidth(1)

plt.savefig(f'{pseudotb_analysis_dir}/heatmap_califf_branches_uncollapsed_production.svg',bbox_inches='tight')
#########################
# %%
# checking correlation between age and estimated gene content
indices_lnba = np.where(np.isin(core_genome_outliers_removed_breadth.columns, lnba_sample_names))[0]
rescaled_breadth_lnba=core_genome_outliers_removed_breadth[core_genome_outliers_removed_breadth.columns[indices_lnba]][core_genome_outliers_removed_breadth.index.isin(pseudo_core_90)]

date=[]
sum_breadth=[]
samples_included=[]
for x,y in zip(rescaled_breadth_lnba.columns,np.sum(rescaled_breadth_lnba)):
    query_result=lnba_metadata.query(f'Sample_Name == "{x}"')
    if len(query_result)>0:
        date.append(query_result.median_bce_date_calibrated.iloc[0])
        sum_breadth.append(y)
        samples_included.append(x)
import statsmodels.api as sm

X,y=np.array(sum_breadth).reshape(-1, 1) ,np.array(date).reshape(-1, 1) 

X2 = sm.add_constant(X)
est = sm.OLS(y, X2)
est2 = est.fit()
print(est2.summary())

#########################
#
# END PSEUDOTUBERCULOSIS INVESTIGATIONS
#
#########################

# %%
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PESTIS VIRULENCE GENES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def get_chr_start_end_conversion_virulence(gff_path):
    conversion_locus_tag={}
    with open(gff_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if not line.startswith('#'):
                chr_name_start_stop_this_line=line.strip().split('\t')
                chr_name_start_stop=chr_name_start_stop_this_line[0]+'_'+chr_name_start_stop_this_line[1]+'_'+chr_name_start_stop_this_line[2]
                conversion_locus_tag[chr_name_start_stop]=chr_name_start_stop_this_line[3]
    return conversion_locus_tag

def generate_collapsed_heatmap_data(breadth,include_probes=True):
    branches_1_to_4=(re.compile(r'^[1-4]{1}\.'),'Branches 1-4')
    pe2=(re.compile(r'^0{1}\.PE2'),'0.PE2')
    pe4=(re.compile(r'^0{1}\.PE4'),'0.PE4')
    pe5=(re.compile(r'^0{1}\.PE5'),'0.PE5')
    pe7=(re.compile(r'^0{1}\.PE7'),'0.PE7')
    ant0=(re.compile(r'^0{1}\.ANT'),'0.ANT')
    collapsed_indices=np.full(breadth.columns.shape,False)
    data_to_collapse=[]
    for collapse_pattern,pattern_name in [branches_1_to_4,ant0,pe2,pe4,pe5,pe7]:
        set_to_analyze = np.in1d(breadth.columns, np.array(list(filter(collapse_pattern.search, list(breadth.columns)))))
        data_to_collapse.append(pd.DataFrame({pattern_name:np.max(breadth[breadth.columns[set_to_analyze]],axis=1)}))
        collapsed_indices[set_to_analyze] = True

    ancient_indices_post_lnba = np.where(np.isin(breadth.columns,ancient_samples_non_lnba_samples[~np.isin(ancient_samples_non_lnba_samples,['RV2039','Gok2','Ypseudo_FDA_ARGOS_665','Probeset','RT5','RT6','I2470'])]))[0]
    first_pandemic_genomes_indices = np.where(np.isin(breadth.columns,breadth.columns[ancient_indices_post_lnba][np.isin(breadth.columns[ancient_indices_post_lnba],['EDI001.A','AE1175'])]))
    second_pandemic_genomes_indices = np.where(np.isin(breadth.columns,breadth.columns[ancient_indices_post_lnba][~np.isin(breadth.columns[ancient_indices_post_lnba],['EDI001.A','AE1175'])]))
    ancient_indices_contemp_lnba = np.where(np.isin(breadth.columns,ancient_samples_non_lnba_samples[np.isin(ancient_samples_non_lnba_samples,['RT5','RT6','I2470'])]))[0]
    Probeset_index = np.where(np.isin(breadth.columns,['Probeset']))
    ip32953_index = np.where(np.isin(breadth.columns,['IP32953']))

    # mark as collapsed
    collapsed_indices[Probeset_index] = True
    collapsed_indices[ip32953_index] = True
    collapsed_indices[ancient_indices_post_lnba] = True
    collapsed_indices[ancient_indices_contemp_lnba] = True

    data_to_collapse = data_to_collapse[:1]  + [pd.DataFrame({'2nd Pandemic':np.max(breadth[breadth.columns[second_pandemic_genomes_indices]],axis=1)})]+ data_to_collapse[1:2] + [pd.DataFrame({'1st Pandemic':np.max(breadth[breadth.columns[first_pandemic_genomes_indices]],axis=1)})]+ data_to_collapse[2:-1] + [pd.DataFrame({'RT5/RT6/I2470':np.max(breadth[breadth.columns[ancient_indices_contemp_lnba]],axis=1)})] +  data_to_collapse[-1:]
    lnba_uncollapsed_samples=breadth[breadth.columns[~collapsed_indices]].transpose().iloc[::-1].transpose()
    data_to_collapse.append(lnba_uncollapsed_samples)
    data_to_collapse.append(pd.DataFrame({'IP32953':np.max(breadth[breadth.columns[ip32953_index]],axis=1)}))
    if include_probes:
        data_to_collapse.append(pd.DataFrame({'Probeset':np.max(breadth[breadth.columns[Probeset_index]],axis=1)}))
    concatenated_collapsed_breadth=pd.concat(data_to_collapse,axis=1)
    return concatenated_collapsed_breadth


conversion_file_virulence='/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/pangeome/pestis_virulence_analysis/co92_virulence_genes_for_bedtools_mapping.tsv'
virulence_name_conversion=get_chr_start_end_conversion_virulence(conversion_file_virulence)
# load breadth files
breadth_bed_files=glob.glob('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/pangeome/pestis_virulence_analysis/bed_files/*breadth*')
depth_bed_files=glob.glob('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/pangeome/pestis_virulence_analysis/bed_files/*depth*')

breadth=parse_breadth_bed_files(breadth_bed_files)
depth=parse_depth_bed_files(depth_bed_files)
tree_order=np.loadtxt('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/production_run/tree_order.txt', dtype=str)

tree_order = ['pestis_probes'] + [x for x in tree_order if x in breadth.columns]
breadth=breadth[tree_order]
depth=depth[tree_order]
# rename YAC --> IP32953
breadth=breadth.rename({'YAC': 'IP32953'},axis=1)
breadth=breadth.rename({'pestis_probes': 'Probeset'},axis=1)


# %% QC noramlize breadth data
rescaled_breadth = breadth / np.mean(breadth,axis=0) 
rescaled_breadth[rescaled_breadth>1] = 1
# %%

heatmap_data=generate_collapsed_heatmap_data(rescaled_breadth.rename(index=virulence_name_conversion),False)
fig, ax = plt.subplots(facecolor='white', figsize=(40, 10))
cbar_ax = sns.heatmap(heatmap_data.transpose(), cmap=sns.dark_palette((10, 60, 100), input="husl", as_cmap=True), linewidths=1, linecolor='lightgray', clip_on=False, ax=ax, cbar=True)
#plt.title('Pseudotuberculosis gene differentiation between LNBA and modern lineages',fontsize=25)
plt.xlabel('')
plt.xticks(rotation=90,fontsize=8)
plt.yticks(fontsize=12)

# Adjust the color bar size
cbar = ax.collections[0].colorbar
cbar.ax.set_position([0.75, 0.2, 0.01, 0.5])
cbar.set_ticks([0, 0.5, 1])
cbar.set_ticklabels(['0', '0.5', '1'])
cbar.ax.tick_params(labelsize=15)

cbar.ax.spines['top'].set_visible(True)
cbar.ax.spines['right'].set_visible(True)
cbar.ax.spines['bottom'].set_visible(True)
cbar.ax.spines['left'].set_visible(True)
cbar.ax.spines['top'].set_color('lightgray')
cbar.ax.spines['right'].set_color('lightgray')
cbar.ax.spines['bottom'].set_color('lightgray')
cbar.ax.spines['left'].set_color('lightgray')
cbar.ax.spines['top'].set_linewidth(1)
cbar.ax.spines['right'].set_linewidth(1)
cbar.ax.spines['bottom'].set_linewidth(1)
cbar.ax.spines['left'].set_linewidth(1)
cbar.ax.spines['left'].set_linewidth(1)
plt.savefig('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/pangeome/pestis_virulence_analysis/ypestis_virulence_genes_heatmap_collapsed.svg',bbox_inches='tight')
