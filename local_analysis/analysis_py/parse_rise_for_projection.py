import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import gzip
import pickle
import sys

sys.path.append("/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/local_analysis/analysis_py")
import analysispy_modules_pestis as apyp

os.chdir('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/production_run/projections')

# load cmt for adding projection
cmtFile_sub_gp = 'candidate_mutation_table_processed_gp.pickle.gz'

with gzip.open('../final_cmt/' + cmtFile_sub_gp, 'rb') as pickle_file:
    pickle_dictionary = pickle.load(pickle_file)
    counts = pickle_dictionary['counts']
    quals = pickle_dictionary['quals']
    coverage_forward_strand = pickle_dictionary['coverage_forward_strand']
    coverage_reverse_strand = pickle_dictionary['coverage_reverse_strand']
    refnti_m = pickle_dictionary['refnti_m']
    p = pickle_dictionary['p']
    refgenome = pickle_dictionary['refgenome']
    sampleNames = pickle_dictionary['sampleNames']
    outgroup_bool = pickle_dictionary['outgroup_bool']
    contig_positions = pickle_dictionary['contig_positions']
    mutantAF = pickle_dictionary['mutantAF']
    maf = pickle_dictionary['maf']
    maNT = pickle_dictionary['maNT']
    minorNT = pickle_dictionary['minorNT']
    minorAF = pickle_dictionary['minorAF']
    calls = pickle_dictionary['calls']
    hasmutation = pickle_dictionary['hasmutation']
    goodpos_raw_cmt = pickle_dictionary['goodpos_raw_cmt']
    goodpos = pickle_dictionary['goodpos_cleaned_cmt']
    analysis_params_output_name = pickle_dictionary['analysis_params_output_name']

# set ref nt folder
ref_genome_folder='/Users/ad_loris/Nextcloud/keylab/reference_genomes/Ypestis_ASM906v1'

# read LNBA sample names
lnba_sample_names=np.loadtxt('../lnba_samples_base_name.txt',dtype=str)
lnba_sample_names_with_rise=np.append(lnba_sample_names,'RISE386')
lnba_sample_names_with_rise_c90=np.append(lnba_sample_names,['RISE386','C90'])

tree_sample_order=np.loadtxt('../tree_order.txt',dtype=str)
lnba_order=[x for x in tree_sample_order if x in lnba_sample_names]
# read in variant chrom-positions 
goodpos_positions=pd.read_csv('rise386_projection/chrom_pos_goodpos.tsv',sep='\t',header=None)
# %% 
# parise information from reference genome
[chrStarts, genomeLength, scafNames] = apyp.genomestats(ref_genome_folder)


# %%
# parse pileup minimally for agreement in SNP calls that overlap variable positions already ID'd
# read in pileup for RISE386
rise386_pileup=pd.read_csv('rise386_projection/variant_pos.pileup',sep='\t',header=None) 

c90_pileup=pd.read_csv('c90_projection/C90_variant_pos.pileup',sep='\t',header=None) 

def get_projection_from_pileup(pileup,goodpos_positions):
    projection=[]
    for index,position in goodpos_positions.iterrows():
        chrom,pos=position[0],position[1]
        index_to_subset=(pileup[0]==chrom) & (pileup[1]==pos)
        if np.sum(index_to_subset)==1:
            pileup_indexed=pileup.loc[(pileup[0]==chrom) & (pileup[1]==pos)].reset_index()
            if '.' in pileup_indexed[4][0] or  ',' in pileup_indexed[4][0]:
                projection.append(pileup_indexed[2][0])
            elif '*' in pileup_indexed[4][0]:
                projection.append('?')
            else:
                site_calls=[x.upper() for x in pileup_indexed[4][0] if x in ['a','t','c','g','A','T','C','G']]
                if np.all(np.array(site_calls)==site_calls[0]):
                    projection.append(site_calls[0])
        else:
            projection.append('?')
    return projection


c90_projection=get_projection_from_pileup(c90_pileup,goodpos_positions)
rise386_projection=get_projection_from_pileup(rise386_pileup,goodpos_positions)

# %% 
# add projection to calls and samplenames
rise386_projection_calls=apyp.nts2idx(np.array(rise386_projection))
c90_projection_calls=apyp.nts2idx(np.array(c90_projection))

calls_projection=np.append(calls,np.array([rise386_projection_calls,c90_projection_calls]).transpose(),axis=1)
sampleNames_projection=np.append(sampleNames,['RISE386','C90'])

# %% 
# get indices of rise, all LNBA
LNBA_clade_indices = np.in1d(sampleNames , lnba_sample_names)

LNBA_clade_indices_with_rise = np.in1d(sampleNames_projection , lnba_sample_names_with_rise)

YAC_clade_indices = np.in1d(sampleNames , ['YAC'])

# %%
# output updated MSA for tree building
# output full projection
calls_for_tree = apyp.idx2nts(calls_projection) # ATCGN translation

refgenome_nts = apyp.extract_outgroup_mutation_positions(ref_genome_folder, apyp.p2chrpos(p[goodpos],chrStarts));
refgenome_nts_for_tree=refgenome_nts
calls_for_tree_ref_outgroup = np.concatenate((apyp.idx2nts(refnti_m[:, 0].reshape(4444,1)),calls_for_tree),axis=1)
treesampleNamesLong = sampleNames_projection
treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)

apyp.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,f'{analysis_params_output_name}_full_qc_with_projection')


# output SNVs which are covered on RISE projection
calls_for_tree = apyp.idx2nts(calls_projection) # ATCGN translation

refgenome_nts = apyp.extract_outgroup_mutation_positions(ref_genome_folder, apyp.p2chrpos(p[goodpos],chrStarts));
refgenome_nts_for_tree=refgenome_nts
calls_for_tree_ref_outgroup = np.concatenate((apyp.idx2nts(refnti_m[:, 0].reshape(4444,1)),calls_for_tree),axis=1)
treesampleNamesLong = sampleNames_projection
treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)

called_in_projection=np.array(rise386_projection)!='?'
calls_for_tree_ref_outgroup_projection_only=calls_for_tree_ref_outgroup[called_in_projection]

apyp.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup_projection_only,treesampleNamesLong_ref_outgroup,f'{analysis_params_output_name}_full_qc_with_rise386_projection_only_calls')


# output SNVs which are covered on C90 projection
calls_for_tree = apyp.idx2nts(calls_projection) # ATCGN translation

refgenome_nts = apyp.extract_outgroup_mutation_positions(ref_genome_folder, apyp.p2chrpos(p[goodpos],chrStarts));
refgenome_nts_for_tree=refgenome_nts
calls_for_tree_ref_outgroup = np.concatenate((apyp.idx2nts(refnti_m[:, 0].reshape(4444,1)),calls_for_tree),axis=1)
treesampleNamesLong = sampleNames_projection
treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)

called_in_projection=np.array(c90_projection)!='?'
calls_for_tree_ref_outgroup_projection_only=calls_for_tree_ref_outgroup[called_in_projection]

apyp.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup_projection_only,treesampleNamesLong_ref_outgroup,f'{analysis_params_output_name}_full_qc_with_c90_projection_only_calls')


# output SNVs which are covered on ARK (no RISE386 projection)
# sanity check
calls_for_tree = apyp.idx2nts(calls) # ATCGN translation

refgenome_nts = apyp.extract_outgroup_mutation_positions(ref_genome_folder, apyp.p2chrpos(p[goodpos],chrStarts));
refgenome_nts_for_tree=refgenome_nts
calls_for_tree_ref_outgroup = np.concatenate((apyp.idx2nts(refnti_m[:, 0].reshape(4444,1)),calls_for_tree),axis=1)
treesampleNamesLong = sampleNames
treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)

called_in_projection=np.array(calls_for_tree[:,-1])!='?'
calls_for_tree_ref_outgroup_projection_only=calls_for_tree_ref_outgroup[called_in_projection]

apyp.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup_projection_only,treesampleNamesLong_ref_outgroup,f'{analysis_params_output_name}_full_qc_with_ark017_only_calls')


# output SNVs which are variable and covered on at least 2 LNBA samples, add call of YAC
variable_calls_lnba=[]
indices_to_get_yac_calls=[]
for idx,r in enumerate(calls_projection[:,LNBA_clade_indices_with_rise]):
    if len([x for x in np.unique(r) if x != 4])>1:
        variable_calls_lnba.append(r)
        indices_to_get_yac_calls.append(idx)

calls_for_tree = apyp.idx2nts(np.append(np.array(variable_calls_lnba),calls[indices_to_get_yac_calls][:,YAC_clade_indices],axis=1)) # ATCGN translation

treesampleNamesLong = sampleNames_projection[LNBA_clade_indices_with_rise]
treesampleNamesLong_outgroup = np.append(treesampleNamesLong,['YAC'])

apyp.write_calls_sampleName_to_fasta(calls_for_tree,treesampleNamesLong_outgroup,f'{analysis_params_output_name}_full_qc_with_projection_only_lnba_subset_calls')

# output SNVs which are variable and covered on at least 2 LNBA samples (NOT INCLUDING RISE), add call of YAC
variable_calls_lnba=[]
indices_to_get_yac_calls=[]
for idx,r in enumerate(calls_projection[:,LNBA_clade_indices_with_rise]):
    if len([x for x in np.unique(r) if x != 4])>1:
        variable_calls_lnba.append(r)
        indices_to_get_yac_calls.append(idx)

calls_for_tree = apyp.idx2nts(np.append(np.array(variable_calls_lnba),calls[indices_to_get_yac_calls][:,YAC_clade_indices],axis=1)) # ATCGN translation

treesampleNamesLong = sampleNames_projection[LNBA_clade_indices_with_rise]
treesampleNamesLong_outgroup = np.append(treesampleNamesLong,['YAC'])

apyp.write_calls_sampleName_to_fasta(calls_for_tree,treesampleNamesLong_outgroup,f'{analysis_params_output_name}_full_qc_with_projection_only_lnba_subset_calls')


# %%
# pairwise differences
def get_pairwise_genetic_distance(sample1,sample2,calls,sampleNames):
    col_index_to_compare=np.nonzero(np.in1d(sampleNames, sample1))[0]
    row_index_to_compare=np.nonzero(np.in1d(sampleNames, sample2))[0]
    calls_same_for_comparissons=np.sum(((calls[:,col_index_to_compare] == calls[:,row_index_to_compare]) & (calls[:,col_index_to_compare] != 4) & (calls[:,row_index_to_compare] != 4)))
    denominator_for_comparissons=np.sum(((calls[:,col_index_to_compare] !=4 ) & (calls[:,row_index_to_compare] != 4)))
    normalized_genetic_distance=(denominator_for_comparissons-calls_same_for_comparissons)/denominator_for_comparissons
    return normalized_genetic_distance,denominator_for_comparissons

# within LNBA
pairwise_dists=[]
pairwise_dists_comp_sample=[]
pairwise_comps=[]
for s in lnba_order:
    if s != 'RISE386':
        pairwise_dist_this_comp,pairwise_comps_this_comp=get_pairwise_genetic_distance(s,'RISE386',calls_projection,sampleNames_projection)
        pairwise_dists.append(pairwise_dist_this_comp)
        pairwise_comps.append(pairwise_comps_this_comp)
        pairwise_dists_comp_sample.append(s)
pairwise_comparisons_rise386_in_lnba=pd.DataFrame({'sample':pairwise_dists_comp_sample,'num_comparisons':pairwise_comps,'pairwise_dist':pairwise_dists,'sample_comps':[f'{x}, {y}' for x,y in zip(pairwise_dists_comp_sample,pairwise_comps)]})

log_data = np.log10(pairwise_comparisons_rise386_in_lnba[['sample_comps','pairwise_dist']].set_index('sample_comps')) # Use np.log1p to handle zero values

plt.subplots(facecolor='white',figsize=(1,20))
sns.heatmap(log_data,cmap='Blues',cbar_kws={'label': '$Log_{10}$ Pairwise Genetic Distance'})
plt.xticks([])
plt.ylabel('')
plt.savefig('rise386_pairwise_lnba.svg')


# outside LNBA
not_lnba_order=[x for x in tree_sample_order if x not in lnba_sample_names and x in sampleNames]

pairwise_dists=[]
pairwise_dists_comp_sample=[]
pairwise_comps=[]
for s in not_lnba_order:
    if s != 'RISE386':
        pairwise_dist_this_comp,pairwise_comps_this_comp=get_pairwise_genetic_distance(s,'RISE386',calls_projection,sampleNames_projection)
        pairwise_dists.append(pairwise_dist_this_comp)
        pairwise_comps.append(pairwise_comps_this_comp)
        pairwise_dists_comp_sample.append(s)
pairwise_comparisons_rise386_not_in_lnba=pd.DataFrame({'sample':pairwise_dists_comp_sample,'num_comparisons':pairwise_comps,'pairwise_dist':pairwise_dists,'sample_comps':[f'{x}, {y}' for x,y in zip(pairwise_dists_comp_sample,pairwise_comps)]})

log_data = np.log10(pairwise_comparisons_rise386_not_in_lnba[['sample_comps','pairwise_dist']].set_index('sample_comps')) # Use np.log1p to handle zero values

plt.subplots(facecolor='white',figsize=(1,40))
sns.heatmap(log_data,cmap='Blues',cbar_kws={'label': '$Log_{10}$ Pairwise Genetic Distance'})
plt.xticks([])
plt.ylabel('')
plt.savefig('rise386_pairwise_non_lnba.svg')


# %%
# ARK017 pairwise comparisons
# within LNBA
pairwise_dists=[]
pairwise_dists_comp_sample=[]
pairwise_comps=[]

lnba_order_with_rise=lnba_order[:lnba_order.index('ARK017')]+['RISE386']+lnba_order[lnba_order.index('ARK017'):]
for s in lnba_order_with_rise:
    if s != 'ARK017':
        pairwise_dist_this_comp,pairwise_comps_this_comp=get_pairwise_genetic_distance(s,'ARK017',calls_projection,sampleNames_projection)
        pairwise_dists.append(pairwise_dist_this_comp)
        pairwise_comps.append(pairwise_comps_this_comp)
        pairwise_dists_comp_sample.append(s)
pairwise_comparisons_ark017_in_lnba=pd.DataFrame({'sample':pairwise_dists_comp_sample,'num_comparisons':pairwise_comps,'pairwise_dist':pairwise_dists,'sample_comps':[f'{x}, {y}' for x,y in zip(pairwise_dists_comp_sample,pairwise_comps)]})

log_data = np.log10(pairwise_comparisons_ark017_in_lnba[['sample_comps','pairwise_dist']].set_index('sample_comps')) # Use np.log1p to handle zero values

plt.subplots(facecolor='white',figsize=(1,20))
sns.heatmap(log_data,cmap='Blues',cbar_kws={'label': '$Log_{10}$ Pairwise Genetic Distance'})
plt.xticks([])
plt.ylabel('')
plt.savefig('ark017_pairwise_lnba.svg')

# %%
# C90
not_lnba_order=[x for x in tree_sample_order if x not in lnba_sample_names and x in sampleNames]

pairwise_dists=[]
pairwise_dists_comp_sample=[]
for s in not_lnba_order:
    if s != 'C90':
        pairwise_dists.append(get_pairwise_genetic_distance(s,'C90',calls_projection,sampleNames_projection))
        pairwise_dists_comp_sample.append(s)

plt.subplots(facecolor='white',figsize=(1,20))
log_data = np.log10(pd.DataFrame.from_dict({x:y for x,y in zip(pairwise_dists_comp_sample,pairwise_dists)},orient='index'))  # Use np.log1p to handle zero values

sns.heatmap(log_data)
plt.xticks([])
plt.savefig('C90_pairwise_non_lnba.svg')

lnba_order=[x for x in tree_sample_order if x in lnba_sample_names]
pairwise_dists=[]
pairwise_dists_comp_sample=[]
for s in lnba_order:
    if s != 'C90':
        pairwise_dists.append(get_pairwise_genetic_distance(s,'C90',calls_projection,sampleNames_projection))
        pairwise_dists_comp_sample.append(s)

plt.subplots(facecolor='white',figsize=(1,20))
log_data = np.log10(pd.DataFrame.from_dict({x:y for x,y in zip(pairwise_dists_comp_sample,pairwise_dists)},orient='index'))  # Use np.log1p to handle zero values

sns.heatmap(log_data)
plt.xticks([])
plt.savefig('C90_pairwise_lnba.svg')


# pairwise LNBA and before
lnba_order=[x for x in tree_sample_order if x in np.append(lnba_sample_names,['RV2039','Gok2'])]
pairwise_dists=[]
pairwise_dists_comp_sample=[]
for s in lnba_order:
    if s != 'C90':
        pairwise_dists.append(get_pairwise_genetic_distance(s,'C90',calls_projection,sampleNames_projection))
        pairwise_dists_comp_sample.append(s)

plt.subplots(facecolor='white',figsize=(1,20))
log_data = np.log10(pd.DataFrame.from_dict({x:y for x,y in zip(pairwise_dists_comp_sample,pairwise_dists)},orient='index'))  # Use np.log1p to handle zero values

sns.heatmap(log_data)
plt.xticks([])
plt.savefig('C90_pairwise_lnba_and_prelnba.svg')

# %% 
# assessing number of pairwise SNVs from C90
# create ancestral reconstruction
calls_for_tree = apyp.idx2nts(calls_projection) # ATCGN translation

refgenome_nts = apyp.extract_outgroup_mutation_positions(ref_genome_folder, apyp.p2chrpos(p[goodpos],chrStarts));
refgenome_nts_for_tree=refgenome_nts
calls_for_tree_ref_outgroup = np.concatenate((apyp.idx2nts(refnti_m[:, 0].reshape(4444,1)),calls_for_tree),axis=1)
treesampleNamesLong = sampleNames_projection
treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)

tree_name_output=f'joint_projection/full_tree_alignment_1000bs_ML_all_production_run_fullqc_projected_pileup.raxml.bestTree'


apyp.create_ancestral_reconstruction(tree_name_output,treesampleNamesLong_ref_outgroup,calls_for_tree_ref_outgroup,f'joint_projection/{analysis_params_output_name}_ancestral_reconstruction')

tree=apyp.parse_tree(f'joint_projection/{analysis_params_output_name}_ancestral_reconstruction/annotated_tree.nexus')
ancestral_reconstruction_fasta=f'joint_projection/{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta'

from Bio import SeqIO
def get_internal_node_calls(internal_node_clade,ancestral_reconstruction_fasta):
    calls_pos={}
    for record in SeqIO.parse(ancestral_reconstruction_fasta, 'fasta'):
        calls_pos[record.id] = str(record.seq)
    return calls_pos[internal_node_clade.name]

# RISE386 differences to parent node (ARK017 parent)
ark_rise_parent=tree.common_ancestor(['ARK017','RISE386'])
ark_rise_parent_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(ark_rise_parent,ancestral_reconstruction_fasta)]))

rise386_calls=np.array([x for x in get_internal_node_calls([x for x in tree.find_elements('RISE386')][0],ancestral_reconstruction_fasta)])
rise386_ns=rise386_calls=='N'

print('Pairwise SNVs RISE386 and ARK017/RISE386 parent:', np.sum(rise386_calls[~(rise386_ns)]!=apyp.idx2nts(ark_rise_parent_clade_calls[~(rise386_ns)])))
# finding single SNP homoplasy:
print('single SNP homoplasy on RISE386:',p[np.where(~rise386_ns)[0][np.where(rise386_calls[~(rise386_ns)]!=apyp.idx2nts(ark_rise_parent_clade_calls[~(rise386_ns)]))[0]]]+1)

lnba_polytomy=tree.common_ancestor(['RK1001','RISE509'])
lnba_polytomy_parent_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(lnba_polytomy,ancestral_reconstruction_fasta)]))

pre_lnba_polytomy=tree.common_ancestor(['RK1001','RISE509','Gok2'])
pre_lnba_polytomy_parent_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(pre_lnba_polytomy,ancestral_reconstruction_fasta)]))
print('Pairwise SNVs pre-LNBA polytomy and LNBA polytomy:', np.sum(pre_lnba_polytomy_parent_clade_calls!=lnba_polytomy_parent_clade_calls))

post_lnba_poltomy=tree.common_ancestor(['VLI092','KNK001'])
post_lnba_poltomy_parent_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(post_lnba_poltomy,ancestral_reconstruction_fasta)]))
print('Pairwise SNVs post-LNBA polytomy and LNBA polytomy:', np.sum(post_lnba_poltomy_parent_clade_calls!=lnba_polytomy_parent_clade_calls))

post_lnba_internal_node=tree.common_ancestor(['0.PE7b_620024','I2470'])
post_lnba_internal_node_parent_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(post_lnba_internal_node,ancestral_reconstruction_fasta)]))
print('Pairwise SNVs post-LNBA internal node and LNBA polytomy:', np.sum(post_lnba_internal_node_parent_clade_calls!=lnba_polytomy_parent_clade_calls))


c90_calls=np.array([x for x in get_internal_node_calls([x for x in tree.find_elements('C90')][0],ancestral_reconstruction_fasta)])
c90_ns=c90_calls=='N'

print('Pairwise SNVs C90 and LNBA polytomy:', np.sum(c90_calls[~(c90_ns)]!=apyp.idx2nts(lnba_polytomy_parent_clade_calls[~(c90_ns)])))

rk1001_calls=np.array([x for x in get_internal_node_calls([x for x in tree.find_elements('RK1001')][0],ancestral_reconstruction_fasta)])
rk1001_ns=rk1001_calls=='N'
print('Pairwise SNVs RK1001 and C90:', np.sum(rk1001_calls[~(c90_ns|rk1001_ns)]!=c90_calls[~(c90_ns|rk1001_ns)]))

VLI092_calls=np.array([x for x in get_internal_node_calls([x for x in tree.find_elements('VLI092')][0],ancestral_reconstruction_fasta)])
VLI092_ns=VLI092_calls=='N'
print('Pairwise SNVs VLI092 and C90:', np.sum(VLI092_calls[~(c90_ns|VLI092_ns)]!=c90_calls[~(c90_ns|VLI092_ns)]))

Gyvakarai1_calls=np.array([x for x in get_internal_node_calls([x for x in tree.find_elements('Gyvakarai1')][0],ancestral_reconstruction_fasta)])
Gyvakarai1_ns=Gyvakarai1_calls=='N'
print('Pairwise SNVs Gyvakarai1 and C90:', np.sum(Gyvakarai1_calls[~(c90_ns|Gyvakarai1_ns)]!=c90_calls[~(c90_ns|Gyvakarai1_ns)]))


# %% pairwise across internal nodes
pairwise_dists=[]
pairwise_dists_comp_sample=[]

for n in tree.get_nonterminals():
    internal_call=apyp.nts2idx(np.array([x for x in get_internal_node_calls(n,ancestral_reconstruction_fasta)]))
    pairwise_dists.append(np.sum(internal_call==calls_projection[:,sampleNames_projection=='C90'].flatten()))
    pairwise_dists_comp_sample.append(n.name)

plt.subplots(facecolor='white',figsize=(1,40))


sns.heatmap(pd.DataFrame.from_dict({x:y/4444 for x,y in zip(pairwise_dists_comp_sample,pairwise_dists)},orient='index'))
plt.xticks([])
print(log_data.index[np.where(log_data==log_data.max())[0]])
