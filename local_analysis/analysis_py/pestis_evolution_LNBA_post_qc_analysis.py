# %%
# =============================================================================
## ARK017-LNBA analysis, Ian Light-Maka
# =============================================================================
# SNP-based analyses:
# =============================================================================
# - Most analysis from apy (analysis.py) module
# - parsimony treeas
# - SNP-specific tree coloring
# - parallel evolution analysis
# - molecular clock analysis

# =============================================================================

## import libraries not gotten from analysispy_module
from datetime import date
from importlib import reload
import sys
import os
import re
import glob
import pickle
from tokenize import group
import numpy as np
from scipy import io
from scipy import stats
import pandas as pd
from Bio import SeqIO, Seq, SeqRecord
import matplotlib.pyplot as plt
import math
from statsmodels.stats.multitest import multipletests
import datetime
import time
import upsetplot
import seaborn as sns
import gzip
from matplotlib.patches import Rectangle

## declare paths for finding analysispy_modules and reference genome folder

os.chdir('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/production_run')

ref_genome_folder='/Users/ad_loris/Nextcloud/keylab/reference_genomes/Ypestis_ASM906v1'

sys.path.append("/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/analysis_py")
import analysispy_modules_pestis as apyp

cmtFile_sub_gp = 'candidate_mutation_table_processed_gp.pickle.gz'

ancient_sample_names=['GEN72','GRS004','KLE031','KLE048','GRH001','VEL003','HGC068','GLZ001','GLZ002','KunilaII','Barcelona_3031','EDI001.A','ES_11972','ES_8291','1343UnTal85','6Post','ARS007','CHC004','Gok2','Gyvakarai1','HGC009','HOP001','HOP004','I2470','I5884','KLZ001','KNK001','KZL002','MIB054','RISE505','RISE509','RK1001','RT5','RT6','RV2039','VLI092','AE1175','BED024.A0102','BED028.A0102','BED030.A0102','BED034.A0102','Bolgar_2370','BRA001.A0101','Ellwangen_549_O','LAI009.A0101','MAN008.B0101','NMS002.A0101','OBS107','OBS110','OBS116','OBS124','OBS137','STA001.A0101','STN002.A0101','STN007.A0101','STN008.A0101','STN014.A0101','STN019.A0101','STN020.A0101','STN021.A0101','OOH003','ARK017', 'DSH008', 'DSH025','C10091', 'C10098', 'C10928']
promotersize=250
NTs = np.array(['A','T','C','G'],dtype=object) # NTs='ATCG'

ts=time.time()
timestamp=datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M')

# %% 
# load data from QC file
with gzip.open('final_cmt/' + cmtFile_sub_gp, 'rb') as pickle_file:
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
num_samples=len(sampleNames)

# Load relevant metadata 
#   Isolation source for modern genomes
samples_to_isolation=pd.read_csv('/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/metadata_tables/samples_to_isolation_source.tsv',sep='\t',header=0)

#   LNBA and pre-LNBA metadata
lnba_metadata_path='/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/metadata_tables/lnba_prelnba_metadata.tsv'
lnba_metadata=pd.read_csv(lnba_metadata_path,sep='\t')

#   Info from genome files
[chrStarts, genomeLength, scafNames] = apyp.genomestats(ref_genome_folder);

# %% 
# Extract refnt and define out/in-group bools
# =============================================================================

## Note ancnti/outs_nti defined below after filtered calls has been generated!
## get reference allele for all p; NOTE: analysis.m stored ATCG as double-digit-numeric
# use ref nt for mutation calls. important if multiple outgroups called 
refnt = apyp.extract_outgroup_mutation_positions(ref_genome_folder, apyp.p2chrpos(p,chrStarts));
refnti = apyp.nts2idx(refnt.copy())
refnti_m = np.tile(refnti,(len(sampleNames),1)).transpose() # build 2D matrix with outgroup (ancestral) allele
# print(np.unique(refnt)) # sanity check

# When no outgroup defined: refnt ~= ancnt:
#ancnt = refnt   

## Estimate outgroup (ancestral allele) from ALL samples added as outgroup to SM pipeline (ancnti* == major allele)
##  NOTE: NOTE: NOTE:

## TODO: update how this is done
##  NOTE: 
# CAN ALSO JUST GRAB OUTGROUP BY OUTPUTTING THE SAMPLES CSV from case
outgroup_bool = np.isin(sampleNames,['YAC','Ypseudo_FDA_ARGOS_665'])
outgroup_name = sampleNames[outgroup_bool]
outgroup_idx=np.nonzero(outgroup_bool)[0]

# ingroup array (bool, idx) used later
ingroup_bool = np.invert(outgroup_bool)
ingroup_idx = np.nonzero(ingroup_bool)[0]

# %%
# Generate ancestral reconstruction based on ML tree
# ancestral reconstruction on goodpos, MRCA of all modern, LNBA samples (node)
tree_name_output=f'raxml_final/full_tree_alignment_1000bs_ML_all_production_run_fullqc.raxml.bestTree'
treesampleNamesLong_ref_outgroup = np.append(['Sref'] ,sampleNames )
calls_for_tree_ref_outgroup =  np.concatenate((refnt[:,None], apyp.idx2nts(calls)),axis=1)


apyp.create_ancestral_reconstruction(tree_name_output,treesampleNamesLong_ref_outgroup,calls_for_tree_ref_outgroup,f'{analysis_params_output_name}_ancestral_reconstruction')

tree=apyp.parse_tree(f'{analysis_params_output_name}_ancestral_reconstruction/annotated_tree.nexus')

LNBA_samples=['RK1001', 'RISE509', 'VLI092', 'KNK001', 'Gyvakarai1', 'GEN72',
       'GRS004', 'VEL003', 'I5884', 'KunilaII', 'GLZ002', 'GLZ001',
       '1343UnTal85', 'HOP004', 'HOP001', 'KLZ001', 'CHC004', '6Post',
       'HGC009', 'KLE031', 'MIB054', 'KLE048', 'OOH003', 'RISE505',
       'ARS007', 'KZL002', 'GRH001','ARK017', 'DSH008', 'DSH025','C10091', 'C10098', 'C10928']

MRCA_all_pestis=tree.common_ancestor(list(sampleNames[(sampleNames != 'YAC') & (sampleNames != 'Ypseudo_FDA_ARGOS_665')]))
MRCA_lnba=tree.common_ancestor(list(np.intersect1d(sampleNames,LNBA_samples)))
modern_plus_lnba=list(filter(re.compile(r'^[0-4]{1}\.').search,sampleNames))+LNBA_samples
modern_samples=list(filter(re.compile(r'^[0-4]{1}\.').search,sampleNames))
MRCA_pestis_lnba=tree.common_ancestor(list(np.intersect1d(sampleNames,modern_plus_lnba)))
LNBA_samples=sampleNames[np.isin(sampleNames,LNBA_samples)]

MRCA_big_bang=tree.common_ancestor(list(filter(re.compile(r'^[1-4]{1}\.').search,sampleNames)))

mrca_big_bang_reconstruction_nti_goodpos=apyp.parse_ancestra_fasta_treetime(f'{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta',MRCA_big_bang.name)
ancestral_reconstruction_nti_goodpos=apyp.parse_ancestra_fasta_treetime(f'{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta',MRCA_pestis_lnba.name)

ancestral_reconstruction_nti=refnti.copy()
ancestral_reconstruction_nti[goodpos]=ancestral_reconstruction_nti_goodpos
ancestral_reconstruction_nti_m=np.tile(ancestral_reconstruction_nti,(num_samples,1)).transpose()

mrca_all_pestis_reconstruction_nti_goodpos=apyp.parse_ancestra_fasta_treetime(f'{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta',MRCA_all_pestis.name)

mrca_all_pestis_reconstruction_nti=refnti.copy()
mrca_all_pestis_reconstruction_nti[goodpos]=mrca_all_pestis_reconstruction_nti_goodpos
mrca_all_pestis_reconstruction_nti_m=np.tile(mrca_all_pestis_reconstruction_nti,(num_samples,1)).transpose()

mrca_big_bang_reconstruction_nti_goodpos=apyp.parse_ancestra_fasta_treetime(f'{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta',MRCA_big_bang.name)

mrca_big_bang_reconstruction_nti=refnti.copy()
mrca_big_bang_reconstruction_nti[goodpos]=mrca_big_bang_reconstruction_nti_goodpos
mrca_big_bang_reconstruction_nti_m=np.tile(mrca_big_bang_reconstruction_nti,(num_samples,1)).transpose()

ancestral_reconstruction_tree=f'{analysis_params_output_name}_ancestral_reconstruction/annotated_tree.nexus'
ancestral_reconstruction_fasta=f"{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta"

#apyp.create_full_genome_ancestral_reconstruction_fasta(ancestral_reconstruction_nti,p,chrStarts,ref_genome_folder)


def get_external_nodes_nodes_bfs(tree,internal_node_name):
    labelled_internal_node_tree=apyp.parse_tree(tree)
    nodes_to_check=list(labelled_internal_node_tree.find_clades(internal_node_name))
    unique_clades_to_output=[]
    while len(nodes_to_check)>0:
        checking_node = nodes_to_check[0]
        nodes_to_add = [x for x in checking_node.clades]
        nodes_to_check += nodes_to_add
        nodes_to_check = nodes_to_check[1:]
        if checking_node.is_terminal():
            unique_clades_to_output.append(checking_node)
    return unique_clades_to_output

post_big_bang_samples = [x.name for x in  get_external_nodes_nodes_bfs(f'{analysis_params_output_name}_ancestral_reconstruction/annotated_tree.nexus',MRCA_big_bang.name)]

ancient_sample_bool = np.in1d(sampleNames,ancient_sample_names)

tree_sample_order=[x.name for x in tree.get_terminals()]
tree_sample_order_array=np.array(tree_sample_order)

LNBA_tree_order = np.array([x for x in tree_sample_order_array if x in LNBA_samples])
np.savetxt('tree_order.txt',tree_sample_order,fmt='%s')


clades={}
for sample in LNBA_samples:
    clades[sample]=apyp.find_clade_terminals(sample,tree_name_output)
lnba_samples_base_name=[i.split('_')[0] for i in LNBA_samples]

np.savetxt('lnba_samples_base_name.txt',lnba_samples_base_name,fmt='%s')

LNBA_clade_indices = np.in1d(sampleNames , LNBA_samples)
LNBA_and_before = np.in1d(sampleNames , np.concatenate([LNBA_samples,['YAC','Ypseudo_FDA_ARGOS_665','Gok2','RV2039']]))

all_samples_post_MRCA_modern_LNBA=~np.in1d(sampleNames,['YAC','Ypseudo_FDA_ARGOS_665','Gok2','RV2039'])

# %%
# Finding homoplasic sites, QC 
# =============================================================================
# find number of mutational events based on tree structure
os.makedirs(os.getcwd() + f"/pdf/qc/homoplasies", exist_ok=True)
num_mutational_events_goodpos_terminal_nodes,node_tracking_terminal = apyp.count_number_mutational_events(f'{analysis_params_output_name}_ancestral_reconstruction/annotated_tree.nexus',f'{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta',False,True)
num_mutational_events_goodpos_internal_nodes,node_tracking_internal = apyp.count_number_mutational_events(f'{analysis_params_output_name}_ancestral_reconstruction/annotated_tree.nexus',f'{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta',True,True)

# %%
# Make a tree for each SNP location, coloring each tip with given SNP call

## create empty tree_counting folder and add for_tree_labeling.csv

## add tree for each mutation dsiplaying mutation (generate_mutation_colored_tree.py)
apyp.build_table_for_tree_labeling(apyp.p2chrpos(p[goodpos], chrStarts),treesampleNamesLong_ref_outgroup,calls_for_tree_ref_outgroup)
os.chdir('tree_counting/snps')

apyp.generate_mutation_colored_tree(f"../../{tree_name_output}","for_tree_labeling.csv", outgroup_name[0], num_mutational_events_goodpos_terminal_nodes, True, True, False)
os.chdir('../..')

## END TREE PRODUCTION, BEGIN SNP-LEVEL ANALYSIS
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# %%  
# Store SNP table
# =============================================================================
## SOM SNP Table. Should be a function.
## get NT,sampleName df
# build sampleNames w/ metainfo
# snp_table.csv contains 'nt_anc' > nt of explicit single outgroup isolate. '.' == NA
sampleNames_ingroup = sampleNames[ingroup_bool] # remove outgroup

# translate index to nucleotide
calls_for_treei = calls; 
calls_for_treei = calls_for_treei[ goodpos, : ]
calls_for_treei_ingroup = calls_for_treei[ : , ingroup_bool ]
calls_for_tree = apyp.idx2nts(calls_for_treei_ingroup) # ATCGN translation
snp_data = pd.DataFrame(calls_for_tree,columns=sampleNames_ingroup)

# store full tree alignment, LNBA and before alignment
apyp.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,f'{analysis_params_output_name}_full_tree_alignment_ref_outgroup')
apyp.write_calls_sampleName_to_fasta(apyp.idx2nts(calls_for_treei),sampleNames,f'{analysis_params_output_name}_full_tree_alignment')

# %%
# annotate mutations
calls_outgroup = calls[:,outgroup_bool] 
ancnti = apyp.major_allele(calls_outgroup) 

ancnti_m = np.tile(ancnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    

outsplnti_m=ancnti_m

[mutQual, mutQualIsolates] = apyp.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 

## generate mutation annotation from gff
annotation_genes = apyp.parse_gff(ref_genome_folder,scafNames,forceReDo=True) # ref_genome_folder+"/annotation_genes.pandas.py.pk1"
annotation_mutations = apyp.annotate_mutations(annotation_genes , p , refnti_m, outsplnti_m, calls, counts, hasmutation, mutQual, promotersize , ref_genome_folder, goodpos) # extract relevant annotation info for each SNP    

hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis = (calls != ancestral_reconstruction_nti_m) & (np.tile(mutQual,(1,num_samples)) >= 1) & (calls != 4)
hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,outgroup_bool]=False
goodpos_relative_ancestral_reconstruction_mrca_lnba = np.intersect1d(goodpos,np.where( np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis, axis=1) > 0 )[0]) # NOTE: candpos/goodpos is INDEX of good positions for p!

hasmutation_relative_ancestral_MRCA_all_pestis =  (calls != mrca_all_pestis_reconstruction_nti_m) & (np.tile(mutQual,(1,num_samples)) >= 1) & (calls != 4)
hasmutation_relative_ancestral_MRCA_all_pestis[:,outgroup_bool]=False
goodpos_relative_ancestral_reconstruction = np.intersect1d(goodpos,np.where( np.sum(hasmutation_relative_ancestral_MRCA_all_pestis, axis=1) > 0 )[0]) # NOTE: candpos/goodpos is INDEX of good positions for p!

annotation_mutations_relative_MRCA = apyp.annotate_mutations(annotation_genes , p[goodpos_relative_ancestral_reconstruction] , refnti_m[goodpos_relative_ancestral_reconstruction,:] , ancestral_reconstruction_nti_m[goodpos_relative_ancestral_reconstruction,:] , calls[goodpos_relative_ancestral_reconstruction,:] , counts[:,:,goodpos_relative_ancestral_reconstruction] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[goodpos_relative_ancestral_reconstruction,:], mutQual[goodpos_relative_ancestral_reconstruction,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    

# %% get snp metadata
# Contig	Location	Cluster ID	NCTC9343 homolog	Assembly locus	Gene	Protein annotation (Prokka)	Type*	Ancestor allele	Amino acid changes	Nucleotide position in the gene**
snp_metadata = annotation_mutations.copy()
# for I/P mutations: turn all to I (P somewhat random); add up/downstream info to 0tag and gene description to product
for i,row in snp_metadata.iterrows(): # I/P SNV
    if row.isnull()['locustag']:
        if row.isnull()['locustag1'] and not row.isnull()['locustag2']: # no preceding gene. start of contig
            snp_metadata.at[i,'locustag'] = "NA;"+str(int(row['distance2']))+":"+row['locustag2']
            snp_metadata.at[i,'gene'] = "NA;"+row['gene2']
            snp_metadata.at[i,'product'] = "NA;"+row['product2']
            snp_metadata.at[i,'type'] = 'I'
        elif row.isnull()['locustag2'] and not row.isnull()['locustag1']: # no subsequent gene. end of contig
            snp_metadata.at[i,'locustag'] = str(row['distance1'])+":"+str(row['locustag1'])+";NA"
            snp_metadata.at[i,'gene'] = row['gene1']+";NA"
            snp_metadata.at[i,'product'] = row['product1']+";NA"
            snp_metadata.at[i,'type'] = 'I'
        elif row.isnull()['locustag1'] and row.isnull()['locustag2']: # no annotated gene on contig
            snp_metadata.at[i,'locustag'] = "NA"
            snp_metadata.at[i,'gene'] = "NA"
            snp_metadata.at[i,'product'] = "NA"
            snp_metadata.at[i,'type'] = 'I'
        else: # intergenic SNV with preceding/subsequent gene                
            snp_metadata.at[i,'locustag'] = str(row['distance1'])+":"+str(row['locustag1'])+";"+str(int(row['distance2']))+":"+row['locustag2']
            snp_metadata.at[i,'gene'] = row['gene1']+";"+row['gene2']
            snp_metadata.at[i,'product'] = row['product1']+";"+row['product2']
            snp_metadata.at[i,'type'] = 'I'

snp_metadata = snp_metadata[['chr','pos','type','muts','locustag','gene','loc1','loc2','strand','product','nt_pos','aa_pos','nt_ref','nt_alt','nt_anc']]
snp_metadata[['nt_anc']] = snp_metadata[['nt_anc']].replace('.','?') # turn NA to '?', similar to NT data 
snp_table = pd.concat([snp_metadata,snp_data.reset_index(drop=True)], axis=1)

## store
with open(f'{analysis_params_output_name}_snp_table.csv', 'w') as file:
    snp_table.to_csv(file, header=True, index=False)

# %% 
# Estimate substitution rate
# =============================================================================
allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG'],dtype=object)
# AT, TA  0
# AC, TG  1
# AG, TC  2
# GC, CG  3
# GT, CA  4
# GA, CT  5
allmuts_types = np.array([0,2,1,0,1,2,5,4,3,4,5,3])
mutationalspectrum = [1/12] * 12 # uniform distribution. replace with time stratified observation based on all samples

## loop over all ref/anc nucleotides and count mutation types. count twice for each direction
allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG'],dtype=object)
allmuts_counts = np.zeros( allmuts.shape,dtype=int )
for i,ref in enumerate(apyp.idx2nts(refnti[goodpos])):
    obs_allele = np.unique(calls_for_tree[i,1:])
    obs_allele = obs_allele[ ~(obs_allele == '?') & ~(obs_allele == ref) ]
    for j in obs_allele:
        idx = np.where( (ref+j)==allmuts)
        allmuts_counts[idx] += 1
        idx = np.where( (j+ref)==allmuts)
        allmuts_counts[idx] += 1
mutationalspectrum = allmuts_counts/len(calls_for_tree[:,0])/2
# store
afile = open(f"{analysis_params_output_name}_mutationalspectrum.py.pk1", 'wb')
pickle.dump(mutationalspectrum, afile)
afile.close()
# %% adding homoplasy calls
# Out of python, generate prokka annotate ancestral reconstruction annotation
#anc_reconstruction_folder = '/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/all/snppar/ancestral_reconstruction'
#ancestral_fasta = apyp.extract_outgroup_mutation_positions(anc_reconstruction_folder, apyp.p2chrpos(p,chrStarts));

#annotation_genes_ancestral_reconstruction=apyp.parse_gff('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/all/ancestral_reconstruction_backup/prokka/',scafNames,forceReDo=True)

#annotation_mutations_ancestral=apyp.annotate_mutations(annotation_genes_ancestral_reconstruction , p , refnti_m , ancestral_reconstruction_nti_m , calls , counts , hasmutation, mutQual , promotersize , ref_genome_folder, goodpos) # extract relevant annotation info for each SNP    

## update annotation mutations with this info
annotation_mutations = apyp.annotate_mutations(annotation_genes , p , refnti_m, outsplnti_m, calls, counts, hasmutation, mutQual, promotersize , ref_genome_folder, goodpos) # extract relevant annotation info for each SNP    
annotation_mutations.insert(27,'num_mutational_events',num_mutational_events_goodpos_terminal_nodes)
annotation_mutations.insert(27,'num_mutational_events_internal_nodes',num_mutational_events_goodpos_internal_nodes)
#annotation_mutations.insert(27,'homoplasy_bool',homoplasy_bool)


annotation_mutations_relative_MRCA.insert(27,'num_mutational_events',num_mutational_events_goodpos_terminal_nodes)
# %%
# generate SNP analysis sets
def get_internal_nodes_bfs(tree,internal_node_name):
    labelled_internal_node_tree=apyp.parse_tree(tree)
    nodes_to_check=list(labelled_internal_node_tree.find_clades(internal_node_name))
    unique_clades_to_output=[]
    while len(nodes_to_check)>0:
        checking_node = nodes_to_check[0]
        print(checking_node)
        nodes_to_add = [x for x in checking_node.clades if not x.is_terminal()]
        nodes_to_check += nodes_to_add
        nodes_to_check = nodes_to_check[1:]
        if checking_node.branch_length != 0:
            unique_clades_to_output.append(checking_node)
    return unique_clades_to_output

def get_distance_from_clade_to_clade(parsed_tree,cladeA,cladeB):
    traced=parsed_tree.trace(cladeA,cladeB)
    return len(traced)-1

def get_internal_node_calls(internal_node_clade,ancestral_reconstruction_fasta):
    calls_pos={}
    for record in SeqIO.parse(ancestral_reconstruction_fasta, 'fasta'):
        calls_pos[record.id] = str(record.seq)
    return calls_pos[internal_node_clade.name]

# %% generate SNP distances separating 83bootstrap ARK polytomy vs 17 bootstrap polytomy
ark_83_support_parent_clade=tree.common_ancestor(['ARK017','KZL002'])
ark_17_support_parent_clade=tree.common_ancestor(['KLE031','KZL002'])
ark_83_support_parent_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(ark_83_support_parent_clade,ancestral_reconstruction_fasta)]))
ark_17_support_parent_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(ark_17_support_parent_clade,ancestral_reconstruction_fasta)]))
print('SNVs separating polytomies:',np.sum(ark_83_support_parent_clade_calls!=ark_17_support_parent_clade_calls))

# TRUE BACKBONE SNPS
# lnba backbone
last_internal_node_on_lnba_backbone=tree.common_ancestor(['GRH001','KZL002'])
lnba_backbone_clades=tree.trace(MRCA_pestis_lnba,last_internal_node_on_lnba_backbone)

backbone_snps=set()
for clade in lnba_backbone_clades:
    current_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(clade,ancestral_reconstruction_fasta)]))
    has_mut_current_clade=(ancestral_reconstruction_nti[goodpos]!=current_clade_calls)
    backbone_snps.update(set(np.where(has_mut_current_clade)[0]))
backbone_snps_all_pos=goodpos[np.array(sorted(backbone_snps))]
non_singletons_lnba=np.intersect1d(goodpos,np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,LNBA_clade_indices],axis=1)>1)[0])
backbone_snps_all_pos_non_singletons=np.intersect1d(non_singletons_lnba,backbone_snps_all_pos)
annotation_mutations_lnba_true_backbone_relative_MRCA=apyp.annotate_mutations(annotation_genes , p[backbone_snps_all_pos_non_singletons] , refnti_m[np.ix_(backbone_snps_all_pos_non_singletons,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(backbone_snps_all_pos_non_singletons,np.full((num_samples), True))] , calls[np.ix_(backbone_snps_all_pos_non_singletons,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),backbone_snps_all_pos_non_singletons)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(backbone_snps_all_pos_non_singletons,np.full((num_samples), True))], mutQual[backbone_snps_all_pos_non_singletons,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
#annotation_mutations_lnba_true_backbone_relative_MRCA=annotation_mutations_lnba_true_backbone_relative_MRCA[~annotation_mutations_lnba_true_backbone_relative_MRCA['product'].isin(genes_to_drop)]


# INTERNAL LNBA SNPS
candidate_backbone_snps=set()
for clade_key in clades:
    clade_sample_names=clades[clade_key]
    clade_indices=np.in1d(sampleNames , clade_sample_names)
    # find mutations that occur on all the samples in the clade (except one to allow for uncalled), add to growing set of positings in p to investigate
    backbone_mutations=np.where( (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,clade_indices], axis=1) + np.sum((calls[:,clade_indices]==4),axis=1)) == len(clade_sample_names))[0]
    backbone_mutations_not_all_ns_lnba=backbone_mutations[np.sum((calls[backbone_mutations][:,clade_indices] < 4),axis=1) > 1]
    candidate_backbone_snps.update(backbone_mutations_not_all_ns_lnba)
candidate_backbone_snps_in_goodpos_relative_MRCA=np.intersect1d(np.array(list(candidate_backbone_snps)),goodpos)

candidate_backbone_snps_in_goodpos_relative_MRCA=candidate_backbone_snps_in_goodpos_relative_MRCA
annotation_mutations_lnba_internal_calls_all_relative_MRCA=apyp.annotate_mutations(annotation_genes , p[candidate_backbone_snps_in_goodpos_relative_MRCA] , refnti_m[np.ix_(candidate_backbone_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(candidate_backbone_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , calls[np.ix_(candidate_backbone_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),candidate_backbone_snps_in_goodpos_relative_MRCA)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(candidate_backbone_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))], mutQual[candidate_backbone_snps_in_goodpos_relative_MRCA,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
candidate_backbone_snps_in_goodpos_relative_MRCA_no_ark=np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,np.isin(sampleNames,LNBA_samples[~np.isin(LNBA_samples,['ARK017'])])],axis=1)>1)[0]
annotation_mutations_lnba_internal_calls_all_relative_MRCA_no_ark=apyp.annotate_mutations(annotation_genes , p[candidate_backbone_snps_in_goodpos_relative_MRCA_no_ark] , refnti_m[np.ix_(candidate_backbone_snps_in_goodpos_relative_MRCA_no_ark,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(candidate_backbone_snps_in_goodpos_relative_MRCA_no_ark,np.full((num_samples), True))] , calls[np.ix_(candidate_backbone_snps_in_goodpos_relative_MRCA_no_ark,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),candidate_backbone_snps_in_goodpos_relative_MRCA_no_ark)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(candidate_backbone_snps_in_goodpos_relative_MRCA_no_ark,np.full((num_samples), True))], mutQual[candidate_backbone_snps_in_goodpos_relative_MRCA_no_ark,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
#annotation_mutations_lnba_internal_calls_all_relative_MRCA_no_ark=annotation_mutations_lnba_internal_calls_all_relative_MRCA_no_ark[~annotation_mutations_lnba_internal_calls_all_relative_MRCA_no_ark['product'].isin(genes_to_drop)]

# EXCLUSIVE ARK SNVS
candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_only_ark=  np.where( ((np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,np.isin(sampleNames,LNBA_samples[np.isin(LNBA_samples,['ARK017'])])], axis=1) ) == 1 ) & ( (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,~np.isin(sampleNames,LNBA_samples[~np.isin(LNBA_samples,['ARK017'])])],axis=1)==0)))[0]
annotation_mutations_lnba_singletons_relative_MRCA_only_ark=apyp.annotate_mutations(annotation_genes , p[candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_only_ark] , refnti_m[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_only_ark,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_only_ark,np.full((num_samples), True))] , calls[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_only_ark,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_only_ark)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_only_ark,np.full((num_samples), True))], mutQual[candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_only_ark,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    


np.where( (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,LNBA_clade_indices], axis=1) >1) & ( (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,~LNBA_clade_indices],axis=1)==0)) )[0]


## exclusive LNBA
not_on_non_lnba=np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,:],axis=1)-np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,:][:,LNBA_clade_indices],axis=1)==0)[0]
on_lnba=np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,LNBA_clade_indices],axis=1) > 0)[0]
exclusive_lnba=np.intersect1d(goodpos,np.intersect1d(on_lnba,not_on_non_lnba))
annotation_mutations_exclusive_lnba=apyp.annotate_mutations(annotation_genes , p[exclusive_lnba] , refnti_m[np.ix_(exclusive_lnba,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(exclusive_lnba,np.full((num_samples), True))] , calls[np.ix_(exclusive_lnba,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),exclusive_lnba)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(exclusive_lnba,np.full((num_samples), True))], mutQual[exclusive_lnba,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
lnba_mutation_counts_for_pos=[]
for index,row in annotation_mutations_exclusive_lnba.iterrows():
    row_chr_start=chrStarts[row['chr']-1]
    pos=row_chr_start+row['pos']
    index_in_hasmutation=np.where(p==pos-1)[0]
    lnba_mutation_counts_for_pos.append(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[index_in_hasmutation,:][:,LNBA_clade_indices]))
annotation_mutations_exclusive_lnba.insert(27,'num_LNBA_samples_mutated',lnba_mutation_counts_for_pos)
#annotation_mutations_exclusive_lnba=annotation_mutations_exclusive_lnba[~annotation_mutations_exclusive_lnba['product'].isin(genes_to_drop)]

# general LNBA snps
has_mut_lnba=np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,:][:,LNBA_clade_indices],axis=1)>0)[0]
annotation_mutations_lnba=apyp.annotate_mutations(annotation_genes , p[has_mut_lnba] , refnti_m[np.ix_(has_mut_lnba,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(has_mut_lnba,np.full((num_samples), True))] , calls[np.ix_(has_mut_lnba,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),has_mut_lnba)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(has_mut_lnba,np.full((num_samples), True))], mutQual[has_mut_lnba,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    


## exluding LNBA
has_mut_not_lnba=np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,:][:,~LNBA_clade_indices],axis=1)>0)[0]
has_no_mut_lnba=np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,:][:,LNBA_clade_indices],axis=1)==0)[0]
exclusive_non_lnba=np.intersect1d(has_mut_not_lnba,has_no_mut_lnba)
annotation_mutations_exclusive_non_lnba=apyp.annotate_mutations(annotation_genes , p[exclusive_non_lnba] , refnti_m[np.ix_(exclusive_non_lnba,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(exclusive_non_lnba,np.full((num_samples), True))] , calls[np.ix_(exclusive_non_lnba,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),exclusive_non_lnba)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(exclusive_non_lnba,np.full((num_samples), True))], mutQual[exclusive_non_lnba,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
#annotation_mutations_exclusive_non_lnba=annotation_mutations_exclusive_non_lnba[~annotation_mutations_exclusive_non_lnba['product'].isin(genes_to_drop)]

# modern mutations
#annotation_mutations_exclusive_non_lnba.insert(27,'num_mutational_events',multiply_mutated_pos[exclusive_modern])
has_mut_modern=np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,:][:,np.isin(sampleNames,modern_samples)],axis=1)>0)[0]
annotation_mutations_modern=apyp.annotate_mutations(annotation_genes , p[has_mut_modern] , refnti_m[np.ix_(has_mut_modern,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(has_mut_modern,np.full((num_samples), True))] , calls[np.ix_(has_mut_modern,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),has_mut_modern)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(has_mut_modern,np.full((num_samples), True))], mutQual[has_mut_modern,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    



## adding to goodpos all relative MRCA
lnba_mutation_counts_for_pos=[]
total_mutation_counts_for_pos=[]
non_lnba_mutation_counts_for_pos=[]
for index,row in annotation_mutations_lnba_internal_calls_all_relative_MRCA.iterrows():
    row_chr_start=chrStarts[row['chr']-1]
    pos=row_chr_start+row['pos']
    index_in_hasmutation=np.where(p==pos-1)[0]
    total_mutation_counts_for_pos.append(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[index_in_hasmutation,:][:,~outgroup_bool]))
    lnba_mutation_counts_for_pos.append(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[index_in_hasmutation,:][:,LNBA_clade_indices]))
    non_lnba_mutation_counts_for_pos.append(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[index_in_hasmutation,:][:,~outgroup_bool])-np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[index_in_hasmutation,:][:,LNBA_clade_indices]))

annotation_mutations_lnba_internal_calls_all_relative_MRCA.insert(27,'num_LNBA_samples_mutated',lnba_mutation_counts_for_pos)
annotation_mutations_lnba_internal_calls_all_relative_MRCA.insert(27,'num_all_samples_mutated',total_mutation_counts_for_pos)
annotation_mutations_lnba_internal_calls_all_relative_MRCA.insert(27,'num_non_LNBA_samples_mutated',non_lnba_mutation_counts_for_pos)
#annotation_mutations_lnba_internal_calls_all_relative_MRCA=annotation_mutations_lnba_internal_calls_all_relative_MRCA[~annotation_mutations_lnba_internal_calls_all_relative_MRCA['product'].isin(genes_to_drop)]

num_genes=0
for g in annotation_genes:
    num_genes+=len(g)

# non_lnba_backbone to big-bang
modern_backbone_clades=tree.trace(MRCA_pestis_lnba.name,MRCA_big_bang.name)

modern_backbone_snps=set()
for clade in modern_backbone_clades:
    current_clade_calls=apyp.nts2idx(np.array([x for x in get_internal_node_calls(clade,ancestral_reconstruction_fasta)]))
    has_mut_current_clade=(ancestral_reconstruction_nti[goodpos]!=current_clade_calls)
    modern_backbone_snps.update(set(np.where(has_mut_current_clade)[0]))
modern_backbone_snps_all_pos=goodpos[np.array(sorted(modern_backbone_snps))]
annotation_mutations_modern_true_backbone_relative_MRCA=apyp.annotate_mutations(annotation_genes , p[modern_backbone_snps_all_pos] , refnti_m[np.ix_(modern_backbone_snps_all_pos,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(modern_backbone_snps_all_pos,np.full((num_samples), True))] , calls[np.ix_(modern_backbone_snps_all_pos,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),modern_backbone_snps_all_pos)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(modern_backbone_snps_all_pos,np.full((num_samples), True))], mutQual[modern_backbone_snps_all_pos,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
#annotation_mutations_modern_true_backbone_relative_MRCA=annotation_mutations_modern_true_backbone_relative_MRCA[~annotation_mutations_modern_true_backbone_relative_MRCA['product'].isin(genes_to_drop)]


total_mutation_counts_for_pos=[]
for index,row in annotation_mutations.iterrows():
    row_chr_start=chrStarts[row['chr']-1]
    pos=row_chr_start+row['pos']
    index_in_hasmutation=np.where(p==pos-1)[0]
    total_mutation_counts_for_pos.append(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[index_in_hasmutation,:][:,~outgroup_bool]))
annotation_mutations.insert(27,'num_all_samples_mutated',total_mutation_counts_for_pos)


# modern internal
modern_internal = np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,:][:,np.isin(sampleNames,modern_samples)],axis=1)>1)[0]
annotation_mutations_modern_internals = apyp.annotate_mutations(annotation_genes , p[modern_internal] , refnti_m[np.ix_(modern_internal,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(modern_internal,np.full((num_samples), True))] , calls[np.ix_(modern_internal,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),modern_internal)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(modern_internal,np.full((num_samples), True))], mutQual[modern_internal,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
#annotation_mutations_modern_internals=annotation_mutations_modern_internals[~annotation_mutations_modern_internals['product'].isin(genes_to_drop)]

modern_singletons = np.where(np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,:][:,np.isin(sampleNames,modern_samples)],axis=1) == 1)[0]
annotation_mutations_modern_singletons = apyp.annotate_mutations(annotation_genes , p[modern_singletons] , refnti_m[np.ix_(modern_singletons,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(modern_singletons,np.full((num_samples), True))] , calls[np.ix_(modern_singletons,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),modern_singletons)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(modern_singletons,np.full((num_samples), True))], mutQual[modern_singletons,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
#annotation_mutations_modern_singletons=annotation_mutations_modern_singletons[~annotation_mutations_modern_singletons['product'].isin(genes_to_drop)]

# %% generating singletons check for blast (subsequent analysis in ancient_singleton_validation_blast.sh )
## singletons LNBA SNPs
candidate_lnba_singleton_snps_in_goodpos_relative_MRCA = np.where( (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,LNBA_clade_indices], axis=1) ==1) )[0]
candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark=  np.where( ((np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,np.isin(sampleNames,LNBA_samples[~np.isin(LNBA_samples,['ARK017'])])], axis=1) ) == 1 ) & ( (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,~LNBA_clade_indices],axis=1)==0)))[0]

annotation_mutations_lnba_singletons_relative_MRCA=apyp.annotate_mutations(annotation_genes , p[candidate_lnba_singleton_snps_in_goodpos_relative_MRCA] , refnti_m[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , calls[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),candidate_lnba_singleton_snps_in_goodpos_relative_MRCA)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))], mutQual[candidate_lnba_singleton_snps_in_goodpos_relative_MRCA,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
annotation_mutations_lnba_singletons_relative_MRCA_no_ark=apyp.annotate_mutations(annotation_genes , p[candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark] , refnti_m[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark,np.full((num_samples), True))] , calls[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark,np.full((num_samples), True))], mutQual[candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    

# all singletons
candidate_singleton_snps_in_goodpos_relative_MRCA=np.where( (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis, axis=1) ) == 1 )[0]
annotation_mutations_singletons_relative_MRCA=apyp.annotate_mutations(annotation_genes , p[candidate_singleton_snps_in_goodpos_relative_MRCA] , refnti_m[np.ix_(candidate_singleton_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(candidate_singleton_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , calls[np.ix_(candidate_singleton_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),candidate_singleton_snps_in_goodpos_relative_MRCA)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(candidate_singleton_snps_in_goodpos_relative_MRCA,np.full((num_samples), True))], mutQual[candidate_singleton_snps_in_goodpos_relative_MRCA,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    

# non-LNBA singletons
non_lnba_singletons = np.where((np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,~LNBA_clade_indices],axis=1)==1) & (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,LNBA_clade_indices],axis=1)==0))[0]
annotation_mutations_non_lnba_singletons = apyp.annotate_mutations(annotation_genes , p[modern_singletons] , refnti_m[np.ix_(non_lnba_singletons,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(non_lnba_singletons,np.full((num_samples), True))] , calls[np.ix_(non_lnba_singletons,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),non_lnba_singletons)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(non_lnba_singletons,np.full((num_samples), True))], mutQual[non_lnba_singletons,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    



# %% generate rates of mutation accumulation
# %% production figure 4d
def generate_proportion_of_mutation_type(annotation_mutation, type):
    subset=annotation_mutation.loc[~annotation_mutation.locustag.isna()]
    num_observed=len(subset.query(f'NonSyn == {type}'))
    return num_observed/len(annotation_mutation)

def generate_bootstrap_external_internal_rates(external,internal,type,num_bootstraps):
    output_bootstrapped_observations=np.zeros(num_bootstraps)
    for i in range(num_bootstraps):
        bootstrapped_external=np.random.randint(0,len(external),len(external))
        bootstrapped_internal=np.random.randint(0,len(internal),len(internal))
        observed_external=generate_proportion_of_mutation_type(external.iloc[bootstrapped_external], type)
        observed_internal=generate_proportion_of_mutation_type(internal.iloc[bootstrapped_internal], type)
        output_bootstrapped_observations[i]=(observed_external+1)/(observed_internal+1)
    return output_bootstrapped_observations
os.makedirs(f'pdf/mutation_rates/{analysis_params_output_name}', exist_ok=True)


boostrapped_rate_LNBA_nonsyn=generate_bootstrap_external_internal_rates(annotation_mutations_lnba_singletons_relative_MRCA,annotation_mutations_lnba_internal_calls_all_relative_MRCA,True,bootstraps)
boostrapped_rate_modern_nonsyn=generate_bootstrap_external_internal_rates(annotation_mutations_modern_singletons,annotation_mutations_modern_internals,True,bootstraps)


fig, ax = plt.subplots()
plt.hist(boostrapped_rate_LNBA_nonsyn,color='#2b8cbe',alpha=1,label='LNBA')
#plt.hist(boostrapped_rate_modern_syn,color='red',alpha=0.5,label='Modern Synonymous')
plt.hist(boostrapped_rate_modern_nonsyn,color='#d0d1e6',alpha=0.75,label='Extant')
plt.legend(loc='upper left')
ymax_new=ax.get_ylim()[1]*1.1
plt.ylim(0,ymax_new)
y_arrow_height=0.95*ymax_new
plt.axvline(x=1,ymin=0,ymax=y_arrow_height/ax.get_ylim()[1],color='black',linestyle='--',alpha=0.85)
plt.xlabel('Rate External/Rate Internal')
plt.ylabel('Number of bootstrap realizations')
plt.xlim(0.925,1.125)

# Arrow dimensions
x = 1
dx, dy = 0.005, 0  # dx = 0.25 to the right, dy = 0 means no vertical movement
# Draw a red arrow
ax.arrow(x, y_arrow_height, dx, dy, head_width=ax.get_ylim()[1]*0.025, head_length=0.005, fc='black', ec='black',alpha=0.85)

plt.savefig(f'pdf/mutation_rates/{analysis_params_output_name}/production_lnba_vs_modern.svg',bbox_inches='tight')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Para evo Analysis # # Para evo Analysis # # Para evo Analysis # # Para evo Analysis # # Para evo Analysis # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# %% Parallel evolution module
# =============================================================================
# Internal LNBA SNPs
# define dict with parameters 
# production figure S9b
for i in range(np.max(np.unique(annotation_mutations_lnba_internal_calls_all_relative_MRCA.locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'lnba_internal_min_muts_'+str(parameters['Min_num_mutations_cand'])
    ## only SNPs that appear and spread down LNBA lineage (including those shared by all after MRCA of all LNBA + Modern)

    [res_cand_nummut,annotation_mutation_paraSignal_lnba_internal] = apyp.parallel_evo_module_production_pestis( candidate_backbone_snps_in_goodpos_relative_MRCA, contig_positions , annotation_mutations_lnba_internal_calls_all_relative_MRCA , annotation_genes , parameters, True, False,"Parallel evolution: Internal LNBA SNVs")
    #candidate_genes_lnba_backbone=annotation_mutation_paraSignal_lnba_internal[['chr','pos','type','product','gene','locustag','num_LNBA_samples_mutated','num_all_samples_mutated','num_non_LNBA_samples_mutated']]
    #sig_annotation_mutation_paraSignal_lnba_internal=annotation_mutation_paraSignal_lnba_internal[annotation_mutation_paraSignal_lnba_internal['chr_locustag'].isin(res_cand_nummut[np.where(res_cand_nummut[:,5]<0.05)[0]][:,0])]

for i in range(np.max(np.unique(annotation_mutations_lnba_internal_calls_all_relative_MRCA.query('type == "N"').locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'lnba_internal_nonsyn_min_muts_'+str(parameters['Min_num_mutations_cand'])
    nonsyn_singletons=np.where(np.in1d(p,annotation_mutations_lnba_internal_calls_all_relative_MRCA.query('type == "N"').p))[0]

    apyp.parallel_evo_module_production_pestis( nonsyn_singletons, contig_positions , annotation_mutations_lnba_internal_calls_all_relative_MRCA.query('type == "N"').reset_index() , annotation_genes , parameters, True, False,"Multiply mutated genes in LNBA-internal nonsynonymous positions")

# %% 
# all LNBA SINGLETONS 
for i in range(np.max(np.unique(annotation_mutations_lnba_singletons_relative_MRCA.locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'lnba_singletons_all_min_muts_'+str(parameters['Min_num_mutations_cand'])

    ## only LNBA singleton SNPs
    [res_cand_nummut_LNBA_singletons,annotation_mutation_paraSignal_LNBA_singletons] = apyp.parallel_evo_module_production_pestis( candidate_lnba_singleton_snps_in_goodpos_relative_MRCA, contig_positions , annotation_mutations_lnba_singletons_relative_MRCA , annotation_genes , parameters, True, False,"Multiply mutated genes in LNBA-singletons positions")
    sig_annotation_mutation_paraSignal_LNBA_singletons=annotation_mutation_paraSignal_LNBA_singletons[annotation_mutation_paraSignal_LNBA_singletons['chr_locustag'].isin(res_cand_nummut_LNBA_singletons[np.where(res_cand_nummut_LNBA_singletons[:,5]<0.05)[0]][:,0])]

# # # # 
    # subsetting to only nonysn
for i in range(np.max(np.unique(annotation_mutations_lnba_singletons_relative_MRCA.query('type == "N"').locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'TEST_lnba_singletons_nonsyn_all_min_muts_'+str(parameters['Min_num_mutations_cand'])
    nonsyn_singletons=np.where(np.in1d(p,annotation_mutations_lnba_singletons_relative_MRCA.query('type == "N"').p))[0]
    apyp.parallel_evo_module_production_pestis( nonsyn_singletons, contig_positions , annotation_mutations_lnba_singletons_relative_MRCA.query('type == "N"').reset_index() , annotation_genes , parameters, True, False,"Multiply mutated genes in LNBA-singletons nonsynonymous positions")

# Production figure 4e (no title)
for i in range(np.max(np.unique(annotation_mutations_lnba_singletons_relative_MRCA.query('type == "N"').locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'production_lnba_singletons_nonsyn_all_min_muts_'+str(parameters['Min_num_mutations_cand'])
    nonsyn_singletons=np.where(np.in1d(p,annotation_mutations_lnba_singletons_relative_MRCA.query('type == "N"').p))[0]
    para_evo_out1_lnba_singletons_nonsyn,para_evo_out2_lnba_singletons_nonsyn=apyp.parallel_evo_module_production_pestis( nonsyn_singletons, contig_positions , annotation_mutations_lnba_singletons_relative_MRCA.query('type == "N"').reset_index() , annotation_genes , parameters, True, False,"")

# manually get this from the file with 
    # tail -n23 parallel_evo_results/strand_1_maf_9_heterozygosityAncient_recombinationAll_aidaexclusionsRepeatRNA_minorafcovper95/parallel_evo_module_results_lnba_singletons_nonsyn_all_min_muts_2.txt | cut -f1    
output_accession_order=['WP_002212253.1',
'WP_002214644.1',
'WP_002210477.1',
'WP_002228103.1',
'WP_002210389.1',
'WP_002357277.1',
'WP_002211388.1',
'WP_002224905.1',
'WP_002211638.1',
'WP_002210726.1',
'WP_002214197.1',
'WP_002211888.1',
'WP_002218214.1',
'WP_002210878.1',
'WP_002212060.1',
'WP_002211250.1',
'WP_002227945.1',
'WP_002210215.1',
'WP_002212209.1',
'WP_002229547.1',
'WP_002209276.1',
'WP_002210160.1',
'WP_002211776.1']

for locustag in output_accession_order:
    #print(locustag)
    locustag_to_search=f'cds-{locustag}'
    muts_to_output=[]
    for position in para_evo_out2_lnba_singletons_nonsyn.query(f'locustag == "{locustag_to_search}"').p:
        sample_with_mut=sampleNames[hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[p==position][0]][0]
        mut_to_output=para_evo_out2_lnba_singletons_nonsyn.query(f'p == {position}').iloc[0].muts[0]
        muts_to_output.append(f'{sample_with_mut}:{mut_to_output}')
    print(','.join(muts_to_output))

# get YPO names with:
    #tail -n23 parallel_evo_results/strand_1_maf_9_heterozygosityAncient_recombinationAll_aidaexclusionsRepeatRNA_minorafcovper95/parallel_evo_module_results_lnba_singletons_nonsyn_all_min_muts_2.txt | cut -f1 | while read accession
    #do
    #grep "${accession}" -B1 /Users/ad_loris/Nextcloud/keylab/reference_genomes/Ypestis_ASM906v1/GCF_000009065.1_ASM906v1_genomic_gffbackup | head -n1 | awk -F'old_locus_tag=' '{print $2}'
    #done

# get chrom_start_stop with
    #tail -n23 parallel_evo_results/strand_1_maf_9_heterozygosityAncient_recombinationAll_aidaexclusionsRepeatRNA_minorafcovper95/parallel_evo_module_results_lnba_singletons_nonsyn_all_min_muts_2.txt | cut -f1 | while read accession
    #do
    #grep "${accession}" -B1 /Users/ad_loris/Nextcloud/keylab/reference_genomes/Ypestis_ASM906v1/GCF_000009065.1_ASM906v1_genomic_gffbackup | head -n1 | awk -F'\t' '{print $1":"$4"-"$5}'
    #done
# %%
# production figure S9a
only_human_modern_samples=np.in1d(sampleNames,samples_to_isolation['sampleName'][samples_to_isolation['Host/vector']=='Patient'])
has_mut_modern_human_singletons=np.where(np.sum(hasmutation_relative_ancestral_MRCA_all_pestis[:,:][:,only_human_modern_samples],axis=1)==1)[0]
has_no_mut_non_modern_human=np.where(np.sum(hasmutation_relative_ancestral_MRCA_all_pestis[:,:][:,~only_human_modern_samples],axis=1)==0)[0]
exclusive_modern_human_singletons=np.intersect1d(has_mut_modern_human_singletons,has_no_mut_non_modern_human)
annotation_mutations_exclusive_modern_humans_singletons=apyp.annotate_mutations(annotation_genes , p[exclusive_modern_human_singletons] , refnti_m[np.ix_(exclusive_modern_human_singletons,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(exclusive_modern_human_singletons,np.full((num_samples), True))] , calls[np.ix_(exclusive_modern_human_singletons,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),exclusive_modern_human_singletons)] , hasmutation_relative_ancestral_MRCA_all_pestis[np.ix_(exclusive_modern_human_singletons,np.full((num_samples), True))], mutQual[exclusive_modern_human_singletons,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    


for i in range(2,3):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':'0',
            'output_name':'modern_singletons',
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1",
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name']='modern_human_min_muts_'+str(parameters['Min_num_mutations_cand'])

    [res_cand_nummut_modern_human_singletons,annotation_mutation_paraSignal_modern_human_singletons] = apyp.parallel_evo_module_production_pestis( exclusive_modern_human_singletons, contig_positions , annotation_mutations_exclusive_modern_humans_singletons , annotation_genes, parameters, True, False,"Multiply mutated genes in modern singletons")
# # # # 
    # subsetting to only nonsyn
for i in range(2,3):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':'0',
            'output_name':'modern_singletons',
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1",
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name']='modern_human_nonsyn_min_muts_'+str(parameters['Min_num_mutations_cand'])
    nonsyn_singletons=np.where(np.in1d(p,annotation_mutations_exclusive_modern_humans_singletons.query('type == "N"').p))[0]

    apyp.parallel_evo_module_production_pestis( nonsyn_singletons, contig_positions , annotation_mutations_exclusive_modern_humans_singletons.query('type == "N"').reset_index() , annotation_genes , parameters, True, False,"Multiply mutated genes in extant human infections-singleton nonsynonymous positions")

# %%
# LNBA SINGLETONS NO ARK
for i in range(np.max(np.unique(annotation_mutations_lnba_singletons_relative_MRCA_no_ark.locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'lnba_singletons_no_ark_all_min_muts_'+str(parameters['Min_num_mutations_cand'])

    ## only LNBA singleton SNPs
    [res_cand_nummut_LNBA_no_ark_singletons,annotation_mutation_paraSignal_LNBA_no_ark_singletons] = apyp.parallel_evo_module_production_pestis( candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark, contig_positions , annotation_mutations_lnba_singletons_relative_MRCA_no_ark , annotation_genes , parameters, True, False,"Multiply mutated genes in LNBA-singletons positions")
    sig_annotation_mutation_paraSignal_LNBA_no_ark_singletons=annotation_mutation_paraSignal_LNBA_no_ark_singletons[annotation_mutation_paraSignal_LNBA_no_ark_singletons['chr_locustag'].isin(res_cand_nummut_LNBA_no_ark_singletons[np.where(res_cand_nummut_LNBA_no_ark_singletons[:,5]<0.05)[0]][:,0])]

for i in range(np.max(np.unique(annotation_mutations_lnba_singletons_relative_MRCA_no_ark.query('type == "N"').locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'lnba_singletons_no_ark_nonsyn_all_min_muts_'+str(parameters['Min_num_mutations_cand'])
    nonsyn_singletons=np.where(np.in1d(p,annotation_mutations_lnba_singletons_relative_MRCA_no_ark.query('type == "N"').p))[0]
    apyp.parallel_evo_module_production_pestis( nonsyn_singletons, contig_positions , annotation_mutations_lnba_singletons_relative_MRCA_no_ark.query('type == "N"').reset_index() , annotation_genes , parameters, True, False,"Multiply mutated genes in LNBA-singletons positions")


# %%
# LNBA SINGLETONS NO KZL002 (last sample)
candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl=  np.where( ((np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,np.isin(sampleNames,LNBA_samples[~np.isin(LNBA_samples,['KZL002'])])], axis=1) ) == 1 ) & ( (np.sum(hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[:,~np.isin(sampleNames,LNBA_samples[~np.isin(LNBA_samples,['KZL002'])])],axis=1)==0)))[0]
annotation_mutations_lnba_singletons_relative_MRCA_no_kzl=apyp.annotate_mutations(annotation_genes , p[candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl] , refnti_m[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl,np.full((num_samples), True))] , ancestral_reconstruction_nti_m[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl,np.full((num_samples), True))] , calls[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl,np.full((num_samples), True))] , counts[np.ix_(np.full((num_samples), True),np.full((8), True),candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl)] , hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[np.ix_(candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl,np.full((num_samples), True))], mutQual[candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    

for i in range(np.max(np.unique(annotation_mutations_lnba_singletons_relative_MRCA_no_kzl.locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'lnba_singletons_no_kzl_all_min_muts_'+str(parameters['Min_num_mutations_cand'])

    ## only LNBA singleton SNPs
    [res_cand_nummut_LNBA_no_ark_singletons,annotation_mutation_paraSignal_LNBA_no_kzl_singletons] = apyp.parallel_evo_module_production_pestis( candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_kzl, contig_positions , annotation_mutations_lnba_singletons_relative_MRCA_no_kzl , annotation_genes , parameters, True, False,"Multiply mutated genes in LNBA-singletons positions")
for i in range(np.max(np.unique(annotation_mutations_lnba_singletons_relative_MRCA_no_kzl.query('type == "N"').locustag.dropna(),return_counts=True)[1]),1,-1):
    parameters = {
            'NumTrialsSim':10000,
            'Min_num_mutations_cand':i, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':0, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum': f"{analysis_params_output_name}_mutationalspectrum.py.pk1", # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 50, # within subject 15 is fine
            'timestamp':timestamp,
            'analsysis_params_output_name_folder': analysis_params_output_name
            }
    parameters['output_name'] = 'lnba_singletons_no_kzl_nonsyn_all_min_muts_'+str(parameters['Min_num_mutations_cand'])
    nonsyn_singletons=np.where(np.in1d(p,annotation_mutations_lnba_singletons_relative_MRCA_no_kzl.query('type == "N"').p))[0]
    apyp.parallel_evo_module_production_pestis( nonsyn_singletons, contig_positions , annotation_mutations_lnba_singletons_relative_MRCA_no_kzl.query('type == "N"').reset_index() , annotation_genes , parameters, True, False,"Multiply mutated genes in LNBA-singletons positions")

print('proportion of KZL singletons of all LNBA singletons:', (len(annotation_mutations_lnba_singletons_relative_MRCA)-len(annotation_mutations_lnba_singletons_relative_MRCA_no_kzl))/len(annotation_mutations_lnba_singletons_relative_MRCA))
 
# %% 
# validation of KZL002 non outlier status for number of double hits

def simulate_kzl_double_draws_distribution(para_evo_annotation_mutation_set,singleton_set):
    num_double_picks=0
    num_genes_to_pick=len(np.unique(para_evo_annotation_mutation_set.locustag))
    for locus in np.unique(para_evo_annotation_mutation_set.locustag):
        query_result=para_evo_annotation_mutation_set.query(f'locustag == "{locus}"')
        picks=[]
        for position_in_locus in query_result.p:
            sample_with_mutaiton=sampleNames[hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[p == position_in_locus][0]][0]
            print(f'position: {position_in_locus}, sample: {sample_with_mutaiton}')
            picks.append(sample_with_mutaiton)
        if picks == ['KZL002']*2:
            num_double_picks+=1

    # calculate proprotion of mutations from KZL
    non_kzl=0
    kzl=0
    for s_idx in np.where(LNBA_clade_indices)[0]:
        if sampleNames[s_idx] != 'KZL002':
            non_kzl+=hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[singleton_set].sum(axis=0)[s_idx]
        else:
            kzl+=hasmutation_relative_ancestral_MRCA_LNBA_modern_pestis[singleton_set].sum(axis=0)[s_idx]

    # run simulations of how likely it is to get two draws from KZL, given the amount of mutations it has
    trial_results=[]
    for num_trials in range(10000):
        num_both_draws_from_kzl=0
        for x in range(num_genes_to_pick): #number of genes with multiple mutations that we need to find mutations for
            random_chooser=np.random.default_rng()
            draw1,draw2 = random_chooser.choice(kzl+non_kzl, size=2, replace=False)
            if draw1 < kzl and draw2 < kzl: #164 total SNVs are available from KZL, set these as first 163 draws
                num_both_draws_from_kzl+=1
        trial_results.append(num_both_draws_from_kzl)
    print('Percentile of observed double picks in KZL002 vs simulated:',stats.percentileofscore(trial_results,num_double_picks))

simulate_kzl_double_draws_distribution(annotation_mutation_paraSignal_LNBA_no_ark_singletons,candidate_lnba_singleton_snps_in_goodpos_relative_MRCA_no_ark)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# dN/dS Analysis # # dN/dS Analysis # # dN/dS Analysis # # dN/dS Analysis # # dN/dS Analysis # # dN/dS Analysis # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# %% 
# data gathering/parsing functions

# helper functions
def convert_dn_ds_prop_to_value(output_from_calc_dnds_across_genes_with_CI,ds_dn=False):
    total_muts=output_from_calc_dnds_across_genes_with_CI[5]
    observed_pN=output_from_calc_dnds_across_genes_with_CI[0]
    if ds_dn:
        expectation_pS = 1 - output_from_calc_dnds_across_genes_with_CI[4]
        expected_denom=(expectation_pS*total_muts)/((1-expectation_pS)*total_muts)
        observed_pS=1-observed_pN
        observed_stat=(observed_pS*total_muts)/((1-observed_pS)*total_muts)
        low_pN=output_from_calc_dnds_across_genes_with_CI[2][0]
        high_pN=output_from_calc_dnds_across_genes_with_CI[2][1]
        CI_low,CI_high=((1-low_pN)*total_muts)/(low_pN*total_muts),((1-high_pN)*total_muts)/(high_pN*total_muts)
    else:
        expectation_pN=output_from_calc_dnds_across_genes_with_CI[4]
        expected_denom=(expectation_pN*total_muts)/((1-expectation_pN)*total_muts)
        observed_stat=(observed_pN*total_muts)/((1-observed_pN)*total_muts)
        low_pN=output_from_calc_dnds_across_genes_with_CI[2][0]
        high_pN=output_from_calc_dnds_across_genes_with_CI[2][1]
        CI_low,CI_high=(low_pN*total_muts)/((1-low_pN)*total_muts),(high_pN*total_muts)/((1-high_pN)*total_muts)
    return (observed_stat/expected_denom,CI_low/expected_denom,CI_high/expected_denom)

def parse_data_dict(data_collector,dn_ds=True, ds_dn=False):
    data_dict_out={'value':[],'high':[],'low':[]}
    if dn_ds and not ds_dn:
        estimate_df=[1 for x in data_collector]
        for value_list in data_collector:
            if value_list[0]==1:
                data_dict_out['value'].append(value_list[0])
                data_dict_out['high'].append(value_list[0])
                data_dict_out['low'].append(value_list[0])
            elif type(value_list[2]) == tuple:
                data_dict_out['value'].append(convert_dn_ds_prop_to_value(value_list)[0])
                data_dict_out['high'].append(convert_dn_ds_prop_to_value(value_list)[2])
                data_dict_out['low'].append(convert_dn_ds_prop_to_value(value_list)[1])
            else:
                data_dict_out['value'].append(value_list[0])
                data_dict_out['high'].append(value_list[0])
                data_dict_out['low'].append(value_list[0])
    elif dn_ds and ds_dn:
        estimate_df=[1 for x in data_collector]
        for value_list in data_collector:
            if type(value_list[2]) == tuple:
                data_dict_out['value'].append(convert_dn_ds_prop_to_value(value_list,True)[0])
                data_dict_out['high'].append(convert_dn_ds_prop_to_value(value_list,True)[2])
                data_dict_out['low'].append(convert_dn_ds_prop_to_value(value_list,True)[1])
            else:
                data_dict_out['value'].append(value_list[0])
                data_dict_out['high'].append(value_list[0])
                data_dict_out['low'].append(value_list[0])
    else:
        for value_list in data_collector:
            if type(value_list[2]) != tuple: ## na returned for CI, therefore NAN value
                data_dict_out['value'].append(value_list[0])
                data_dict_out['high'].append(value_list[0])
                data_dict_out['low'].append(value_list[0])
            else:
                data_dict_out['value'].append(value_list[0])
                data_dict_out['high'].append(value_list[2][1])
                data_dict_out['low'].append(value_list[2][0])
        estimate_df=[x[4] for x in data_collector]
    return data_dict_out,estimate_df

def generate_dn_ds_stats_across_sets(group_collector,sets_to_test,substitution_spectrum,annotation_genes):
    data_collector=[]
    for index,this_set_mutation_annotations in enumerate(sets_to_test):
        params={'subjectID':group_collector[index],'substitution_spectrum': substitution_spectrum}
        data_collector.append(apyp.calculate_dn_ds_across_genes(params,annotation_genes,this_set_mutation_annotations,True))
    return data_collector

# %%
# production run dN/dS ratio figure 4c
os.makedirs(os.getcwd() + f'/pdf/dnds/{analysis_params_output_name}', exist_ok=True)
os.makedirs(os.getcwd() + f'/pdf/adaptive_evo/{analysis_params_output_name}', exist_ok=True)

group_collector=['LNBA','Extant']
# all
sets_to_test=[annotation_mutations_exclusive_lnba,annotation_mutations_modern]

GROUP_TO_GET=[x for x in range(len(group_collector))]
group_collector_final=[]
sets_to_test_final=[]
for x in GROUP_TO_GET:
    group_collector_final.append(group_collector[x])
    sets_to_test_final.append(sets_to_test[x])

# collect + plot
dn_ds_all_sets=generate_dn_ds_stats_across_sets(group_collector_final,sets_to_test_final,f"{analysis_params_output_name}_mutationalspectrum.py.pk1",annotation_genes)

data_observed_df,estimate_df=parse_data_dict(data_collector=dn_ds_all_sets,dn_ds=True)
df=pd.DataFrame.from_dict(data_observed_df,orient='index',columns=group_collector_final)


fig,axs=plt.subplots(facecolor='white')
plt.axhline(y = 1,xmin=0,xmax=0.7, color = 'black', linestyle = (0, (5, 2.5)),linewidth=2,alpha=0.55)
axs.set_xlim(0.15,2.25)
axs.set_ylim(0.5,1.5)

#sns.pointplot(y=df.loc['value'],x=group_collector_final,join=False,markers='D')
lnba_low=df.loc['low']['LNBA']
lnba_high=df.loc['high']['LNBA']
extant_low=df.loc['low']['Extant']
extant_high=df.loc['high']['Extant']
central_location_lnba_error_bar=(0.25,df.loc['value']['LNBA']-0.5*(lnba_high-lnba_low))
central_location_extant_error_bar=(1,df.loc['value']['Extant']-0.5*(extant_high-extant_low))

axs.add_patch(Rectangle(central_location_lnba_error_bar, 0.5, (lnba_high-lnba_low),alpha=0.75,color='#2b8cbe'))
axs.add_patch(Rectangle(central_location_extant_error_bar, 0.5, (extant_high-extant_low),alpha=0.75,color='#d0d1e6'))
plt.scatter(x=0.5,y=df.loc['value']['LNBA'],marker='_',c='black',s=6500,linewidths=2.5)
plt.scatter(x=1.25,y=df.loc['value']['Extant'],marker='_',c='black',s=6500,linewidths=2.5)
axs.set_xticks([0.5,1.25],['LNBA','Extant'])
plt.ylabel('dN/dS Ratio')

plt.savefig(f'pdf/dnds/{analysis_params_output_name}/production_dnds_lnba_vs_extant.svg',bbox_inches="tight")


# %% 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Pairwise Comparisons analysis # # Pairwise Comparisons analysis # # Pairwise Comparisons analysis # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# %% 
# helper functions
def calculate_se_slope(X,y,residuals,return_value=False):
    # Sum of squared residuals
    SSR = np.sum(residuals**2)

    # Number of observations and features
    n = len(y)
    k = X.shape[1]

    # Mean of X
    mean_X = np.mean(X)

    # Sum of the squared differences between X and its mean
    sum_sq_diff = np.sum((X - mean_X)**2)

    # Standard error of the coefficient
    SE_beta = np.sqrt(SSR / (n - 2) / sum_sq_diff)
    if return_value:
        return SE_beta
    print("Standard Error of the Coefficient:", SE_beta)

def generate_upper_and_lower_age_differences(row1,row2):
    min = row1['End_Range_C14'] - row2['Start_Range_C14']
    max = row1['Start_Range_C14'] - row2['End_Range_C14']
    return min,max

def get_pairwise_genetic_distance(sample1,sample2,calls,sampleNames):
    col_index_to_compare=np.nonzero(np.in1d(sampleNames, sample1))[0]
    row_index_to_compare=np.nonzero(np.in1d(sampleNames, sample2))[0]
    calls_same_for_comparissons=np.sum(((calls[:,col_index_to_compare] == calls[:,row_index_to_compare]) & (calls[:,col_index_to_compare] != 4) & (calls[:,row_index_to_compare] != 4)))
    denominator_for_comparissons=np.sum(((calls[:,col_index_to_compare] !=4 ) & (calls[:,row_index_to_compare] != 4)))
    normalized_genetic_distance=(denominator_for_comparissons-calls_same_for_comparissons)/denominator_for_comparissons
    return normalized_genetic_distance

# %%
# generating min distance travelled per day ranges scatter plot for next nearest contemporaneous samples
pairwise_geographic_distances=np.zeros((len(lnba_metadata.query('Branch == "LNBA"')[['Sample_Name', 'Start_Range_C14', 'End_Range_C14']]),len(lnba_metadata.query('Branch == "LNBA"')[['Sample_Name', 'Start_Range_C14', 'End_Range_C14']])))


for r_index,row in lnba_metadata.query('Branch == "LNBA"').reset_index().iterrows():
    for r_index2,row2 in lnba_metadata.query('Branch == "LNBA"').reset_index().iterrows():
        header_lat=float(row['Latitude'])
        header_lon=float(row['Longitude'])
        row_lat=float(row2['Latitude'])
        row_lon=float(row2['Longitude'])
        pairwise_geographic_distances[r_index,r_index2]=apyp.measure_dist_from_coords(header_lat,header_lon,row_lat,row_lon)

pairwise_genetic_distances=np.zeros((np.sum(LNBA_clade_indices),np.sum(LNBA_clade_indices)))

# calculating if any differences in estimated slope?


sampleNames_in_metadata=[s for s in sampleNames[LNBA_clade_indices] if s in list(lnba_metadata.Sample_Name)]
sampleNames_in_metadata=sorted(sampleNames_in_metadata,key=lambda x: tree_sample_order.index(x))

sampleNames_in_metadata_no_ark=[x for x in sampleNames_in_metadata if x != 'ARK017']

# output data for itol heatmap for C14 dates
investigated=[]
for i,x in lnba_metadata.query('Branch == "LNBA"').query('Dating == "C14"').iterrows():
    investigated.append(x.Sample_Name)
    print(str(x.Sample_Name)+','+str(x.median_bce_date_calibrated))
for x in sampleNames:
    if x not in investigated:
        print(f'{x},X')

# %% Plotting categorical variable ranges
plt.subplots(facecolor='white',figsize=(20,15))
lnba_metadata_only_sorted=lnba_metadata.query('Branch == "LNBA" & Dating =="C14"').sort_values('Median_C14',ascending=False).reset_index()
plt.scatter(lnba_metadata_only_sorted['Sample_Name'],lnba_metadata_only_sorted['Median_C14'],color='black',marker='d')
error_above=lnba_metadata_only_sorted['Start_Range_C14']-lnba_metadata_only_sorted['Median_C14']
error_below=lnba_metadata_only_sorted['Median_C14']-lnba_metadata_only_sorted['End_Range_C14']
plt.errorbar(lnba_metadata_only_sorted['Sample_Name'],lnba_metadata_only_sorted['Median_C14'],yerr=(error_above,error_below),fmt='none', elinewidth=30,alpha=0.5)
plt.tick_params(axis='x',labelrotation=90)

# %% generating distribution of tree distances from all LNBA pairwise comparisons:
pairwise_tree_distances_non_ark=[]
pairwise_tree_distances_ark=[]
parsed=[]
for s in LNBA_samples:
    parsed.append(s)
    other_samples = [x for x in LNBA_samples if x not in parsed]
    for s2 in other_samples:
        if s == 'ARK017' or s2 == 'ARK017':
            pairwise_tree_distances_ark.append(tree.distance(s,s2))
        else:
            pairwise_tree_distances_non_ark.append(tree.distance(s,s2))

data_non_ark = pd.DataFrame(pairwise_tree_distances_non_ark, columns=['Distance'])
data_non_ark['Type'] = 'Human-Human'

data_ark = pd.DataFrame(pairwise_tree_distances_ark, columns=['Distance'])
data_ark['Type'] = 'Sheep-Human'

# Combine the data
data_combined = pd.concat([data_non_ark, data_ark])

# Step 2: Plotting
# Create the jitter plot
fig, ax = plt.subplots(facecolor='white',figsize=(3,5))

sns.stripplot(x='Type', y='Distance', data=data_combined, jitter=True,alpha=0.5,hue='Type',palette=['#2b8cbe', '#c90076'],legend=False)

# calculate pvalue:
ark_non_ark_pairwise_out=stats.mannwhitneyu(data_non_ark['Distance'],data_ark['Distance'])
ark_non_ark_pairwise_out_pvalue=ark_non_ark_pairwise_out.pvalue
if ark_non_ark_pairwise_out_pvalue > 0.05:
    ark_non_ark_pairwise_out_pvalue='ns'
x1,x2=0,1
y=np.max(data_combined.Distance)+0.25*np.max(data_combined.Distance)
h=y*.05
ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.25, c='black')
ax.text((x1+x2)*0.5,y+h, f'{ark_non_ark_pairwise_out_pvalue}', ha='center', va='bottom', color='black')

ax.set_ylim(0,0.3)
# Optional: Add titles and labels
plt.xlabel('Pairwise Comparisons')
plt.ylabel('Genetic Distance')
plt.savefig('pdf/pairwise_tree_distances_ark_vs_non_ark.svg',bbox_inches='tight')

# %% generating root to tip
root_to_tip_data=np.zeros((len(sampleNames_in_metadata),2))
for index,sample1 in enumerate(sampleNames_in_metadata):
    dist_to_root = tree.distance(sample1,MRCA_pestis_lnba)
    sample1_metadata=lnba_metadata.query('Sample_Name == @sample1').iloc[0]
    temporal_signal=sample1_metadata['median_bce_date_calibrated']
    root_to_tip_data[index,0]=dist_to_root
    root_to_tip_data[index,1]=-float(temporal_signal)

plot = sns.regplot(x=root_to_tip_data[:,0],y=root_to_tip_data[:,1])
plt.title('Root-to-tip from MRCA LNBA')
plt.xlabel('Genetic distance')
plt.ylabel('C14 date estimate')

slope, intercept, r_value, p_value, std_err = stats.linregress(x=root_to_tip_data[:,0],y=root_to_tip_data[:,1])
plt.savefig('pdf/root_to_tip_lnba.png')
plt.close()

# generating root to tip, highlighting ark
root_to_tip_data=np.zeros((len(sampleNames_in_metadata),2))
ark_root_to_tip=np.zeros((1,2))
for index,sample1 in enumerate(sampleNames_in_metadata):
    dist_to_root = tree.distance(sample1,MRCA_pestis_lnba)
    sample1_metadata=lnba_metadata.query('Sample_Name == @sample1').iloc[0]
    temporal_signal=sample1_metadata['median_bce_date_calibrated']
    root_to_tip_data[index,0]=dist_to_root
    root_to_tip_data[index,1]=-float(temporal_signal)
    if sample1 == 'ARK017':
        dist_to_root = tree.distance(sample1,MRCA_pestis_lnba)
        sample1_metadata=lnba_metadata.query('Sample_Name == @sample1').iloc[0]
        temporal_signal=sample1_metadata['median_bce_date_calibrated']
        ark_root_to_tip[0,0]=dist_to_root
        ark_root_to_tip[0,1]=-float(temporal_signal)
plot = sns.regplot(x=root_to_tip_data[:,0],y=root_to_tip_data[:,1])
plt.title('Root-to-tip from MRCA LNBA, ARK=red')
plt.xlabel('Genetic distance')
plt.ylabel('C14 date estimate')
slope, intercept, r_value, p_value, std_err = stats.linregress(x=root_to_tip_data[:,0],y=root_to_tip_data[:,1])
plt.scatter(x=ark_root_to_tip[0,0],y=ark_root_to_tip[0,1],c='r')
plt.savefig('pdf/root_to_tip_lnba_ark_highlighted.png')
plt.close()

# %%
# tree dist vs c14 difference linear regression
calculated_c14_diff=[]
calculated_genetic_dist=[]
calculated_geo_dist=[]
median_c14_based_rate=[]
parsed=[]
comparisons=[]
for s in sampleNames_in_metadata_no_ark:
    parsed.append(s)
    other_samples = [x for x in sampleNames_in_metadata_no_ark if x not in parsed]
    for s2 in other_samples:
        sample1_metadata=lnba_metadata.query('Sample_Name == @s').iloc[0]
        sample2_metadata=lnba_metadata.query('Sample_Name == @s2').iloc[0]

        calculated_c14_diff.append(abs(sample1_metadata['median_bce_date_calibrated'] - sample2_metadata['median_bce_date_calibrated']))
        calculated_genetic_dist.append(tree.distance(s,s2))
        calculated_geo_dist.append(apyp.measure_dist_from_coords(sample1_metadata['Latitude'],sample1_metadata['Longitude'],sample2_metadata['Latitude'],sample2_metadata['Longitude'])/1000)
        median_c14_based_rate.append((apyp.measure_dist_from_coords(sample1_metadata['Latitude'],sample1_metadata['Longitude'],sample2_metadata['Latitude'],sample2_metadata['Longitude'])/1000)/abs(sample1_metadata['median_bce_date_calibrated'] - sample2_metadata['median_bce_date_calibrated']))
        comparisons.append((s,s2))

calculated_c14_diff_ark=[]
calculated_genetic_dist_ark=[]
calculated_geo_dist_ark=[]
median_c14_based_rate_ark=[]
comparisons_ark=[]
for s in ['ARK017']:
    other_samples = [x for x in sampleNames_in_metadata if x != s]
    for s2 in other_samples:
        sample1_metadata=lnba_metadata.query('Sample_Name == @s').iloc[0]
        sample2_metadata=lnba_metadata.query('Sample_Name == @s2').iloc[0]

        calculated_c14_diff_ark.append(abs(sample1_metadata['median_bce_date_calibrated'] - sample2_metadata['median_bce_date_calibrated']))
        calculated_genetic_dist_ark.append(tree.distance(s,s2))
        calculated_geo_dist_ark.append(apyp.measure_dist_from_coords(sample1_metadata['Latitude'],sample1_metadata['Longitude'],sample2_metadata['Latitude'],sample2_metadata['Longitude'])/1000)
        median_c14_based_rate_ark.append((apyp.measure_dist_from_coords(sample1_metadata['Latitude'],sample1_metadata['Longitude'],sample2_metadata['Latitude'],sample2_metadata['Longitude'])/1000)/abs(sample1_metadata['median_bce_date_calibrated'] - sample2_metadata['median_bce_date_calibrated']))
        comparisons_ark.append((s,s2))
    
combined_calc_genetic_dist,combined_calc_c14_diff=calculated_genetic_dist+calculated_genetic_dist_ark,calculated_c14_diff+calculated_c14_diff_ark

# %%
# production figure 4b
## with intercept and p-value calc:
import statsmodels.api as sm

X,y=np.array(calculated_genetic_dist).reshape(-1, 1) ,np.array(calculated_c14_diff).reshape(-1, 1) 

X2 = sm.add_constant(X)
est = sm.OLS(y, X2)
est2 = est.fit()
print(est2.summary())

intercept,slope=est2.params

X_ark,y=np.array(calculated_genetic_dist_ark).reshape(-1, 1) ,np.array(calculated_c14_diff_ark).reshape(-1, 1) 

X2 = sm.add_constant(X_ark)
est_ark = sm.OLS(y, X2)
est2_ark = est_ark.fit()
print(est2_ark.summary())

intercept_ark,slope_ark=est2_ark.params

def calculate_se_slope(X,y,residuals,return_value=False):
    # Sum of squared residuals
    SSR = np.sum(residuals**2)

    # Number of observations and features
    n = len(y)
    k = X.shape[1]

    # Mean of X
    mean_X = np.mean(X)

    # Sum of the squared differences between X and its mean
    sum_sq_diff = np.sum((X - mean_X)**2)

    # Standard error of the coefficient
    SE_beta = np.sqrt(SSR / (n - 2) / sum_sq_diff)
    if return_value:
        return SE_beta
    print("Standard Error of the Coefficient:", SE_beta)

residuals = np.array(calculated_c14_diff)- (np.array(calculated_genetic_dist)*slope)+intercept
std_dev_residuals =  stats.tstd(residuals)
std_err_reg_coef=calculate_se_slope(np.array(calculated_genetic_dist).reshape(-1, 1),np.array(calculated_c14_diff).reshape(-1, 1) ,residuals, return_value=True)

fig,ax = plt.subplots(facecolor='white',)
ax.fill_between([np.min(combined_calc_genetic_dist),np.max(combined_calc_genetic_dist)],np.array(([np.min(combined_calc_genetic_dist)*slope+intercept,np.max(combined_calc_genetic_dist)*slope+intercept])+2*std_dev_residuals),np.array(([np.min(combined_calc_genetic_dist)*slope+intercept,np.max(combined_calc_genetic_dist)*slope+intercept])-2*std_dev_residuals),alpha=0.2,color='gray')
ax.scatter(x=calculated_genetic_dist,y=calculated_c14_diff,color='#2b8cbe',alpha=0.25,label='Human-Human')
ax.scatter(x=calculated_genetic_dist_ark,y=calculated_c14_diff_ark,color='#c90076',label='Sheep-Human')
ax.plot([np.min(combined_calc_genetic_dist),np.max(combined_calc_genetic_dist)],[np.min(combined_calc_genetic_dist)*slope+intercept,np.max(combined_calc_genetic_dist)*slope+intercept],color='#2b8cbe')

ax.plot([np.min(combined_calc_genetic_dist),np.max(calculated_genetic_dist_ark)],[np.min(combined_calc_genetic_dist)*slope_ark+intercept_ark,np.max(calculated_genetic_dist_ark)*slope_ark+intercept_ark],color='#c90076')


plt.ylim(0,2500)
plt.xlim(0,0.21)
plt.xlabel('Genetic Distsance')
plt.ylabel('Median cal. BCE date difference')
plt.legend(loc='upper left')
plt.savefig('pdf/production_pairwise_dist_vs_c14_2SD_residuals_with_legend_both_regressions.svg',bbox_inches='tight')

# %%
# finding if ARK comparisons fall outside the comparison?
outlier_comparisons={x:0 for x in sampleNames_in_metadata}
for i,c14 in enumerate(np.array(calculated_c14_diff_ark)):
    if c14 < ((np.array(calculated_genetic_dist_ark)*slope)-2*std_dev_residuals)[i] or c14 > ((np.array(calculated_genetic_dist_ark)*slope)+2*std_dev_residuals)[i]:
        print(c14, (np.array(calculated_genetic_dist_ark)*slope-2*std_dev_residuals)[i],(np.array(calculated_genetic_dist_ark)*slope+2*std_dev_residuals)[i],comparisons_ark[i])
        for x in comparisons_ark[i]:
            outlier_comparisons[x]+=1


for i,c14 in enumerate(np.array(calculated_c14_diff)):
    if c14 < ((np.array(calculated_genetic_dist)*slope)-2*std_dev_residuals)[i] or c14 > ((np.array(calculated_genetic_dist)*slope)+2*std_dev_residuals)[i]:
        print(c14, (np.array(calculated_genetic_dist)*slope-2*std_dev_residuals)[i],(np.array(calculated_genetic_dist)*slope+2*std_dev_residuals)[i], np.array(comparisons)[i])
        for x in comparisons[i]:
            outlier_comparisons[x]+=1

# ranking of samples with number of outlier samples for estimated vs observed C14
list_of_tuples_outlier_comparisons = list(outlier_comparisons.items())
outlier_comparisons_sorted = sorted(list_of_tuples_outlier_comparisons, key=lambda x: x[1],reverse=True)
print(outlier_comparisons_sorted)





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# END Pairwise Comparisons analysis # # Pairwise Comparisons analysis # # END Pairwise Comparisons analysis # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# END ALL analysis # # END ALL analysis # # END ALL # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #