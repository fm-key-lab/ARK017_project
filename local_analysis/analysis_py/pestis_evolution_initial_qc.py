# %%
# =============================================================================
## Pestis evolution project, Ian Light
# =============================================================================
# SNP-based analyses - QC File:
# =============================================================================
# Parse and filter raw output from Candidate Mutation Table generation(Snakemake)
# # =============================================================================

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
from Bio import Phylo
import matplotlib.pyplot as plt
import math
from statsmodels.stats.multitest import multipletests
import datetime
import time
import upsetplot
import seaborn as sns
import gzip

## declare paths for finding analysispy_modules and reference genome folder

os.chdir('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/production_run')

ref_genome_folder='/Users/ad_loris/Nextcloud/keylab/reference_genomes/Ypestis_ASM906v1'

sys.path.append("/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/analysis_py")
import analysispy_modules_pestis as apyp

ancient_sample_names=['GEN72','GRS004','KLE031','KLE048','GRH001','VEL003','HGC068','GLZ001','GLZ002','KunilaII','Barcelona_3031','EDI001.A','ES_11972','ES_8291','1343UnTal85','6Post','ARS007','CHC004','Gok2','Gyvakarai1','HGC009','HOP001','HOP004','I2470','I5884','KLZ001','KNK001','KZL002','MIB054','RISE505','RISE509','RK1001','RT5','RT6','RV2039','VLI092','AE1175','BED024.A0102','BED028.A0102','BED030.A0102','BED034.A0102','Bolgar_2370','BRA001.A0101','Ellwangen_549_O','LAI009.A0101','MAN008.B0101','NMS002.A0101','OBS107','OBS110','OBS116','OBS124','OBS137','STA001.A0101','STN002.A0101','STN007.A0101','STN008.A0101','STN014.A0101','STN019.A0101','STN020.A0101','STN021.A0101','OOH003','ARK017', 'DSH008', 'DSH025','C10091', 'C10098', 'C10928']
# %% Filtering parameters
# =============================================================================

# Much of this will need to vary according to your reference genome,
# coverage, and particular samples
# consider to add SD filter in addition to coverage. Isolates of other species that harbor plasmid present in pangenome at very high cov can otherwise meet the coverage threshold even though FP

# variable for outputting files names:

# for finding fixed mutations between samples
filter_parameter_sample_across_sites = {\
                                        'min_average_coverage_to_include_sample': 2,  # low, since aDNA
}

filter_parameter_site_per_sample = {\
                                    'min_maf_for_call' : 0.9, #on individual samples, calls
                                    'min_cov_per_strand_for_call' : 1, 
                                    'min_qual_for_call' : 30,  #on individual samples, calls
                                    'min_cov_on_pos' : 2,
                                    'min_indel_gl_diff': 20,
                                    'min_indel_af': 0.8,
                                    'max_prop_0_covg_ancient': 0.9,
                                    'max_percentile_cov_ancient': 0.95}

filter_parameter_site_across_samples = {\
                                        'max_fraction_ambiguous_samples' : 0.25, #across samples per position
                                        'min_median_coverage_position' : 2, #across samples per position
                                        'correlation_threshold_recombination' : 0.98,
                                        'distance_for_recombination_checking': 500}

indel_filtering_params={'min_covg_on_pos': 3,
                        'major_allele_freq_indel_call': 0.85,
                        'min_gl_diff': 5,
                        'max_fraction_ambiguous_samples':0.25}

optional_filtering = {
    'heterozygosity': 'Ancient',
    'recombination': 'All',
    'aida_filtering': ['Repeat','RNA'],
    'minoraf_covariance': 95
}

heterozygosity_status=optional_filtering['heterozygosity']
recombination_status=optional_filtering['recombination']
aida_excluded=''.join(optional_filtering['aida_filtering'])
minAF_per=optional_filtering['minoraf_covariance']
strand_status=str(filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
maf_status=str(filter_parameter_site_per_sample['min_maf_for_call']).split('.')[-1]

analysis_params_output_name=f'strand_{strand_status}_maf_{maf_status}_heterozygosity{heterozygosity_status}_recombination{recombination_status}_aidaexclusions{aida_excluded}_minorafcovper{minAF_per}'

filtering_dicts=[filter_parameter_sample_across_sites,filter_parameter_site_per_sample,filter_parameter_site_across_samples,indel_filtering_params,optional_filtering]

with open(f'{analysis_params_output_name}_full_parameters.txt', 'w') as f:
    f.write(f'Full parameters for files named {analysis_params_output_name}')
    f.write('\n')
    for x in filtering_dicts:
        for key,value in x.items():
            f.write(f'{key}: {value}')
            f.write('\n')

#FQ_cutoff=60; #min -FQ that samples supporting both alleles must have
#min_average_coverage_to_include_sample = 8;
#min_maf_for_call = .85; #on individual samples, calls
#min_cov_per_strand_for_call = 2;  #on individual samples, calls # off since paired end and single end data
#min_qual_for_call = 30;  #on individual samples, calls
#max_fraction_ambiguous_samples = .25; #across samples per position
#min_median_coverage_position = 3; #across samples per position

## how far upstream of the nearest gene to annotate something a promoter
## mutation (not used if no annotation)
promotersize=250;

# %% 
# Define variables to store subject-specific results 
# =============================================================================
NTs = np.array(['A','T','C','G'],dtype=object) # NTs='ATCG'

para_evo_cand = np.array([],dtype=object)
mut_genes_df = {}
snp_freq_shifts_anno_lod = []
regression_res_lod = [] # stores all the results of the regression, and turned to pandas DF and saved
annotation_mutation_allParaSignal = {} # store all annotation_mutation for all detected candidates
mol_clock_data_dc = {} # dc to store for each subject the calculated data for mol clock estimation
allP_gene_paraevo = {} # store all para evo candidates that fullfill bonferroni

######################################################
### SETUP DONE ###### SETUP DONE ###### SETUP DONE ###
######################################################

######################################################
## START PROCESSING DATA #### START PROCESSING DATA ##
######################################################
# %% 
# Read in genome information
# =============================================================================
[chrStarts, genomeLength, scafNames] = apyp.genomestats(ref_genome_folder);
refgenome='Yersinia_pestis_CO92'
# # %% Load data from candidate_mutation_
# =============================================================================
#py_file='candidate_mutation_table.pickle.gz'
#[quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats] = apy.read_candidate_mutation_table_pickle_gzip(py_file)
    # # %% Breakpoint: candidate_mutation_table problems
    # # currenty only called in function when counts lacks three dimensions
    # if failFlag:
    #     continue

py_file_indels='final_cmt/candidate_mutation_table.pickle.gz'
[quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats,indel_p,indel_depth,indel_support,indel_identities,indel_index_for_identites] = apyp.read_candidate_mutation_table_pickle_gzip(py_file_indels)

p_chr_indices=[0]+[np.max(np.where(p < x + 1))+1 for x in chrStarts[1:]]

# %% 
# Mean cov per samples based on all positions in counts (aka p)
# =============================================================================

a = np.column_stack( (sampleNames,np.mean( np.sum(counts,axis=1),axis=1)) )
mean_cov_p_per_sample = np.array(['sample','mean_p_cov'])

mean_cov_p_per_sample = np.vstack( (mean_cov_p_per_sample,a) )    
np.savetxt('./mean_cov_per_sample.csv', mean_cov_p_per_sample, delimiter=',',fmt='%s')
#    # # analysis/plot /Users/u_key/Documents/mit/stapAD/readme/staphAD_figures_tables.r l.1040+

# %% 
# convert outputs to more permanent datatypes
# =============================================================================

sampleNames_all=np.asarray(sampleNames,dtype=object)
quals_all=-quals;
counts_all=counts;
coverage_all = counts_all.sum(axis=1).transpose() # axis=1 == rows; transpose needed > rows: pos and col: samples
indel_depth_all=indel_depth
indel_support_all=indel_support

# convert indels from pileups
# =============================================================================
# % indel counter reduction to indel count
# indel_counter: The first statistic is the number of reads (at this position and in this sample) that 
# support an indel. The second statistics is the number of reads (at this position and in this sample) that support a deletion.
# for now we need only the count for indels overall 
indel_counter = indel_counter[:,0,:].transpose()
# indel counter >> 50% indel row:1; row2 >> deletion

# %% Display samples the NOT fullfill min_average_coverage_to_include_sample:
# =============================================================================
 
lowcovsamples = sampleNames_all[ np.mean(coverage_all, axis=0) < filter_parameter_sample_across_sites['min_average_coverage_to_include_sample'] ]
print("The following samples will be excluded:", lowcovsamples)

# %% Parse regions to mask from SNP and indel calling
# =============================================================================
if "None" not in optional_filtering['aida_filtering']:
    possibilities_for_exclusion={
        'Morelli':'morelli_exclude_region',
        'Repeat':'repeat_region',
        'Noncore':'non_core_region',
        'Mcmaster':'outside_intersection_with_mcmaster',
        'Referr':'ref_error',
        'RNA':'RNA',
        'False':'SNP'
    }
    aida_filtering=set()
    parsed_exclusion = [possibilities_for_exclusion[x] for x in optional_filtering['aida_filtering']]
    with open('/Users/ad_loris/Documents/key_lab/test/regions2exclude_with-non-core_morelli-exclude_ref-errors_outsideMcMasterIntersect_lnbapaperfalseSNPS_20211125.gff') as f:
        for line in f.readlines():
            entries=line.strip().split('\t')
            if len(entries)>1:
                for exclusion in parsed_exclusion:
                    if exclusion in entries[2]:
                        for i in range(int(entries[3]),int(entries[4])+1):
                            aida_filtering.add(i)

    aida_filtering_not_done=set()
    parsed_exclusion_not_done = [possibilities_for_exclusion[x] for x in possibilities_for_exclusion if x not in parsed_exclusion]
    with open('/Users/ad_loris/Documents/key_lab/test/regions2exclude_with-non-core_morelli-exclude_ref-errors_outsideMcMasterIntersect_lnbapaperfalseSNPS_20211125.gff') as f:
            for line in f.readlines():
                entries=line.strip().split('\t')
                if len(entries)>1:
                    for exclusion in parsed_exclusion_not_done:
                        if exclusion in entries[2]:
                            for i in range(int(entries[3]),int(entries[4])+1):
                                aida_filtering_not_done.add([i])
# %% 
# Define goodsamples and filter data, good samples have > avg coverage
# =============================================================================
goodsamples =  np.all([coverage_all.mean(axis=0) >= filter_parameter_sample_across_sites['min_average_coverage_to_include_sample']],axis=0)

#Breakpoint: Too few samples passed filter, checking that at least 2 samples pass QC
if np.sum(goodsamples) < 2:
    print("Too few samples fullfill filter criteria! >> skip: " + refgenome)

sampleNames = sampleNames_all[goodsamples]
counts = counts_all[goodsamples , : , : ] # keep only level (samples) that fullfil filter!
quals = quals_all[ : , goodsamples ]
coverage = coverage_all[ : ,goodsamples]
indels = indel_counter[:,goodsamples]

indel_depth=indel_depth_all[:,goodsamples,:]
indel_support=indel_support_all[:,goodsamples]
indel_index_for_identites=indel_index_for_identites[:,goodsamples]
indel_total_depth=np.nansum(indel_depth,axis=2)

num_samples = len(sampleNames)

coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()

ancient_sample_bool = np.in1d(sampleNames,ancient_sample_names)

# %% 
# apply filtering parameters on indels
indel_total_depth=np.nansum(indel_depth,axis=2)

## getting goodpos for indels based on freebayes calls
## indel_depth == (pos x samples x (reads supporting highest GL indel x reads supporting all other calls))
has_indel_call_support=(indel_support > indel_filtering_params['min_gl_diff']) \
    & (indel_total_depth >= indel_filtering_params['min_covg_on_pos']) \
    & (indel_depth[:,:,0]/indel_total_depth > indel_filtering_params['major_allele_freq_indel_call']) \
    
has_indel_call_support[(np.sum(np.isnan(indel_index_for_identites),axis=1)/indel_index_for_identites.shape[1] > indel_filtering_params['max_fraction_ambiguous_samples']),:] = False

indel_positions_aida_filtering=[]
for idx,pos_in_p in enumerate(indel_p):
    if pos_in_p[0]+1 in aida_filtering:
        indel_positions_aida_filtering.append(idx)
indel_positions_aida_filtering=np.array(indel_positions_aida_filtering)
if len(indel_positions_aida_filtering) > 0:
    has_indel_call_support[indel_positions_aida_filtering,:] = False

has_indel = has_indel_call_support & ~np.isnan(indel_index_for_identites) & (indel_index_for_identites !=0)

indel_sizes_called = np.empty(has_indel.shape)
for index,row in enumerate(indel_index_for_identites):
    ref_len = len(indel_identities[index][0])
    for sample_index,value in enumerate(row):
        if not np.isnan(value):
            indel_sizes_called[index,sample_index]=len(indel_identities[index][int(value)]) - ref_len
        else: 
            indel_sizes_called[index,sample_index]=0


## define goodpos indels by finding all positions which have >1 indel type segregating (indel(s) or indel(s) and 0)
goodpos_indels=np.where(np.sum(has_indel,axis=1)>0)[0]

candidate_indels=[]
for i in goodpos_indels:
    if len(np.unique(indel_sizes_called[i][~np.isnan(indel_sizes_called[i])])) > 1:
        candidate_indels.append(i)

   
# %% 
# Extract refnt and define out/in-group bools
# =============================================================================

## Note ancnti/outs_nti defined below after filtered calls has been generated!
## get reference allele for all p; NOTE: analysis.m stored ATCG as double-digit-numeric
# use ref nt for mutation calls. important if multiple outgroups called 
refnt = apyp.extract_outgroup_mutation_positions(ref_genome_folder, apyp.p2chrpos(p,chrStarts));
refnti = apyp.nts2idx(refnt)
refnti_m = np.tile(refnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele
# print(np.unique(refnt)) # sanity check

# When no outgroup defined: refnt ~= ancnt:
#ancnt = refnt   

## Estimate outgroup (ancestral allele) from ALL samples added as outgroup to SM pipeline (ancnti* == major allele)
##  NOTE: NOTE: NOTE:

## TODO: update how this is done
##  NOTE: 
# CAN ALSO JUST GRAB OUTGROUP BY OUTPUTTING THE SAMPLES CSV from case
outgroup_bool = np.isin(sampleNames,['YAC'])
outgroup_name = sampleNames[outgroup_bool]
outgroup_idx=np.nonzero(outgroup_bool)[0]

# ingroup array (bool, idx) used later
ingroup_bool = np.invert(outgroup_bool)
ingroup_idx = np.nonzero(ingroup_bool)[0]

## get subsets of LNBA branch (clades) to find mutations
#likely inputs
# %% convert outgroup indels
indel_support[:,outgroup_idx]=np.nan
indel_sizes_called[:,outgroup_idx]=np.nan

# %% Find positions with fixed mutations
# 
# Extract allele and frequencies
contig_positions = apyp.p2chrpos(p,chrStarts) # 1col: chr, 2col: pos on chr; for all p
[maf, maNT, minorNT, minorAF] = apyp.div_major_allele_freq(counts) 
calls = maNT
intermediate_hasmutation=(calls != refnti_m) & (calls < 4)
intermediate_hasmutation[:,outgroup_bool] = False
print("Initial number of candidate variants: ", len(np.where( np.sum(intermediate_hasmutation, axis=1) > 0 )[0]))

# NOTE: function assumes first 8 rows in counts == 4nucl fwd&rev! watch out if extended counts used!
# NOTE: maf==0 -> no data;minorAF==0 -> no minor allele/ or no major allele; NT number corresponds to index in NTs [ATCG] or if maf==0 > NA == 4  
    
#  Make some basic structures for finding mutations
mutantAF = np.zeros(maNT.shape)
mutantAF[maNT != refnti_m] = maf[ maNT != refnti_m]; 
## generate mutantAF --> all positions that are not the reference, fetch major allele frequency, if reference at this pos, put zero

# mutantAF[ (minorNT != refnti_m) & (minorNT != 4) ] = mutantAF[  (minorNT != refnti_m) & (minorNT != 4)] + minorAF[  (minorNT != refnti_m) & (minorNT != 4) ] #this construction allows for positions with two different non-ancestral-mutations (when analysing data of more than one colony...plate sweeps)   

# Define mutations we do not trust in each and across samples.
# goodpos are indices of p that we trust

#
#
#
# Witin sample checks
 ## Filter per mutation
failed_quals = (quals < filter_parameter_site_per_sample['min_qual_for_call'])
failed_maf=(maf < filter_parameter_site_per_sample['min_maf_for_call'])
failed_forward=(coverage_forward_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
failed_reverse=(coverage_reverse_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
failed_cov=(coverage_reverse_strand + coverage_forward_strand < filter_parameter_site_per_sample['min_cov_on_pos'])
failed_indels=(indels > (0.5*coverage) )

# update calls:
failed_within_sample= (failed_quals | failed_maf | failed_forward | failed_reverse | failed_cov | failed_indels)
calls[failed_within_sample]=4

intermediate_hasmutation_post_sample_checks=(calls != refnti_m) & (calls < 4)
intermediate_hasmutation_post_sample_checks[:,outgroup_bool] = False
intermediate_num_variants=len(np.where( np.sum(intermediate_hasmutation_post_sample_checks, axis=1) > 0 )[0])
variants_removed_failed_within_sample=len(p)-intermediate_num_variants
print("Number of variants after this filtering step: ", intermediate_num_variants)
print("Number of variants removed: ", variants_removed_failed_within_sample)

#
#
#
# Metagenomic checks
# aDNA checks (coverage percentile, and )
# generated by file /Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/scripts_master/raw_data_processing/generate_bedcoverage_files.sh
bed_histogram_path='bed_files/final_out_files/*_genome_coverage_hist.tsv.gz'
bed_zero_covg_path = 'bed_files/final_out_files/*_merged_zero_covg_regions.tsv.gz'

failed_genomic_islands = apyp.filter_bed_0_cov_regions(bed_zero_covg_path,p,scafNames,chrStarts,sampleNames,filter_parameter_site_per_sample['max_prop_0_covg_ancient'])
failed_genomic_islands[:,~ancient_sample_bool]=False

failed_coverage_percentile=apyp.filter_bed_cov_hist(bed_histogram_path,p,scafNames,chrStarts,sampleNames,coverage,filter_parameter_site_per_sample['max_percentile_cov_ancient'],two_tailed=False,upper=True)
failed_coverage_percentile[:,~ancient_sample_bool]=False

## Removing SNP calls that are very close to nearby heterozygosity per sample (10bp default)
if optional_filtering['heterozygosity'] == 'All':
    failed_heterozygosity=apyp.find_calls_near_heterozygous_sites(p,minorAF,10,0.1)
elif optional_filtering['heterozygosity'] == 'Ancient':
    failed_heterozygosity=apyp.find_calls_near_heterozygous_sites(p,minorAF,10,0.1)
    failed_heterozygosity[:,~ancient_sample_bool]=False
else: failed_heterozygosity = np.full(calls.shape, False)

failed_metagenomic=(failed_genomic_islands | failed_coverage_percentile | failed_heterozygosity)
calls[(failed_metagenomic | failed_within_sample)]=4

intermediate_hasmutation_post_sample_checks_and_metagenomics=(calls != refnti_m) & (calls < 4)
intermediate_hasmutation_post_sample_checks_and_metagenomics[:,outgroup_bool] = False
intermediate_num_variants=len(np.where( np.sum(intermediate_hasmutation_post_sample_checks_and_metagenomics, axis=1) > 0 )[0])
variants_removed_failed_within_sample_failed_metagenomics=len(p)-variants_removed_failed_within_sample-intermediate_num_variants

print("Number of variants after this filtering step: ", intermediate_num_variants)
print("Number of variants removed: ", variants_removed_failed_within_sample_failed_metagenomics)

# save positions and samples where a metagenomic variant was filtered
np.where(intermediate_hasmutation_post_sample_checks_and_metagenomics != intermediate_hasmutation_post_sample_checks)

with open(f'{analysis_params_output_name}_masked_metagenomic_snvs.tsv', 'w') as f:
    f.write(f'Masked metagenomics snvs on ancient samples')
    f.write('\n')
    for x,y in zip(apyp.p2chrpos(p[np.where(intermediate_hasmutation_post_sample_checks_and_metagenomics != intermediate_hasmutation_post_sample_checks)
[0]],chrStarts),sampleNames[np.where(intermediate_hasmutation_post_sample_checks_and_metagenomics != intermediate_hasmutation_post_sample_checks)
[1]]):
        true_chrom=scafNames[x[0]-1]
        position_on_chrom=x[1]+1
        f.write(f'{y}\t{true_chrom}:{position_on_chrom}\n')
   

#
#
#
# Across sample checks 
## Remove putative recombinants
failed_any_QC = ( failed_metagenomic | failed_within_sample)

if optional_filtering['recombination'] == 'All':
    failed_recombinants = apyp.findrecombinantSNPs(p,mutantAF[ : , ingroup_bool ],filter_parameter_site_across_samples['distance_for_recombination_checking'],.99, failed_any_QC[ : , ingroup_bool ] )[1]
elif optional_filtering['recombination'] == 'Ancient':
    failed_recombinants = apyp.findrecombinantSNPs(p,mutantAF[ : , ancient_sample_bool ],filter_parameter_site_across_samples['distance_for_recombination_checking'],.99, failed_any_QC[ : , ancient_sample_bool ] )[1]


## Filter per site across samples
# Ignore here outgroup samples!
siteFilt = np.any(( (calls[:,ingroup_bool]>3).sum(axis=1) >= ((num_samples-np.sum(outgroup_bool)) * filter_parameter_site_across_samples['max_fraction_ambiguous_samples']) \
                        ,np.median( coverage[:,ingroup_bool], axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'] ),axis=0)
calls[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     

failed_optional_siteFilt = np.full(calls.shape, False)
if optional_filtering['recombination'] != 'None':
    failed_optional_siteFilt[failed_recombinants ,:] = 4

if optional_filtering['minoraf_covariance']:
    cov_scores=[]
    for x in minorAF:
        cov_scores.append(np.cov(x))
    minorAF_covariance_failed=np.where(np.array(cov_scores) > np.percentile(np.array(cov_scores),optional_filtering['minoraf_covariance']))[0]
    failed_optional_siteFilt[minorAF_covariance_failed,:]=4
## Removing SNPs in identified repeat regions, or morelli regions (from aida Valutena 2022)
positions_to_mask_aida_filtering=[]
for idx,pos_in_p in enumerate(p):
    if pos_in_p+1 in aida_filtering:
        positions_to_mask_aida_filtering.append(idx)
positions_to_mask_aida_filtering=np.array(positions_to_mask_aida_filtering)

if len(positions_to_mask_aida_filtering)>0:
    failed_optional_siteFilt[ positions_to_mask_aida_filtering ,:] = True # sites that fail True, for all samples incl. outgroup     
calls[ failed_optional_siteFilt ] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     


# NOTE: func below takes forever with many SNPs...saved below
[mutQual, mutQualIsolates] = apyp.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
mutQual = np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning

# translate filtered calls of ingroup into goodpos. mutations we believe. fixedmutation part removed in v6.
hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!
hasmutation[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

# calculate last set of removed variants
variants_removed_across_sites=intermediate_num_variants-len(candpos)
print("Number of variants after this filtering step: ", len(candpos))
print("Number of variants removed: ", variants_removed_across_sites)

goodpos = candpos
print(goodpos.size,'goodpos found.')

# create copies of to-be-updated values in ancient singleton checks (for validation):
counts_copy=counts.copy()
calls_copy=calls.copy()

# %%
# # # # 
# # # # # # # # # # # # # # # # 
# ancient singleton validation # 
# # # # # # # # # # # # # # # # 
# # # #
ancient_singletons = np.intersect1d(candpos,np.where((np.sum(hasmutation[:,np.in1d(sampleNames,ancient_sample_names)],axis=1)==1) & (np.sum(hasmutation[:,~np.in1d(sampleNames,ancient_sample_names)],axis=1)==0))[0])

for pos in ancient_singletons:
    samples=[x for x in sampleNames[hasmutation[p==p[pos]][0]] if x in ancient_sample_names]
    print(p[pos], samples)
    

sample_pos_to_output={}
for pos in ancient_singletons:
    samples=[x for x in sampleNames[hasmutation[p==p[pos]][0]] if x in ancient_sample_names][0]
    if samples in sample_pos_to_output:
        sample_pos_to_output[samples].append(pos)
    else:
        sample_pos_to_output[samples] = [pos]

regions_to_inspect=[]
for s in sample_pos_to_output:
    for singletons_position in sample_pos_to_output[s]:
        chrom,position_on_chrom=apyp.p2chrpos([p[singletons_position]],chrStarts)[0]
        true_chrom=scafNames[chrom-1]
        true_position_on_chrom_1_based=position_on_chrom+1
        regions_to_inspect.append(f'{s},{true_chrom}:{true_position_on_chrom_1_based}-{true_position_on_chrom_1_based}')
np.savetxt('singletons_to_extract_mappinginfo.txt', np.array(regions_to_inspect,dtype=str), fmt='%s')

# # # # # # # # # # # # # # # # # # # # # # # # #
# %% BREAKPOINT FOR RUNNING BLAST+
# # # # # # # # # # # # # # # # # # # # # # # # #
# %% BREAKPOINT FOR RUNNING BLAST+

# # # # # # # # # # # # # # # # # # # # # # # # #
# %% BREAKPOINT FOR RUNNING BLAST+

# # # # # # # # # # # # # # # # # # # # # # # # #
# %% BREAKPOINT FOR RUNNING BLAST+

# run blast on identified ancient singleton reads
# local_analysis/analysis_py/blast_execution_ancient_singletons.sh

# # # # # # # # # # # # # # # # # # # # # # # # #
# %%

# LOAD BACK IN RESULTS
def reparse_singleton_mpileup_calls_to_reads(blast_result_file):
    # TODO: move to apy
    fields=['Chrom','Pos','Ref','Covg','Calls','BQ','MQ','Read_Pos','Reads']
    pileup_file=pd.read_csv(blast_result_file.replace('_out.tsv.gz','.pileup').replace('blast_results','mpileup_with_readnames'),sep='\t', names=fields)
    reads=pileup_file.Reads.iloc[0].split(',')
    reference_call=pileup_file.Ref.iloc[0]
    # reparse calls
    raw_calls=pileup_file.Calls.iloc[0].replace('.',reference_call).replace(',',reference_call.lower())
    cleaned_calls=[]
    indices_to_skip=set()
    for idx,call in enumerate(raw_calls):
        if call in ['+', '-']: # skip indexing for all indel info
            size_of_indel=raw_calls[idx+1]
            indices_to_skip.add([i for i in range(idx,idx+size_of_indel+2)]) # exclude indices for eg +,2,A,A
        if call in ['^','$']: # skip indexing for all start/end info
            indices_to_skip.add(idx)
        if idx not in indices_to_skip:
            cleaned_calls.append(call)
    if len(cleaned_calls) != len(reads):
        raise ValueError
    return np.array(cleaned_calls),np.array(reads)

def update_counts_for_singleton(counts,call_support_remaining,sample_index_this_query,p_index_this_query):
    # update counts, calls, maf, and cov with cleaned read support for ancient singletons
    nt_order=['A','T','C','G','a','t','c','g']
    for idx,nt in enumerate(nt_order):
        cleaned_counts_this_nt=np.sum(call_support_remaining==nt)
        counts[sample_index_this_query,idx,p_index_this_query]=cleaned_counts_this_nt
    return counts


def update_matrices_post_singleton_counts_update(counts_updated,calls_previous,maNT_previous,refnti_m,filter_parameter_site_per_sample,sample_pos_to_output):
    # Recalculate entire datastructures    
    coverage_forward_strand_updated = counts_updated[:,0:4,:].sum(axis=1).transpose()
    coverage_reverse_strand_updated = counts_updated[:,4:8,:].sum(axis=1).transpose()
    [maf_updated, maNT_updated, minorNT_updated, minorAF_updated] = apyp.div_major_allele_freq(counts_updated) 
    calls_updated = maNT_updated
    mutantAF_updated = np.zeros(maNT_updated.shape)
    mutantAF_updated[maNT_updated != refnti_m] = maf_updated[ maNT_updated != refnti_m]; 
    failed_maf_updated=(maf_updated < filter_parameter_site_per_sample['min_maf_for_call'])
    failed_forward_updated=(coverage_forward_strand_updated < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
    failed_reverse_updated=(coverage_reverse_strand_updated < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
    failed_cov_updated=(coverage_reverse_strand_updated + coverage_forward_strand_updated < filter_parameter_site_per_sample['min_cov_on_pos'])
    failed_any_QC_updated = ( failed_maf_updated | failed_forward_updated | failed_reverse_updated | failed_cov_updated )
    calls_updated[failed_any_QC_updated] = 4
    # get indices to subset to only ancient singleton positions investigated:
    sample_indices_to_update_from_dict=[]
    pos_indices_to_update_from_dict=[]
    for sample in sample_pos_to_output:
        sample_index=np.where(sampleNames==sample)[0][0]
        sample_indices_to_update_from_dict+=[sample_index]*len(sample_pos_to_output[sample])
        pos_indices_to_update_from_dict+=sample_pos_to_output[sample]
    sample_indices_to_update_from_dict=np.array(sample_indices_to_update_from_dict)
    pos_indices_to_update_from_dict=np.array(pos_indices_to_update_from_dict)

    # only update to calls, maNT for the positions that are investigated
    calls_previous [ pos_indices_to_update_from_dict,sample_indices_to_update_from_dict ] = calls_updated[pos_indices_to_update_from_dict,sample_indices_to_update_from_dict]
    
    return [calls_previous,coverage_forward_strand_updated,coverage_reverse_strand_updated,maf_updated, calls_previous, minorNT_updated, minorAF_updated, mutantAF_updated]



blast_results=glob.glob('ancient_singleton_validation/blast_results/*.tsv.gz')

pseudo_taxids=[1649845,633,109458,273123,349747,502800,502801,748676,748679,1218101,1286084,1286085,1286086,1286087,1286088,1286089,1324587,1324608,1355443,748672]
pestis_taxids=[1035377,1234662,1345701,1345702,1345703,1345704,1345705,1345706,1345707,1345708,1345710,1455696,187410,214092,229193,349746,360102,377628,386656,547048,632,637382,637386,649716,748678]

singletons_to_mask=[]
singletons_to_mask_possibly_fine=[]
singletons_to_mask_sample=[]
counts_should_update=0
reads_failed={}
for blast_result_file in blast_results:
    blast_results_dict={}
    sampleid_for_query,chromid_posid_for_query=blast_result_file.split('/')[2].split('_NC_')
    chromid='NC_'+chromid_posid_for_query.split(':')[0]
    query=int(chromid_posid_for_query.split(':')[1].split('-')[0])
    cleaned_calls,cleaned_reads=reparse_singleton_mpileup_calls_to_reads(blast_result_file)
    pileup_file=pd.read_csv(blast_result_file.replace('_out.tsv.gz','.pileup').replace('blast_results','mpileup_with_readnames'),sep='\t', names=['Chrom','Pos','Ref','Covg','Calls','BQ','MQ','Read_Pos','Reads'])
    pileup_coverage=pileup_file.Covg.iloc[0]
    reads_in_counts=0
    with gzip.open(blast_result_file, 'rt') as f:
        for l in f:
            query_id,hit_id,hit_taxid,bitscore,pident=l.strip().split('\t')
            if query_id in pileup_file.Reads.iloc[0].split(','):
                reads_in_counts+=1
                if query_id not in blast_results_dict:
                    blast_results_dict[query_id]={(bitscore,pident):set([int(hit_taxid)])}
                elif (bitscore,pident) in blast_results_dict[query_id]:
                    blast_results_dict[query_id][(bitscore,pident)].add(int(hit_taxid))
    # parse what the MAF is in the supported read call
    p_index_this_query=np.where(p==apyp.chrpos2p(chromid,int(query)-1,scafNames,chrStarts))[0][0]
    if p_index_this_query not in ancient_singletons:
        print('indexing error')
    sample_index_this_query=np.where(sampleNames==sampleid_for_query)[0][0]
    maf_this_call=maf[p_index_this_query,sample_index_this_query]
    # number of pestis/pseudotb specific reads MUST exceed this proportion of reads to be retained
    reads_passed=[]
    reads_failed[sampleid_for_query]=[]
    reads_assessed=0
    for read_assessed in blast_results_dict:
        reads_assessed+=1
        for q,taxid_hits in blast_results_dict[read_assessed].items():
            taxid_hits_array=np.array([x for x in taxid_hits])
            if len(np.intersect1d(np.array(pseudo_taxids+pestis_taxids),taxid_hits_array)) > 0:
                reads_passed.append(read_assessed)
            else:
                reads_failed[sampleid_for_query].append(read_assessed)
    num_reads_passed=len(reads_passed)
    counts_coverage=np.sum(counts[sample_index_this_query,:,p_index_this_query])
    if num_reads_passed != counts_coverage:
        counts_should_update+=1
    print(f'{sampleid_for_query} {chromid}:{query} -- Passed: {num_reads_passed}, Assessed: {reads_assessed}, Counts Coverage: {counts_coverage}')
    if 0 < num_reads_passed < pileup_coverage:
        singletons_to_mask_possibly_fine.append(p_index_this_query)
        #print(f'Passed: {num_reads_passed}, Assessed: {reads_assessed}, Total Coverage: {pileup_coverage}')
    call_support_remaining=cleaned_calls[np.in1d(cleaned_reads,np.array(reads_passed))]
    counts=update_counts_for_singleton(counts,call_support_remaining,sample_index_this_query,p_index_this_query)

# outputting reads that failed check of origin to pestis or pseudotb
with open('blast_failed_reads_with_sample.tsv','w') as f:
    for s in reads_failed:
        if len(reads_failed[s])>0:
            reads_failed_parsed="\t".join(reads_failed[s])
            f.write(f'{s}\t{reads_failed_parsed}\n')

with open('blast_failed_reads.tsv','w') as f:
    for s in reads_failed:
        if len(reads_failed[s])>0:
            for r in reads_failed[s]:
                f.write(f'{r}\n')

# finding sample + position indices that have been updated
samples_indices_to_update=np.where(hasmutation[np.unique(np.where(counts!=counts_copy)[2])])[1]
indices_to_update=np.unique(np.where(counts!=counts_copy)[2])
# NOTE: ONLY these positions should be possibly updated in following function call

# %%
# Now final update to relevant matrices for saving
[calls, 
coverage_forward_strand, 
coverage_reverse_strand, 
maf, 
maNT,
minorNT, 
minorAF, 
mutantAF ] = update_matrices_post_singleton_counts_update(counts,calls,maNT,refnti_m,filter_parameter_site_per_sample,sample_pos_to_output)

# confirm no positions other than those set above have changed:
print("Only changed positions are within the set of indices to update?:",len(np.where(calls[indices_to_update,samples_indices_to_update]!=calls_copy[indices_to_update,samples_indices_to_update])[0]) == len(np.where(calls!=calls_copy)[0]))

hasmutation_final = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!

hasmutation_final[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
candpos_final = np.where( np.sum(hasmutation_final, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

goodpos_final = candpos_final
print('Number of ancient singletons removed: ',len(goodpos)-len(goodpos_final))
# saving indices of discarded ancient singletons:
with open(f'{analysis_params_output_name}_masked_singletons.tsv', 'w') as f:
    f.write(f'Masked singletons on ancient samples')
    f.write('\n')
    for x,y in zip(apyp.p2chrpos(p[np.where(calls!=calls_copy)[0]],chrStarts),sampleNames[np.where(calls!=calls_copy)[1]]):
        true_chrom=scafNames[x[0]-1]
        position_on_chrom=x[1]+1
        f.write(f'{y}\t{true_chrom}:{position_on_chrom}\n')


print(goodpos_final.size,'goodpos found.')

print("Carry out interactive QC investigations (below) for these positions and adjust parameters as necssary.")

# %% 
# output fully reduced and filtered CMT
cmtFile_sub_gp = 'candidate_mutation_table_processed_gp.pickle.gz'
with gzip.open('final_cmt/' + cmtFile_sub_gp, 'wb') as pickle_file:
    pickle.dump({
        'counts': counts[:, :, goodpos_final],
        'quals': quals[goodpos_final, :],
        'coverage_forward_strand': coverage_forward_strand[goodpos_final, :],
        'coverage_reverse_strand': coverage_reverse_strand[goodpos_final, :],
        'refnti_m': refnti_m[goodpos_final, :],
        'p': p[goodpos_final],
        'refgenome': refgenome,  
        'sampleNames': sampleNames,  
        'outgroup_bool': outgroup_bool,  
        'contig_positions': contig_positions[goodpos_final, :],
        'mutantAF': mutantAF[goodpos_final, :],
        'maf': maf[goodpos_final, :],
        'maNT': maNT[goodpos_final, :],
        'minorNT': minorNT[goodpos_final, :],
        'minorAF': minorAF[goodpos_final, :],
        'calls': calls[goodpos_final, :],
        'hasmutation': hasmutation[goodpos_final, :],
        'goodpos_raw_cmt': goodpos_final,
        'goodpos_cleaned_cmt': np.array(range(len(goodpos_final))),
        'analysis_params_output_name': analysis_params_output_name
        }, pickle_file)

# quickly test that dumping worked correctly:
with gzip.open('final_cmt/' + cmtFile_sub_gp, 'rb') as pickle_file:
    pickle_dictionary = pickle.load(pickle_file)
    counts_loaded = pickle_dictionary['counts' ]
    quals_loaded = pickle_dictionary['quals' ]
    coverage_forward_strand_loaded = pickle_dictionary['coverage_forward_strand' ]
    coverage_reverse_strand_loaded = pickle_dictionary['coverage_reverse_strand' ]
    refnti_m_loaded = pickle_dictionary['refnti_m' ]
    p_loaded = pickle_dictionary['p' ]
    refgenome_loaded = pickle_dictionary['refgenome' ]
    sampleNames_loaded = pickle_dictionary['sampleNames' ]
    outgroup_bool_loaded = pickle_dictionary['outgroup_bool' ]
    contig_positions_loaded = pickle_dictionary['contig_positions' ]
    mutantAF_loaded = pickle_dictionary['mutantAF' ]
    maf_loaded = pickle_dictionary['maf' ]
    maNT_loaded = pickle_dictionary['maNT' ]
    minorNT_loaded = pickle_dictionary['minorNT' ]
    minorAF_loaded = pickle_dictionary['minorAF' ]
    calls_loaded = pickle_dictionary['calls' ]
    hasmutation_loaded = pickle_dictionary['hasmutation' ]
    goodpos_final_raw_cmt_loaded = pickle_dictionary['goodpos_raw_cmt' ]
    goodpos_final_cleaned_cmt_loaded = pickle_dictionary['goodpos_cleaned_cmt' ]
    analysis_params_output_name_loaded = pickle_dictionary['analysis_params_output_name']


print('Checking parsing for counts - - - - - - - - - - - - -  PARSING OK?:',np.all(counts[:, :, goodpos_final] == counts_loaded))
print('Checking parsing for quals - - - - - - - - - - - - - - PARSING OK?:',np.all(quals[goodpos_final, :] == quals_loaded))
print('Checking parsing for coverage_forward_strand - - - - - PARSING OK?:',np.all(coverage_forward_strand[goodpos_final, :] == coverage_forward_strand_loaded))
print('Checking parsing for coverage_reverse_strand - - - - - PARSING OK?:',np.all(coverage_reverse_strand[goodpos_final, :] == coverage_reverse_strand_loaded))
print('Checking parsing for refnti_m - - - - - - - - - - - -  PARSING OK?:',np.all(refnti_m[goodpos_final, :] == refnti_m_loaded))
print('Checking parsing for p - - - - - - - - - - - - - - - - PARSING OK?:',np.all(p[goodpos_final] == p_loaded))
print('Checking parsing for refgenome - - - - - - - - - - - - PARSING OK?:',np.all(refgenome  == refgenome_loaded))
print('Checking parsing for sampleNames - - - - - - - - - - - PARSING OK?:',np.all(sampleNames  == sampleNames_loaded))
print('Checking parsing for outgroup_bool - - - - - - - - - - PARSING OK?:',np.all(outgroup_bool  == outgroup_bool_loaded))
print('Checking parsing for contig_positions - - - - - - - -  PARSING OK?:',np.all(contig_positions[goodpos_final, :] == contig_positions_loaded))
print('Checking parsing for mutantAF - - - - - - - - - - - -  PARSING OK?:',np.all(mutantAF[goodpos_final, :] == mutantAF_loaded))
print('Checking parsing for maf - - - - - - - - - - - - - - - PARSING OK?:',np.all(maf[goodpos_final, :] == maf_loaded))
print('Checking parsing for maNT - - - - - - - - - - - - - -  PARSING OK?:',np.all(maNT[goodpos_final, :] == maNT_loaded))
print('Checking parsing for minorNT - - - - - - - - - - - - - PARSING OK?:',np.all(minorNT[goodpos_final, :] == minorNT_loaded))
print('Checking parsing for minorAF - - - - - - - - - - - - - PARSING OK?:',np.all(minorAF[goodpos_final, :] == minorAF_loaded))
print('Checking parsing for calls - - - - - - - - - - - - - - PARSING OK?:',np.all(calls[goodpos_final, :] == calls_loaded))
print('Checking parsing for hasmutation - - - - - - - - - - - PARSING OK?:',np.all(hasmutation[goodpos_final, :] == hasmutation_loaded))
print('Checking parsing for goodpos_final_raw_cmt - - - - - - PARSING OK?:',np.all(goodpos_final == goodpos_final_raw_cmt_loaded))
print('Checking parsing for goodpos_final_cleaned_cmt - - - - PARSING OK?:',np.all(np.array(range(len(goodpos_final))) == goodpos_final_cleaned_cmt_loaded))
print('Checking parsing for analysis_params_output_name - - - PARSING OK?:', analysis_params_output_name == analysis_params_output_name_loaded)

# %% output fasta file for tree generation
# =============================================================================
goodpos_final2use = goodpos_final; # leave. downsize only goodpos_final2useTree below if necessary
order = np.arange(sampleNames.shape[0])
num_contigs = np.max(contig_positions[:,0]);
contig_lengths = np.append(chrStarts[1:], genomeLength) - chrStarts

# =============================================================================    
# define pos to use in tree
goodpos_final2useTree = goodpos_final #(1:1000); %trim for easier tree view; TDL called quality_positions

# get data and filter for goodpos_final
calls_for_treei = calls; 
calls_for_treei = calls_for_treei[ goodpos_final2useTree, : ]
calls_for_treei_no_outgroup = calls_for_treei[ : , ingroup_bool ]

# build sampleNames  (that passed filter, see above from sampleNames_all --> sampleNames) w/ metainfo
treesampleNamesLong_no_outgroup = sampleNames[ingroup_bool] # remove outgroup samples
treesampleNamesLong = sampleNames # include all samples 

# sampleNamesDnapars : max 10c. Use numeric with 10c (works for up to 10^10 samples! )
sampleNamesDnapars_no_outgroup = [ "{:010d}".format(i) for i in range(len(treesampleNamesLong_no_outgroup))]
sampleNamesDnapars = [ "{:010d}".format(i) for i in range(len(treesampleNamesLong))] # includes all samples

# translate index to nucleotide
calls_for_tree = apyp.idx2nts(calls_for_treei) # ATCGN translation
calls_for_tree_no_outgroup= apyp.idx2nts(calls_for_treei_no_outgroup)
# add reference nucleotide for all positions

## outgroup already in samples, so unnecessary (CO92)
refgenome_nts = apyp.extract_outgroup_mutation_positions(ref_genome_folder, apyp.p2chrpos(p[goodpos_final2useTree],chrStarts));
refgenome_nts_for_tree=refnt[goodpos_final2useTree]
calls_for_tree_ref_no_outgroup = np.concatenate((apyp.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree_no_outgroup),axis=1) # first column now refgenome_nts; refgenome_nts[:, None] to make ndims (2) same for both
calls_for_tree_ref_outgroup = np.concatenate((apyp.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree),axis=1)

sampleNamesDnapars_ref_no_outgroup = np.append(['Sref'],sampleNamesDnapars_no_outgroup) # add name for outgroup
sampleNamesDnapars_ref_outgroup=np.append(['Sref'],sampleNamesDnapars)

treesampleNamesLong_ref_no_outgroup = np.append(['Sref'],treesampleNamesLong_no_outgroup) # add name for outgroup
treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)

apyp.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,f'{analysis_params_output_name}_full_qc')

# RUN HPC-BASED TREE GENERATION BASED ON THIS SCRIPT: 
print('Tree production ready! Run tree production on an HPC cluster using this script: analysis_py/logan_final_tree_generation.sh')
print('After running, move on to the final analysis script: analysis_py/pestis_evolution_LNBA_analysis.py')
print('Paths for the tree, finalized CMT, and metadata tables should be updated as necessary.')
print('Be sure to correctly import the downsampled candidate mutation table output above and the final tree from raxml run from HPC!')
# %%
################################################################################################################
######## MAIN QC DONE, FOLLOWING IS INTERACTIVE QC INVESTIGATIONS ##############################################
################################################################################################################
######## MAIN QC DONE, FOLLOWING IS INTERACTIVE QC INVESTIGATIONS ##############################################
################################################################################################################
######## MAIN QC DONE, FOLLOWING IS INTERACTIVE QC INVESTIGATIONS ##############################################
################################################################################################################
######## MAIN QC DONE, FOLLOWING IS INTERACTIVE QC INVESTIGATIONS ##############################################
################################################################################################################

# %%

# %% custom functions (mostly for plotting QC)
def find_overlapping_filtering(counts, quals, min_qual_for_call, min_maf_for_call, min_cov_per_strand,indels,coverage,indices_to_subset_by=False):
    if isinstance(indices_to_subset_by, bool):
        indices_to_subset_by = np.full((counts.shape[0]),True)
    # generate all filterings
    [maf, maNT, minorNT, minorAF] = apyp.div_major_allele_freq(counts)
    coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
    coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()

    collection_of_snp_filtering={
        'snps_removed_by_quals' : ((quals[:,indices_to_subset_by] < min_qual_for_call) & (maNT[:,indices_to_subset_by] < 4)),
        'snps_removed_by_maf' : ((maf[:,indices_to_subset_by] < min_maf_for_call) & (maNT[:,indices_to_subset_by] < 4)),
        'snps_removed_by_fwdstrand' : ((coverage_forward_strand[:,indices_to_subset_by] < min_cov_per_strand) & (maNT[:,indices_to_subset_by] < 4)),
        'snps_removed_by_revstrand' : ((coverage_reverse_strand[:,indices_to_subset_by] < min_cov_per_strand) & (maNT[:,indices_to_subset_by] < 4)),
        'snps_removed_by_indel' : ((indels[:,indices_to_subset_by] > (0.5*coverage[:,indices_to_subset_by])) & (maNT[:,indices_to_subset_by] < 4))
    }
    
    #generate overlaps matrix
    output_matrix={col:[] for col in collection_of_snp_filtering}
    for r_index,row in enumerate(collection_of_snp_filtering):
        for c_index,column in enumerate(collection_of_snp_filtering):
            overlap=(collection_of_snp_filtering[row] == collection_of_snp_filtering[column]) & (collection_of_snp_filtering[row]==True)
            output_matrix[column].append(np.sum(overlap)/np.sum(collection_of_snp_filtering[column]))
    output_df=pd.DataFrame(output_matrix)
    keys_with_n_removed=[]
    for key in collection_of_snp_filtering:
        n=np.sum(collection_of_snp_filtering[key])
        keys_with_n_removed.append(f'{key} n={n}')
    output_df.insert(0,"Proportion of column filter also removed by row filter", keys_with_n_removed)
    return output_df


def write_fasta_helper(sequences,scaffold_names,sample):
    fasta_to_write=[]
    for index,s in enumerate(sequences):
        fasta_to_write.append(SeqRecord.SeqRecord(Seq.Seq(''.join(s)),id=f'{scaffold_names[index]}'))
    SeqIO.write(fasta_to_write,f'{sample}_mapped_mutaitons.fasta','fasta')

def define_subsection_of_fasta(base_index_start_end,chrStarts,genomeLength):
    if base_index_start_end[1] > genomeLength:
        print('WARNING: End index exceeds length of genome. End point, set at end of genome.')
        base_index_start_end[1] = genomeLength
    chr_index=apyp.p2chrpos(np.array(base_index_start_end), chrStarts)
    if chr_index[0][0] != chr_index[1][0]:
        print(f'WARNING: Subset of genome to output provided spans multiple chromosomes.')
        print(f'End index will be set as the end of chromosome {chr_index[0][0]}')
        base_index_start_end[1] = chrStarts[chr_index[0][0]]
    return apyp.p2chrpos(np.array(base_index_start_end), chrStarts)
        
def map_sample_mutations_onto_reference_helper_snp_mapping(ancestral_reconstruction_tree,ancestral_reconstruction_fasta,sampleName,calls_on_goodpos,sampleNames):
    """
        Helper function for reconstructing region of genome with given set of mutations
        To be used to infer impact of promoter mutations (eg disruption of promoter motifs)
    """
    print(sampleName)
    tree=apyp.parse_tree(ancestral_reconstruction_tree)
    sample_clade_in_tree=[x for x in tree.find_elements(sampleName)]
    if len(sample_clade_in_tree)==0:
        print(f'No sample found in tree with name {sampleName}')
        return np.array([])
    name_of_parent_node=tree.get_path(sample_clade_in_tree[0])[-2].name
    parental_calls_on_goodpos = apyp.parse_ancestra_fasta_treetime(ancestral_reconstruction_fasta, name_of_parent_node)
    calls_this_sample = calls_on_goodpos[:,(sampleNames == sampleName)].flatten()
    to_update=(calls_this_sample != parental_calls_on_goodpos) & (calls_this_sample != 4)
    parental_calls_on_goodpos[to_update] = calls_this_sample[to_update]
    return parental_calls_on_goodpos


def main_map_sample_nt_mutations_onto_reference(ancestral_reconstruction_tree,ancestral_reconstruction_fasta,sampleName,sampleNames,p_on_goodpos,calls_on_goodpos,ref_genome_folder):
    [chrStarts, genomeLength, scafNames] = apyp.genomestats(ref_genome_folder)
    goodpos_calls=map_sample_mutations_onto_reference_helper_snp_mapping(ancestral_reconstruction_tree,ancestral_reconstruction_fasta,sampleName,calls_on_goodpos,sampleNames)
    if len(goodpos_calls)==0:
        return []
    chrpos_matrix_to_update = apyp.p2chrpos(p_on_goodpos, chrStarts)
    ref_genome_as_list_of_seqs=list(list(SeqIO.parse(f'{ref_genome_folder}/genome.fasta', format='fasta')))
    ref_genome_as_list_of_str = [list(str(x.seq)) for x in ref_genome_as_list_of_seqs]
    for index,chrpos in enumerate(chrpos_matrix_to_update):
        ref_genome_as_list_of_str[chrpos[0]-1][chrpos[1]] = NTs[goodpos_calls[index]]
    return ref_genome_as_list_of_str



"""
Not currently run, or implemented fully
def main_map_sample_mutations_onto_reference(ancestral_reconstruction_tree,ancestral_reconstruction_fasta,sampleName,sampleNames,p_on_goodpos,calls_on_goodpos,indel_p,has_indel_call_support,indel_identities,indel_calls_this_sample,ref_genome_folder,base_index_start_end=[],write_fasta=False):
    [chrStarts, genomeLength, scafNames] = apy.genomestats(ref_genome_folder)
    goodpos_calls=map_sample_mutations_onto_reference_helper_snp_mapping(ancestral_reconstruction_tree,ancestral_reconstruction_fasta,sampleName,calls_on_goodpos,sampleNames)
    chrpos_matrix_to_update = apyp.p2chrpos(p_on_goodpos, chrStarts)
    ref_genome_as_list_of_seqs=list(list(SeqIO.parse(f'{ref_genome_folder}/genome.fasta', format='fasta')))
    ref_genome_as_list_of_str = [list(str(x.seq)) for x in ref_genome_as_list_of_seqs]    
    for index,chrpos in enumerate(chrpos_matrix_to_update):
        ref_genome_as_list_of_str[chrpos[0]-1][chrpos[1]] = goodpos_calls[index]
    ref_genome_as_list_of_str,difference_ref_alt=map_sample_mutations_onto_reference_helper_indel_mapping(indel_p,has_indel_call_support,indel_identities,indel_calls_this_sample,sampleName,sampleNames,ref_genome_as_list_of_str,chrStarts,base_index_start_end)
    print(difference_ref_alt)
    if write_fasta:
        write_fasta_helper(ref_genome_as_list_of_str,scafNames,sample)
    else:
        if len(base_index_start_end) == 2:
            chr_pos_subset=define_subsection_of_fasta(base_index_start_end,chrStarts,genomeLength)
            return ref_genome_as_list_of_str[chr_pos_subset[0][0]-1][chr_pos_subset[0][1]:chr_pos_subset[1][1]+difference_ref_alt[chr_pos_subset[0][0]]]
        else:
            return ref_genome_as_list_of_str
"""
"""
# run for comparisons sake:
[calls_test, 
coverage_forward_strand_test, 
coverage_reverse_strand_test, 
maf_test, 
maNT_test,
minorNT_test, 
minorAF_test, 
mutantAF_test ] = update_matrices_post_singleton_counts_update(counts,calls,maNT,refnti_m,filter_parameter_site_per_sample,sample_pos_to_output)

hasmutation_test = (calls_test != refnti_m) & (calls_test < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!

hasmutation_test[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
candpos_test = np.where( np.sum(hasmutation_test, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

print(candpos_test.size,'goodpos found.')

# investigating positions that had counts adjusted but call stayed same
p_with_counts_adjusted=np.array([x[0] for x in singletons_to_mask_possibly_fine]).flatten()
samples_with_counts_adjusted=np.where(hasmutation[p_with_counts_adjusted])[1]
print('Singletons with adjusted but not removed counts:',len(p_with_counts_adjusted))
calls[p_with_counts_adjusted,samples_with_counts_adjusted]
pre_updated_singleton_indices,ancient_singleton_samples = np.where( hasmutation[ancient_singletons] )
pre_updated_singletons=ancient_singletons[pre_updated_singleton_indices]
print('\tSingletons with adjusted but not removed counts, but not adjusted call:',np.sum(calls[p_with_counts_adjusted,samples_with_counts_adjusted] == calls_test[p_with_counts_adjusted,samples_with_counts_adjusted]))
print('\tSingletons with adjusted but not removed counts and adjusted call updated total:', np.sum(calls[p_with_counts_adjusted,samples_with_counts_adjusted] != calls_test[p_with_counts_adjusted,samples_with_counts_adjusted]))
# changed to ref?
refnt_singleton_changed=refnti[p_with_counts_adjusted[(calls[p_with_counts_adjusted,samples_with_counts_adjusted] != calls_test[p_with_counts_adjusted,samples_with_counts_adjusted])]]
num_refnt_changes=np.sum(refnt_singleton_changed==calls_test[p_with_counts_adjusted,samples_with_counts_adjusted][calls[p_with_counts_adjusted,samples_with_counts_adjusted] != calls_test[p_with_counts_adjusted,samples_with_counts_adjusted]])
num_non_refnt_changes=np.sum((refnt_singleton_changed!=calls_test[p_with_counts_adjusted,samples_with_counts_adjusted][calls[p_with_counts_adjusted,samples_with_counts_adjusted] != calls_test[p_with_counts_adjusted,samples_with_counts_adjusted]])&(calls_test[p_with_counts_adjusted,samples_with_counts_adjusted][calls[p_with_counts_adjusted,samples_with_counts_adjusted] != calls_test[p_with_counts_adjusted,samples_with_counts_adjusted]] != 4))
num_changes_to_n=np.sum(calls_test[p_with_counts_adjusted,samples_with_counts_adjusted][calls[p_with_counts_adjusted,samples_with_counts_adjusted] != calls_test[p_with_counts_adjusted,samples_with_counts_adjusted]] == 4)
print('\t\tSingleton call changed to ref:', num_refnt_changes)
print('\t\tSingleton call change to other non-ref allele:', num_non_refnt_changes)
print('\t\tSingleton call change to N:', num_changes_to_n)

# call changed 
updated_singleton_indices,updated_samples = np.where( calls[ancient_singletons] != calls_test[ancient_singletons])
updated_singletons=ancient_singletons[updated_singleton_indices]
# investigating positions that were adjusted but call changed

# investigating positions that were removed
num_total_removed_singletons=np.sum(calls_test[updated_singletons,updated_samples]==4)
print('Overall singletons masked (call updated to 4):', num_total_removed_singletons)
print('\tMasked calls due to no reads passing: ', num_total_removed_singletons-num_changes_to_n)

"""


# %% 
# Explore impact of various filtering steps
# set parameters
filter_parameter_site_per_sample_QC = {'min_cov_on_pos': 4, 'min_qual_for_call': 30, 'min_maf_for_call': 0.9, 'max_prop_0_covg_ancient':0.9,'max_percentile_cov_ancient':0.95}

plot_upset=True

## ONLY EDIT ABOVE FILTERING PARAMETERS:
# set of mutations should no filtering be done:
[maf_copy, maNT_copy, minorNT_copy, minorAF_copy] = apy.div_major_allele_freq(counts) 
calls_copy=maNT_copy
hasmutation_no_filtering = (calls_copy != refnti_m) & (calls_copy < 4)
hasmutation_no_filtering[:,outgroup_bool] = False
already_non_mutated_no_filtering=np.where(np.sum(hasmutation_no_filtering,axis=1)==0)[0]

# failed coverage
failed_cov=(coverage_reverse_strand + coverage_forward_strand < filter_parameter_site_per_sample_QC['min_cov_on_pos'])
removed_cov=np.where((np.sum(hasmutation_no_filtering,axis=1)-np.sum((failed_cov & hasmutation_no_filtering),axis=1))==0)[0]


# failed quals
failed_quals= (quals < filter_parameter_site_per_sample_QC['min_qual_for_call'])
removed_quals=np.where((np.sum(hasmutation_no_filtering,axis=1)-np.sum((failed_quals & hasmutation_no_filtering),axis=1))==0)[0]

# failed MAF
failed_MAF= (maf_copy < filter_parameter_site_per_sample_QC['min_maf_for_call'])
removed_MAF=np.where((np.sum(hasmutation_no_filtering,axis=1)-np.sum((failed_MAF & hasmutation_no_filtering),axis=1))==0)[0]

# Failed genomic islands
failed_genomic_islands=apyp.filter_bed_0_cov_regions(bed_zero_covg_path,p,scafNames,chrStarts,sampleNames,filter_parameter_site_per_sample_QC['max_prop_0_covg_ancient'])
removed_genomic_islands=np.where((np.sum(hasmutation_no_filtering,axis=1)-np.sum((failed_genomic_islands & hasmutation_no_filtering),axis=1))==0)[0]
num_removed_genomic_island=len(removed_genomic_islands)-len(already_non_mutated_no_filtering)
# Failed coverage percentile
failed_coverage_percentile=apyp.filter_bed_cov_hist(bed_histogram_path,p,scafNames,chrStarts,sampleNames,coverage,filter_parameter_site_per_sample_QC['max_percentile_cov_ancient'],two_tailed=False,upper=True)
removed_coverage_percentile=np.where((np.sum(hasmutation_no_filtering,axis=1)-np.sum((failed_coverage_percentile & hasmutation_no_filtering),axis=1))==0)[0]
num_removed_coverage_percentile=len(removed_coverage_percentile)-len(already_non_mutated_no_filtering)

# Failed indel
failed_indel = (indels > (0.5*coverage) )

## add sitefile, all filters, and those removed after all filter by mutqual
# get upset plot data
data_out,filtercomb_out=apyp.combinatorial_filtering_to_goodpos(calls_copy,refnti_m,mutQual_copy,[failed_cov,failed_quals,failed_MAF,failed_genomic_islands,failed_coverage_percentile],[x[0]+'='+str(x[1]) for x in filter_parameter_site_per_sample_QC.items()])

# calculate site filters applied across sites, after applying above filters
calls_copy[(failed_coverage_percentile | failed_genomic_islands | failed_MAF | failed_quals | failed_cov | failed_indel) ] = 4
siteFilt = np.any(( (calls_copy[:,ingroup_bool]>3).sum(axis=1) >= ((num_samples-np.sum(outgroup_bool)) * filter_parameter_site_across_samples['max_fraction_ambiguous_samples']) \
                        ,np.median( coverage[:,ingroup_bool], axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'] ),axis=0)

calls_copy[ siteFilt ,:] = 4

[mutQual_copy, mutQualIsolates_copy] = apy.ana_mutation_quality(calls_copy[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
mutQual_copy = np.nan_to_num(mutQual_copy, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning

hasmut_all_filters_pre_mutqual = (calls_copy != refnti_m) & (calls_copy < 4)
hasmut_all_filters = (calls_copy != refnti_m) & (calls_copy < 4) & (np.tile(mutQual_copy,(1,num_samples)) >= 1)
hasmut_all_filters_pre_mutqual[:,outgroup_bool] = False
hasmut_all_filters[:,outgroup_bool] = False


data_out+=[np.sum(np.sum(failed_indel,axis=1)>0),np.sum(siteFilt),np.sum(np.sum(hasmut_all_filters,axis=1)==0),np.sum(np.sum(hasmut_all_filters,axis=1)==0)-np.sum(np.sum(hasmut_all_filters_pre_mutqual,axis=1)==0)]
filtercomb_out+=[['failed_indels'],['site_filter_across_samples(after all but mutqual)'],filtercomb_out[-1]+['site_filter_across_samples(after all but mutqual)','mutqual','failed_indels'],['mutqual']]

if plot_upset:
    upset_data=upsetplot.from_memberships(filtercomb_out,data=data_out)
    upsetplot_out=upsetplot.UpSet(upset_data,sort_categories_by=None)
    upsetplot_out.plot()
    plt.title('Impact of filters on number of candidate loci removed from goodpos')
    plt.show()

all_filters_individually=np.setxor1d(np.unique(np.concatenate((removed_cov,removed_quals,removed_MAF,removed_genomic_islands,removed_coverage_percentile))),already_non_mutated_no_filtering)

goodpos_removed_all_filters=np.where(np.sum(hasmut_all_filters,axis=1)==0)[0]
num_goodpos_removed_all_filters=np.sum(np.sum(hasmut_all_filters,axis=1)==0)















# %%
# interactive QC plotting, investigation
%matplotlib qt


# investigating different covariance of heterozygous sites filtering
cov_scores_only_not_ancient=[]
for x in minorAF:
    only_ancient_x=x[~ancient_sample_bool]
    cov_scores_only_not_ancient.append(np.cov(only_ancient_x))
minorAF_covariance_failed_only_not_ancient=np.where(np.array(cov_scores_only_not_ancient) > np.percentile(np.array(cov_scores_only_not_ancient),optional_filtering['minoraf_covariance']))[0]
cov_scores_only_ancient=[]
for x in minorAF:
    only_ancient_x=x[ancient_sample_bool]
    cov_scores_only_ancient.append(np.cov(only_ancient_x))
minorAF_covariance_failed_only_ancient=np.where(np.array(cov_scores_only_ancient) > np.percentile(np.array(cov_scores_only_ancient),optional_filtering['minoraf_covariance']))[0]


cov_minor_removed_above_only_non_ancient_cutoff=np.where(np.array(cov_scores_only_ancient) > np.percentile(np.array(cov_scores_only_not_ancient),99))[0]

subset_to_output=np.intersect1d(goodpos,np.setdiff1d(cov_minor_removed_above_only_non_ancient_cutoff,minorAF_covariance_failed))

ts = time.time()
timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
apy.plot_interactive_scatter_barplots(p[subset_to_output],mutQual[subset_to_output],'pos','qual', sampleNames,counts[:,:,subset_to_output],"28.07.2022", filter_parameter_site_across_samples['correlation_threshold_recombination'], "", refgenome, saveplots = False)
matplotlib.use('agg')

# %%














aida_filtering_not_done=set()
parsed_exclusion_not_done = [possibilities_for_exclusion[x] for x in possibilities_for_exclusion if x not in parsed_exclusion]
with open('/Users/ad_loris/Documents/key_lab/test/regions2exclude_with-non-core_morelli-exclude_ref-errors_outsideMcMasterIntersect_lnbapaperfalseSNPS_20211125.gff') as f:
        for line in f.readlines():
            entries=line.strip().split('\t')
            if len(entries)>1:
                for exclusion in parsed_exclusion_not_done:
                    if exclusion in entries[2]:
                        for i in range(int(entries[3]),int(entries[4])+1):
                            aida_filtering_not_done.add(i)

positions_to_mask_aida_filtering_not_done=[]
for idx,pos_in_p in enumerate(p):
    if pos_in_p+1 in aida_filtering_not_done:
        positions_to_mask_aida_filtering_not_done.append(idx)
positions_to_mask_aida_filtering_not_done=np.array(positions_to_mask_aida_filtering_not_done)
annotation_mutations_lnba_singletons_relative_MRCA[annotation_mutations_lnba_singletons_relative_MRCA.p.isin(p[positions_to_mask_aida_filtering_not_done])]


# %%

def find_differences_due_to_filtering(counts, filterset1, filterset2):
    results={}
    for index,filterset in enumerate([filterset1,filterset2]):
        [maf, maNT, minorNT, minorAF] = apy.div_major_allele_freq(counts) 
        my_calls = maNT
        my_calls= maNT.copy()
        my_calls[ filterset ] = 4
        siteFilt = np.any(( (my_calls[:,ingroup_bool]>3).sum(axis=1) >= ((num_samples-np.sum(outgroup_bool)) * filter_parameter_site_across_samples['max_fraction_ambiguous_samples']) \
                                ,np.median( coverage[:,ingroup_bool], axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'] ),axis=0)

        my_calls[ siteFilt ,:] = 4
        failed_optional_siteFilt = np.full(my_calls.shape, False)
        if optional_filtering['recombination'] != 'None':
            failed_optional_siteFilt[failed_recombinants ,:] = 4

        if optional_filtering['minoraf_covariance']:
            cov_scores=[]
            for x in minorAF:
                cov_scores.append(np.cov(x))
            minorAF_covariance_failed=np.where(np.array(cov_scores) > np.percentile(np.array(cov_scores),optional_filtering['minoraf_covariance']))[0]
            failed_optional_siteFilt[minorAF_covariance_failed,:]=4
        ## Removing SNPs in identified repeat regions, or morelli regions (from aida Valutena 2022)
        positions_to_mask_aida_filtering=[]
        for idx,pos_in_p in enumerate(p):
            if pos_in_p+1 in aida_filtering:
                positions_to_mask_aida_filtering.append(idx)
        positions_to_mask_aida_filtering=np.array(positions_to_mask_aida_filtering)

        if len(positions_to_mask_aida_filtering)>0:
            failed_optional_siteFilt[ positions_to_mask_aida_filtering ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     
        my_calls[ failed_optional_siteFilt ] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     


        [mutQual, mutQualIsolates] = apy.ana_mutation_quality(my_calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
        mutQual= np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning

        hasmutation = (my_calls != refnti_m) & (my_calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1)
        hasmutation[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
        results[index] = hasmutation
    return results


[maf, maNT, minorNT, minorAF] = apy.div_major_allele_freq(counts) 

# NOTE: function assumes first 8 rows in counts == 4nucl fwd&rev! watch out if extended counts used!
# NOTE: maf==0 -> no data;minorAF==0 -> no minor allele/ or no major allele; NT number corresponds to index in NTs [ATCG] or if maf==0 > NA == 4  
    
#  Make some basic structures for finding mutations
mutantAF = np.zeros(maNT.shape)
mutantAF[maNT != refnti_m] = maf[ maNT != refnti_m]; 
## generate mutantAF --> all positions that are not the reference, fetch major allele frequency, if reference at this pos, put zero

# mutantAF[ (minorNT != refnti_m) & (minorNT != 4) ] = mutantAF[  (minorNT != refnti_m) & (minorNT != 4)] + minorAF[  (minorNT != refnti_m) & (minorNT != 4) ] #this construction allows for positions with two different non-ancestral-mutations (when analysing data of more than one colony...plate sweeps)   

# Define mutations we do not trust in each and across samples.
# goodpos are indices of p that we trust
## Filter per mutation

failed_quals = (quals < filter_parameter_site_per_sample['min_qual_for_call'])
failed_maf=(maf < filter_parameter_site_per_sample['min_maf_for_call'])
failed_forward=(coverage_forward_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
failed_reverse=(coverage_reverse_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
failed_cov=(coverage_reverse_strand + coverage_forward_strand < filter_parameter_site_per_sample['min_cov_on_pos'])
failed_indels=(indels > (0.5*coverage) )
failed_heterozygosity=apyp.find_calls_near_heterozygous_sites(p,minorAF,10,0.1)

# aDNA checks (coverage percentile, and )
# generated by file /Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/scripts_master/raw_data_processing/generate_bedcoverage_files.sh
bed_histogram_path='bed_files/out_files/*_genome_coverage_hist.tsv.gz'
bed_zero_covg_path = 'bed_files/out_files/*_merged_zero_covg_regions.tsv.gz'

failed_genomic_islands = apyp.filter_bed_0_cov_regions(bed_zero_covg_path,p,scafNames,chrStarts,sampleNames,filter_parameter_site_per_sample['max_prop_0_covg_ancient'])
failed_coverage_percentile=apyp.filter_bed_cov_hist(bed_histogram_path,p,scafNames,chrStarts,sampleNames,coverage,filter_parameter_site_per_sample['max_percentile_cov_ancient'],two_tailed=False,upper=True)
failed_coverage_percentile[:,~ancient_sample_bool]=False



failed_all_but_percentile = (failed_quals | failed_maf | failed_genomic_islands | failed_forward | failed_reverse | failed_cov | failed_indels | failed_heterozygosity )
failed_all_but_island = (failed_quals | failed_maf | failed_coverage_percentile | failed_forward | failed_reverse | failed_cov | failed_indels | failed_heterozygosity )
failed_all = (failed_quals | failed_maf | failed_coverage_percentile | failed_genomic_islands | failed_forward | failed_reverse | failed_cov | failed_indels | failed_heterozygosity )

removed_differences_by_covg_percentile=find_differences_due_to_filtering(counts,failed_all,failed_all_but_percentile)
print(len(np.where( np.sum(removed_differences_by_covg_percentile[0], axis=1) > 0 )[0]))
print(len(np.where( np.sum(removed_differences_by_covg_percentile[1], axis=1) > 0 )[0]))
# find non-overlap:
removed_exclusively_by_covg_percentile=np.setdiff1d(np.where( np.sum(removed_differences_by_covg_percentile[1], axis=1) > 0 )[0],np.where( np.sum(removed_differences_by_covg_percentile[0], axis=1) > 0 )[0])

# print those missed by 
for position,line in zip(p[np.setdiff1d(removed_exclusively_by_covg_percentile,singletons_to_mask)],removed_differences_by_covg_percentile[1][np.setdiff1d(removed_exclusively_by_covg_percentile,singletons_to_mask)]):
    print(position,sampleNames[line])


removed_differences_by_covg_islands=find_differences_due_to_filtering(counts,failed_all,failed_all_but_island)
print(len(np.where( np.sum(removed_differences_by_covg_islands[0], axis=1) > 0 )[0]))
print(len(np.where( np.sum(removed_differences_by_covg_islands[1], axis=1) > 0 )[0]))
removed_exclusively_by_covg_island=np.setdiff1d(np.where( np.sum(removed_differences_by_covg_islands[1], axis=1) > 0 )[0],np.where( np.sum(removed_differences_by_covg_islands[0], axis=1) > 0 )[0])
for position,line in zip(p[np.setdiff1d(removed_exclusively_by_covg_island,singletons_to_mask)],removed_differences_by_covg_islands[1][np.setdiff1d(removed_exclusively_by_covg_island,singletons_to_mask)]):
    print(position,sampleNames[line])



# %%
"""
# POSSIBLE ADJUSTMENT TO GENERATE SINGLETONS TO FILTER PRIOR TO ANY ADDTIOINAL FILTERING
if blast_prior_to_filtering:
    [maf, maNT, minorNT, minorAF] = apy.div_major_allele_freq(counts) 
    calls = maNT
    hasmutation_no_filtering = (calls != refnti_m) & (calls < 4) 

    hasmutation_no_filtering[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
    goodpos_no_filtering = np.where( np.sum(hasmutation_no_filtering, axis=1) > 0 )[0]
    print(goodpos_no_filtering.size,'unfiltered goodpos found.')


    calls_for_treei = calls[goodpos_no_filtering,:]
    calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation
    refgenome_nts_for_tree=refnt[goodpos_no_filtering]
    calls_for_tree_ref_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree),axis=1)

    treesampleNamesLong = sampleNames # include all samples 

    # sampleNamesDnapars : max 10c. Use numeric with 10c (works for up to 10^10 samples! )
    sampleNamesDnapars = [ "{:010d}".format(i) for i in range(len(treesampleNamesLong))] # includes all samples

    # translate index to nucleotide

    ## outgroup already in samples, so unnecessary (CO92)
    refgenome_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));
    refgenome_nts_for_tree=refnt[goodpos_no_filtering]
    calls_for_tree_ref_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree),axis=1)

    treesampleNamesLong_ref_outgroup = np.append(['Sref'],treesampleNamesLong)

    sampleNamesDnapars_ref_outgroup=np.append(['Sref'],sampleNamesDnapars)

    apy.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,'unfiltered_tree')
    # run tree generation in logan
    apyp.create_ancestral_reconstruction(f'unfiltered_tree.raxml.bestTree',treesampleNamesLong_ref_outgroup,calls_for_tree_ref_outgroup,f'unfiltered_tree_initial_ancestral_reconstruction')
    ancestral_reconstruction_tree=f'unfiltered_tree_initial_ancestral_reconstruction/annotated_tree.nexus'
    ancestral_reconstruction_fasta=f"unfiltered_tree_initial_ancestral_reconstruction/ancestral_sequences.fasta"

    ancient_singletons = np.intersect1d(goodpos_no_filtering,np.where((np.sum(hasmutation_no_filtering[:,np.in1d(sampleNames,ancient_sample_names)],axis=1)==1) & (np.sum(hasmutation_no_filtering[:,~np.in1d(sampleNames,ancient_sample_names)],axis=1)==0))[0])

    # intersection with removed genomic islands or removed covg. percentile:
    ancient_singletons=np.concatenate([np.intersect1d(ancient_singletons,removed_genomic_islands),np.intersect1d(ancient_singletons,removed_coverage_percentile)])
    

    for pos in ancient_singletons:
        samples=[x for x in sampleNames[hasmutation_no_filtering[p==p[pos]][0]] if x in ancient_sample_names]
        print(p[pos], samples)
        

    sample_pos_to_output={}
    for pos in ancient_singletons:
        samples=[x for x in sampleNames[hasmutation_no_filtering[p==p[pos]][0]] if x in ancient_sample_names][0]
        if samples in sample_pos_to_output:
            sample_pos_to_output[samples].append(pos)
        else:
            sample_pos_to_output[samples] = [pos]

    os.makedirs('to_blast_tree_singletons',exist_ok=True)
    seq_records_to_write=[]
    for sample in sample_pos_to_output:
        sample_with_mutations=main_map_sample_nt_mutations_onto_reference(ancestral_reconstruction_tree,ancestral_reconstruction_fasta,sample,sampleNames,p[goodpos_no_filtering],calls[goodpos_no_filtering],ref_genome_folder)
        if len(sample_with_mutations)==0:
            pass
        else:
            for pos in sample_pos_to_output[sample]:
                chrom,position_on_chrom=apyp.p2chrpos([p[pos]],chrStarts)[0]
                true_chrom=scafNames[chrom-1]
                true_position_on_chrom_1_based=position_on_chrom+1
                seq_records_to_write.append(SeqRecord.SeqRecord(Seq.Seq(''.join(sample_with_mutations[chrom-1][position_on_chrom-25:position_on_chrom+25])),id=f'{true_chrom}:{true_position_on_chrom_1_based}:{sample}'))
    SeqIO.write(seq_records_to_write,handle=f'to_blast_tree_singletons/singletons_pre_filtering.fasta',format='fasta')



    blast_results_json_file='/Users/ad_loris/Downloads/17M8A76P01R-Alignment.json'
    with open(blast_results_json_file) as json_data:
        blast_results=json.load(json_data)

    pestis_taxids=[1035377,1234662,1345701,1345702,1345703,1345704,1345705,1345706,1345707,1345708,1345710,1455696,187410,214092,229193,349746,360102,377628,386656,547048,632,637382,637386,649716,748678]
    singletons_to_mask=[]
    singletons_to_mask_sample=[]
    for result in blast_results["BlastOutput2"]:
        query_title=result['report']['results']['search']['query_title']
        query_chrom,query_pos_one_based,query_sample=query_title.split(':')
        hits=result['report']['results']['search']['hits']
        max_identity=0
        taxids=[]
        for hit in hits:
            if max_identity==0: 
                max_bitscore=hit['hsps'][0]['bit_score']
                max_identity=hit['hsps'][0]['identity']/50
            if max_identity == hit['hsps'][0]['identity']/50 and max_bitscore == hit['hsps'][0]['bit_score']:
                if 'taxid' in hit['description'][0]:
                    taxids.append(hit['description'][0]['taxid'])
        intersection=np.intersect1d(np.array(taxids),np.array(pestis_taxids))
        if intersection.size == 0:
            p_to_mask = apyp.chrpos2p(query_chrom,int(query_pos_one_based)-1,scafNames,chrStarts)
            index_p_to_mask=np.where(p == p_to_mask)[0]
            singletons_to_mask.append(index_p_to_mask)
            singletons_to_mask_sample.append(query_sample)
    singletons_to_mask = np.array(singletons_to_mask)
    singletons_to_mask_sample = np.array(singletons_to_mask_sample)

    masked_ancient_singletons=set()
    for pos_array,sample_to_mask in zip(singletons_to_mask,singletons_to_mask_sample):
        pos=pos_array[0]
        sample_index=sampleNames==sample_to_mask
        overlapping_snps_25bp_either_side=np.where((p-p[pos]>-25) & (p-p[pos]<25))[0]
        print("P indices removed for sample:", overlapping_snps_25bp_either_side, sampleNames[sample_index])
        for position_masked in overlapping_snps_25bp_either_side:       
            calls[ position_masked,sample_index] = 4
            masked_ancient_singletons.add(position_masked)
            
    # IDENTIFY MANUAL INSPECTION SNPs 
    # ie - those in the same sample which are physcially <100bp away
    ancient_singletons = np.intersect1d(candpos,np.where((np.sum(hasmutation[:,np.in1d(sampleNames,ancient_sample_names)],axis=1)==1) & (np.sum(hasmutation[:,~np.in1d(sampleNames,ancient_sample_names)],axis=1)==0))[0])
    dist_last=0
    for i,pos in enumerate(ancient_singletons):
        if pos not in np.array(list(masked_ancient_singletons)):
            dist_this=p[pos]-dist_last
            if dist_this < 100:
                samples=[x for x in sampleNames[hasmutation[p==p[pos]][0]] if x in ancient_sample_names]
                samples_prev=[x for x in sampleNames[hasmutation[p==p[ancient_singletons[i-1]]][0]] if x in ancient_sample_names]
                if len(np.intersect1d(np.array(samples_prev),np.array(samples)))>0:
                    print(p[ancient_singletons[i-1]],samples_prev,pos_prev)
                    print(p[pos], samples,pos)
            dist_last=p[pos]
            pos_prev=pos

"""



# prior tree based haplotype blasting
def parse_zero_coverage(samplename):
    file_path=f'/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/production_run/bed_files/zerocoverage_regions/{samplename}.zerocoverage.gz'
    output_dict={}
    with gzip.open(file_path, 'rt') as f:
        for l in f:
            entries = l.split('\t')
            if entries[0] not in output_dict:
                output_dict[entries[0]] = set()
            for i in range(int(entries[1]),int(entries[2])):
                output_dict[entries[0]].add(i)
    return output_dict

os.makedirs('to_blast_tree_singletons',exist_ok=True)
seq_records_to_write=[]
for sample in sample_pos_to_output:
    sample_with_mutations=main_map_sample_nt_mutations_onto_reference(ancestral_reconstruction_tree,ancestral_reconstruction_fasta,sample,sampleNames,p[candpos],calls[candpos],ref_genome_folder)
    positions_with_zero_coverage=parse_zero_coverage(sample)
    if len(sample_with_mutations)==0:
        pass
    else:
        for pos in sample_pos_to_output[sample]:
            chrom,position_on_chrom=apyp.p2chrpos([p[pos]],chrStarts)[0]
            true_chrom=scafNames[chrom-1]
            true_position_on_chrom_1_based=position_on_chrom+1
            left_extent=min([x for x in range(position_on_chrom-25,position_on_chrom) if x not in positions_with_zero_coverage])
            right_extent=max([x for x in range(position_on_chrom,position_on_chrom+26) if x not in positions_with_zero_coverage])
            seq_records_to_write.append(SeqRecord.SeqRecord(Seq.Seq(''.join(sample_with_mutations[chrom-1][left_extent:right_extent])),id=f'{true_chrom}:{true_position_on_chrom_1_based}:{sample}'))

max_seqs_for_blast_file=100
for i in range(math.ceil(len(seq_records_to_write)/max_seqs_for_blast_file)):
    if (i+1)*100 > len(seq_records_to_write):
         this_chunk_to_output=seq_records_to_write[i*100:]   
    else:
        this_chunk_to_output=seq_records_to_write[i*100:(i+1)*100]   
    print(len(this_chunk_to_output))
    SeqIO.write(this_chunk_to_output,handle=f'to_blast_tree_singletons/final_concatenated_to_blast_with_samplename_{i}.fasta',format='fasta')

