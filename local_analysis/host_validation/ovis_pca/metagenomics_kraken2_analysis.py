# %% 
import numpy as np
import pandas as pd
import os
import gzip
import glob
import csv
# %% 
os.chdir('/Users/ad_loris/Documents/key_lab/outputs/pestis_evolution/taxonomic_classification')

# %%
def parse_gzipped_tsv_total_reads_assinged(file_path):
    results_dict={}
    with gzip.open(file_path, 'rt') as f:
        reader = csv.reader(f,  delimiter='\t')
        for row in reader:
            if row[7].strip() == 'unclassified':   
                results_dict['unclassified']=row[1]
            if row[7].strip() == 'root':
                results_dict['classified']=row[1]            
            if row[7].strip() == 'Bos':
                results_dict['Bos_classified']=row[1]
            if row[7].strip() == 'Ovis':
                results_dict['Ovis_classified']=row[1]
            if row[7].strip() == 'Capra':
                results_dict['Capra_classified']=row[1]
            if row[7].strip() == 'Yersinia pseudotuberculosis complex':
                results_dict['Yersinia_pseudotuberculosis_complex_classified']=row[1]
            if row[7].strip() == 'Yersinia pestis':
                results_dict['Yersinia_pestis_classified']=row[1]
            if row[7].strip() == 'Bacteria':
                results_dict['Bacteria_classified']=row[1]

    sum_classified_unclassified=int(results_dict['classified'])+int(results_dict['unclassified'])
    classified_items=[x for x in results_dict.keys() if x.endswith('classified') and not x.startswith('un')]
    for data_value in classified_items:
        if data_value == 'classified' or data_value.endswith('_classified'):
            results_dict[f'proportion_{data_value}']=int(results_dict[data_value])/sum_classified_unclassified
    return results_dict

# %% screening data
merged_results={}
for index,file_path in enumerate(glob.glob('screening_data_kraken_results/*.gz')):
    this_file_basename=file_path.split('/')[1].split('_se_')[0]
    merged_results[this_file_basename]=parse_gzipped_tsv_total_reads_assinged(file_path)
to_clean_for_output=pd.DataFrame(data=merged_results,dtype=float)
to_clean_for_output.fillna(0,inplace=True)
to_clean_for_output = to_clean_for_output.reindex(sorted(to_clean_for_output.columns), axis=1)

to_clean_for_output.to_csv('screening_data_kraken_summarized_results.tsv',sep='\t')
# %% capture data
merged_results={}
for index,file_path in enumerate(glob.glob('capture_data_kraken_results/*.gz')):
    this_file_basename=file_path.split('/')[1].split('.YP1')[0]
    merged_results[this_file_basename]=parse_gzipped_tsv_total_reads_assinged(file_path)
to_clean_for_output=pd.DataFrame(data=merged_results,dtype=float)
to_clean_for_output.fillna(0,inplace=True)
to_clean_for_output = to_clean_for_output.reindex(sorted(to_clean_for_output.columns), axis=1)

to_clean_for_output.to_csv('screening_data_kraken_summarized_results.tsv',sep='\t')
