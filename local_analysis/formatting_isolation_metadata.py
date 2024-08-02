# %%
# ## import libraries
import os
import numpy as np
import pandas as pd

## move to metadata folder

os.chdir('/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/metadata_tables')


# %%
# load metadata to format and collapse
strain_metadata_cui_path='associated_study_metadata/cui_2013_strains.csv'
strain_metadata_cui=pd.read_csv(strain_metadata_cui_path, header=1, sep=';').iloc[:133]
strain_metadata_cui['sampleName'] = strain_metadata_cui['Node'] + '_' + strain_metadata_cui['Strain ID']

cui_misformatted={
'0.PE4Cb_M0000002':'0.PE4_M0000002',
'0.ANT3e_42091':'0.ANT3e_42091b',
'1.IN2m_D1964002':'1.IN2m_D1964002b',
'1.IN3c_CMCC84038':'1.IN3c_CMCC84038b',
'1.ORI2c_YN2551':'1.ORI2c_YN2551b',
'1.ORI2i_CMCCK100001':'1.ORI2i_CMCCK100001a',
'1.ORI2b_CMCCK110001':'1.ORI2i_CMCCK110001b',
'2.MED3b_CMCC125002':'2.MED3b_CMCC125002b',
'3.ANT1a_7':'3.ANT1a_7b'}

cui_misformatted_corrected=[]
for s in strain_metadata_cui['sampleName']:
    if s in cui_misformatted:
        cui_misformatted_corrected.append(cui_misformatted[s])
    else:
        cui_misformatted_corrected.append(s)
strain_metadata_cui['sampleName']=cui_misformatted_corrected

foci_code_kislichkina_path='associated_study_metadata/kislichkina_2015_table1_foci_code.tsv'
foci_conversion_kislichkina=pd.read_csv(foci_code_kislichkina_path,sep='\t')

foci_conversion_anisimov_path='associated_study_metadata/anisimov_2004_table1_foci_formatted_with_m_labels.csv'
foci_conversion_anisimov=pd.read_csv(foci_conversion_anisimov_path)

translated_foci=[]
for foci_id in foci_conversion_kislichkina.foci_code_anisimov:
    transated_foci_from_id=foci_conversion_anisimov.query(f'focus_id == "{foci_id}"').main_host.iloc[0]
    translated_foci.append(transated_foci_from_id)
foci_conversion_kislichkina['Host/vector']=translated_foci

samples_to_isolation=pd.concat([strain_metadata_cui[['sampleName', 'Host/vector']],
        foci_conversion_kislichkina[['sampleName', 'Host/vector']]])

# %% export for later use
samples_to_isolation.to_csv('samples_to_isolation_source.tsv',sep='\t',header=True,index=False)
# %%
