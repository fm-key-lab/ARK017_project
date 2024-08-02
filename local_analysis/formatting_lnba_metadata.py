# %%
# ## import libraries
import os
import numpy as np
import pandas as pd

## move to metadata folder

os.chdir('/Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/ark017_lnba_project/metadata_tables')


# gather metadata 
lnba_metadata_path='associated_study_metadata/Metadata_coordinates_dating_sex.csv'
lnba_metadata=pd.read_csv(lnba_metadata_path)

# fix col_names, types
lnba_metadata=lnba_metadata.set_axis(list(map(lambda x: x.replace(' ','_'),lnba_metadata.columns)),axis=1)
lnba_metadata = lnba_metadata.astype({'Median_C14': float})

# update aida publication info
lnba_metadata['Published']=lnba_metadata['Published'].replace('This study','Andrades Valtueña et al. (2020)')

#fix wrong names
aida_names_to_fix=sorted(np.array(lnba_metadata['Sample_Name'])[~np.in1d(lnba_metadata['Sample_Name'],sampleNames)])
our_names_to_fix_in_df=sorted(sampleNames[LNBA_clade_indices][~np.in1d(sampleNames[LNBA_clade_indices],lnba_metadata['Sample_Name'])])

our_names_to_fix_in_df=['1343UnTal85',
 '6Post',
 'Gyvakarai1']
lnba_metadata=lnba_metadata.replace(aida_names_to_fix,our_names_to_fix_in_df)

# add first Crete sample:
crete_sample_metadata=['HGC009','HGC009','Hagios Charalambos','HGC',35.1772505,25.4410963,3922,3944,3859,'M','C14','LNBA','Neumann et al. (2022)','None']
crete_sample=pd.DataFrame({y:x for x,y in zip(crete_sample_metadata,lnba_metadata.columns)},index=[0])

ARK_sample_metadata=['ARK017','ARK017','Arkaim','ARK',52.649261, 59.571443,(1900+1950+1778+1950)/2,1900+1950,1778+1950,'NA','C14','LNBA','This study','None']
ARK_sample=pd.DataFrame({y:x for x,y in zip(ARK_sample_metadata,lnba_metadata.columns)},index=[0])

C10098_sample_metadata=['C10098','C10098','Charterhouse Warren','C',51.299, -2.715,(4145+3910)/2,4145,3910,'?','C14','LNBA','Swali et al. (2023)','None']
C10098_sample=pd.DataFrame({y:x for x,y in zip(C10098_sample_metadata,lnba_metadata.columns)},index=[0])

DSH008_sample_metadata=['DSH008','DSH008','Drasenhofen','C',48.758, 16.651,(2128+1950+1931+1950)/2,2128+1950,1931+1950,'M','C14','LNBA','Neumann et al. (2023)','None']
DSH008_sample=pd.DataFrame({y:x for x,y in zip(DSH008_sample_metadata,lnba_metadata.columns)},index=[0])

DSH025_sample_metadata=['DSH025','DSH025','Drasenhofen','C',48.758, 16.651,(2026+1950+1884+1950)/2,2026+1950,1884+1950,'M','C14','LNBA','Neumann et al. (2023)','None']
DSH025_sample=pd.DataFrame({y:x for x,y in zip(DSH025_sample_metadata,lnba_metadata.columns)},index=[0])


lnba_metadata=pd.concat([ARK_sample,crete_sample,C10098_sample,DSH008_sample,DSH025_sample,lnba_metadata[:]]).reset_index(drop=True)


# Update to include median BCE date in re-calibration (from analysis tables directory)
# recalibarted C14 dates (BCE)
unified_calibration_path='c14_unified_processing/lnba_and_pre_unified_calibration_20.csv'
unified_calibration=pd.read_csv(unified_calibration_path,header=1)
unified_calibration['sample']=[x.split(' ')[1] for x in unified_calibration['Unnamed: 0']]

bce_dates_calibrated=[]
for x in lnba_metadata.Sample_Name:
    if x in ['GLZ002','GLZ001','RV2039']: # dates which had reservoir effect correction in original publication
        x=f'{x}_adj'
    print(x,unified_calibration.query(f'sample=="{x}"')['median'].iloc[0])
    bce_dates_calibrated.append(unified_calibration.query(f'sample=="{x}"')['median'].iloc[0])

lnba_metadata['median_bce_date_calibrated']=bce_dates_calibrated
# %% export for later use
lnba_metadata.to_csv('lnba_prelnba_metadata.tsv',sep='\t',index=False)


# %% Reduce metadata for just hexdec map
subset_lnba_metdata=lnba_metadata.query('Branch =="LNBA" & Dating == "C14"')

# %% export for later use
subset_lnba_metdata[['Sample_Name','Latitude','Longitude','median_bce_date_calibrated']].to_csv('lnba_samplenames_lat_long_age.tsv',sep='\t',index=False)

# %%
