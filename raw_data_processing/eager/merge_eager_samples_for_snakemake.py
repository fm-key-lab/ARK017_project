import pandas as pd
import numpy as np
import subprocess

eager_directory_input_path='/u/iclight/eager/pestis_project'

path_to_gen72_eager_input=f'{eager_directory_input_path}/production_input_files/ancient_pestis_GEN72.tsv'
path_to_other_ancient_input=f'{eager_directory_input_path}/production_input_files/ancient_pestis_non_GEN72.tsv'

path_to_eager_results_gen72='/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_GEN72_production_run'
path_to_eager_results_other_ancient='/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run'

output_path='/ptmp/iclight/pestis_bams_clean_production_run_07_02_2024'

def get_bams_to_merge_helper(eager_input,eager_results_path,sample_basename):
    bam_paths_to_merge=[]
    indices_containing_sample_bool=eager_input.Sample_Name.str.contains(sample_basename)
    eager_input_this_sample=eager_input.loc[indices_containing_sample_bool==True]
    for i,row in eager_input_this_sample.iterrows():
        # get if any trimming would be done
        samplename_this_library=row.Library_Name
        if row.UDG_Treatment.isin(['half','none']):
            # check if path exists
            bam_paths_to_merge.append(f'{eager_results_path}/trimmed_bam/{samplename_this_library}.trimmed.bam')
        else:
            bam_paths_to_merge.append(f'{eager_results_path}/deduplication/{samplename_this_library}/{samplename_this_library}_rmdup.bam')
    return bam_paths_to_merge

def generate_sample_basenames_helper(eager_input):
    list_of_sample_basenames=['_'.join(x.spilt('_')[:-1]) for x in eager_input.Sample_Names]
    return np.unique(np.array(list_of_sample_basenames))

def main(eager_input_path,eager_results_path,output_path):
    eager_input=pd.read_csv(eager_input_path,header=0,sep='\t')
    sample_basenames=generate_sample_basenames_helper(eager_input)
    for sample_basename in sample_basenames:
        list_of_bam_paths=get_bams_to_merge_helper(eager_input,eager_results_path,sample_basename)
        formatted_list_of_bam_paths=' '.join(list_of_bam_paths)
        command=f'samtools merge -b {output_path}/{sample_basename}.bam {formatted_list_of_bam_paths}'
        print(f'Running merging for {sample_basename}')
        subprocess.run([command],shell=True)

if __name__ == "__main__":
    main(path_to_gen72_eager_input,path_to_eager_results_gen72,output_path)
    main(path_to_other_ancient_input,path_to_eager_results_other_ancient,output_path)
    