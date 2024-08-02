#!/bin/bash -l
#SBATCH --array=1-101

# output and error
#SBATCH -o /u/iclight/pestis_work/ancient_singleton_validation/slurm/blast_execution_ancient_singletons_job_%A_%a.out
#SBATCH -e /u/iclight/pestis_work/ancient_singleton_validation/slurm/blast_execution_ancient_singletons_job_%A.err
#
# Job Name:
#SBATCH -J blast_execution_ancient_singletons
#
# Directory:
#SBATCH -D /u/iclight/pestis_work/ancient_singleton_validation

# Node feature:
#SBATCH --nodes=1
#SBATCH --mem=120000
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=fail
#SBATCH --mail-type=time_limit
#SBATCH --mail-user=light@mpiib-berlin.mpg.de
#
# wall clock limit:
#SBATCH --time=24:00:00

module load intel/21.2.0 impi/2021.2 mkl

conda activate basic_genomics

export BLASTDB='/u/iclight/bin/blastdb/nt'

# General I/O variables
outdir=blast_results
queries_dir=blast_queries
queries_info=singletons_to_extract_mappinginfo.txt
dir_of_bams=/ptmp/iclight/pestis_bams_clean_production_run_07_02_2024

# only run 11 lines of the query file (total lenght is just over 1000 singletons so 10*101 is enought)
start_i=$(($((${SLURM_ARRAY_TASK_ID}-1))*11))
end_i=$((${SLURM_ARRAY_TASK_ID}*11))

awk -v start_i="$start_i" -v end_i="$end_i" 'NR > start_i && NR <= end_i {print}' ${queries_info} | while read line
    do 
    # set up singleton specific I/O names 
    sample_name=$(echo $line | cut -d',' -f1)
    query_id=$(echo $line | cut -d',' -f2)
    query_file_name=${queries_dir}/${sample_name}_${query_id}.fasta
    query_output_name=${outdir}/${sample_name}_${query_id}_out.tsv

    # get all (q30) reads supports to singleton variant call, format to a fasta file
    samtools view -q30 ${dir_of_bams}/${sample_name}.bam ${query_id} | cut -f1,10 | sed 's/^/>/g' | sed 's/\t/\n/g' > ${query_file_name}
    
    # blastn all reads supporting singleton variant call
    blastn -task blastn-short -query ${query_file_name} -db nt -evalue 0.1 -num_threads 64 -outfmt "6 qacc sacc staxid bitscore pident" -out "${query_output_name}"
done
