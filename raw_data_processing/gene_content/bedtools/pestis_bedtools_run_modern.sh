#!/bin/bash -l

# output and error
#SBATCH -o /ptmp/iclight/pestis_work/gene_content/array_out/job_%A.out
#SBATCH -e /ptmp/iclight/pestis_work/gene_content/array_out/job_%A.err

# Job Name:
#SBATCH -J standardize_bams_for_bedtools_pangenome_modern
#
# Directory:
#SBATCH -D /ptmp/iclight/pestis_work/gene_content
#
# Node feature:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --mem=500000
#SBATCH --mail-type=fail
#SBATCH --mail-type=time_limit
#SBATCH --mail-user=light@mpiib-berlin.mpg.de
#
# wall clock limit:
#SBATCH --time=12:00:00

conda activate basic_genomics

# relevant paths for I/O
dir_of_run=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_modern_production_run
dir_for_output=/ptmp/iclight/pestis_work/gene_content/pestis_virulence_genes/
anno_file=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/co92_virulence_genes_for_bedtools.tsv
genome_file=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome_file_for_bedtools.txt

tsv_file=/u/iclight/eager/pestis_project/modern_pestis.tsv

# create output dir (if not present already)
mkdir -p ${dir_for_output}

# process genomes
for SLURM_ARRAY_TASK_ID in $(seq 1 135)
    do
    sample_name=$(tail -n +2 ${tsv_file} | cut -f1 | sort | uniq | sed -n "${SLURM_ARRAY_TASK_ID}p")
    mkdir -p temp_for_merging/${SLURM_ARRAY_TASK_ID}_modern
    grep ${sample_name} ${tsv_file} | cut -f2 | while read ERR_acc
    do
        ln -s ${dir_of_run}/results/deduplication/${ERR_acc}/${ERR_acc}_rmdup.bam temp_for_merging/${SLURM_ARRAY_TASK_ID}_modern/${ERR_acc}_rmdup_mapped.bam 
    done

    # run for all reads (no MQ filter)
    bam=/ptmp/iclight/pestis_bams_clean_production_run_07_02_2024/${sample_name}.bam
    samtools merge -@ 72 -o ${bam} temp_for_merging/${SLURM_ARRAY_TASK_ID}_modern/*.bam 
    bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} > ${dir_for_output}/${sample_name}.breadth
    bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} -mean > ${dir_for_output}/${sample_name}.depth
    gzip ${dir_for_output}/${sample_name}.breadth
    gzip ${dir_for_output}/${sample_name}.depth

    rm -rf temp_for_merging/${SLURM_ARRAY_TASK_ID}_modern
done