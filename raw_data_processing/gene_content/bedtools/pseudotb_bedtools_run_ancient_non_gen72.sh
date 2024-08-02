#!/bin/bash -l
#SBATCH --array=1-66

# output and error
#SBATCH -o /ptmp/iclight/pestis_work/gene_content/array_output/job_%A_%a.out
#SBATCH -e /ptmp/iclight/pestis_work/gene_content/array_output/job_%A_%a.err

# Job Name:
#SBATCH -J pseudotb_bedtools_run_ancient_non_gen72
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

dir_of_run=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run_pseudotb
dir_for_output=/ptmp/iclight/pestis_work/gene_content/pseudotuberculosis_bedtools_output/
anno_file=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genic_accession_for_bedtools.txt
genome_file=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome_file_for_bedtools.txt
tsv_file=/u/iclight/eager/pestis_project/ancient_pestis.tsv
ancient_sample_basename_file=/u/iclight/eager/pestis_project/ancient_sample_basename_nongen72_pestis.tsv

# create output dir (if not present already)
mkdir -p ${dir_for_output}

# process genomes
sample_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${ancient_sample_basename_file})
mkdir -p temp_for_merging/${SLURM_ARRAY_TASK_ID}
grep $sample_name ${tsv_file} | cut -f2 | while read ERR_acc
do
    ln -s ${dir_of_run}/results/deduplication/${ERR_acc}/${ERR_acc}_rmdup.bam temp_for_merging/${SLURM_ARRAY_TASK_ID}/${ERR_acc}_rmdup_mapped.bam 
done

# run for all reads (no MQ filter)
bam_merged=temp_for_merging/${SLURM_ARRAY_TASK_ID}/${sample_name}_rmdup_merged.bam
bam=/ptmp/iclight/pestis_bams_pseudotuberculosis_clean_production_run_16_04_2024/${sample_name}.bam

samtools merge -@ 72 -o ${bam_merged} temp_for_merging/${SLURM_ARRAY_TASK_ID}/*.bam 
samtools view -@ 72 -F4 -b ${bam_merged} > ${bam}
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} > ${dir_for_output}/${sample_name}.breadth
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} -mean > ${dir_for_output}/${sample_name}.depth
gzip ${dir_for_output}/${sample_name}.breadth
gzip ${dir_for_output}/${sample_name}.depth

rm -rf temp_for_merging/${SLURM_ARRAY_TASK_ID}

