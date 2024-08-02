#!/bin/bash -l
#SBATCH --array=1-65

# output and error
#SBATCH -o /ptmp/iclight/pestis_work/standardize_bams_for_bedtools_pangenome/array/job_%A_%a.out
#SBATCH -e /ptmp/iclight/pestis_work/standardize_bams_for_bedtools_pangenome/array/job_%A_%a.err

# Job Name:
#SBATCH -J standardize_bams_for_bedtools_pangenome_ancient_non_gen72
#
# Directory:
#SBATCH -D /ptmp/iclight/pestis_work/standardize_bams_for_bedtools_pangenome
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

dir_of_run=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run_Ypestis_pseudotb_enterocolitica_pangenome
dir_for_output=${dir_of_run}/results/bedtools_unified
anno_file=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_pseudotb_enterocolitica_pangenome/genic_accession_for_bedtools.txt
genome_file=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_pseudotb_enterocolitica_pangenome/genome_file_for_bedtools.txt
tsv_file=/u/iclight/eager/pestis_project/ancient_pestis.tsv

# create output dir (if not present already)
mkdir -p ${dir_for_output}

# process genomes
sample_name=$(tail -n +2 ${tsv_file} | cut -f1 | sort | uniq | sed -n "${SLURM_ARRAY_TASK_ID}p")
mkdir -p temp_for_merging/${SLURM_ARRAY_TASK_ID}
grep $sample_name ${tsv_file} | cut -f2 | while read ERR_acc
do
    ln -s ${dir_of_run}/results/deduplication/${ERR_acc}/${ERR_acc}_rmdup.bam temp_for_merging/${SLURM_ARRAY_TASK_ID}/${ERR_acc}_rmdup_mapped.bam 
done

# run for all reads (no MQ filter)
bam_merged=temp_for_merging/${SLURM_ARRAY_TASK_ID}/${sample_name}_rmdup_merged.bam
bam=temp_for_merging/${SLURM_ARRAY_TASK_ID}/${sample_name}_rmdup_mapped_merged.bam

samtools merge -@ 72 -o ${bam_merged} temp_for_merging/${SLURM_ARRAY_TASK_ID}/*.bam 
samtools view -@ 72 -F4 -b ${bam_merged} > ${bam}
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} > ${dir_for_output}/${sample_name}.breadth
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} -mean > ${dir_for_output}/${sample_name}.depth
gzip ${dir_for_output}/${sample_name}.breadth
gzip ${dir_for_output}/${sample_name}.depth

# run for only MQ20 reads (remove all double-mapped reads)
bam=temp_for_merging/${SLURM_ARRAY_TASK_ID}/${sample_name}_rmdup_mapped_qual20_merged.bam
samtools view -@ 72 -F4 -q 20 -b temp_for_merging/${SLURM_ARRAY_TASK_ID}/${sample_name}_rmdup_mapped_merged.bam > ${bam}

bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} > ${dir_for_output}/${sample_name}_qual20.breadth
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} -mean > ${dir_for_output}/${sample_name}_qual20.depth
gzip ${dir_for_output}/${sample_name}_qual20.breadth
gzip ${dir_for_output}/${sample_name}_qual20.depth

rm -rf temp_for_merging/${SLURM_ARRAY_TASK_ID}
