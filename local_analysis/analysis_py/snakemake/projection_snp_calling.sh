#!/bin/bash -l
# specify the indexes (max. 30000) of the job array elements (max. 300 - the default job submit limit per user)
#SBATCH --array=1-201
# Standard output and error:
#SBATCH -o generate_merged_vcf_array_out/job_%A_%a.out        # Standard output, %A = job ID, %a = job array index
#SBATCH -e generate_merged_vcf_array_out/job_%A_%a.err        # Standard error, %A = job ID, %a = job array index
# Initial working directory:
#SBATCH -D /u/iclight/pestis_work
# Job Name:
#SBATCH -J extract_pestis_mapped_reads
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=50000
#
# Wall clock limit (max. is 24 hours):
#SBATCH --time=12:00:00

# Load compiler and MPI modules (must be the same as used for compiling the code)

conda activate /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/7dddcb5f5ec23945954fb8c789584cc9


sampleid=RISE386
bam_dir=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_mapping_rise386/results/trimmed_bam

vcf_dir=/ptmp/iclight/pestis_bams_clean_production_run_07_02_2024/vcf/projection

pileup=${vcf_dir}/${sampleid}_aligned.sorted.pileup
variants=${vcf_dir}/${sampleid}_aligned.sorted.strain.variant.vcf.gz
vcf_strain=${vcf_dir}/${sampleid}_aligned.sorted.strain.vcf.gz
vcf_raw=${vcf_dir}/${sampleid}_aligned.sorted.strain.gz

ref=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta
bam=${bam_dir}/${sampleid}_libmerged.trimmed.bam

samtools mpileup -q30 -x -s -O -d3000 -f ${ref} ${bam} > ${pileup}
samtools mpileup -q30 -t SP -d3000 -vf ${ref} ${bam} > ${vcf_raw}
bcftools call -c -Oz -o ${vcf_strain} ${vcf_raw}
bcftools view -Oz -v snps -q .9 ${vcf_strain} > ${variants}
tabix -p vcf ${variants}



# # # 

sampleid=C90
bam_dir=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_mapping_c90/results/deduplication/C90

vcf_dir=/ptmp/iclight/pestis_bams_clean_production_run_07_02_2024/vcf/projection

pileup=${vcf_dir}/${sampleid}_aligned.sorted.pileup
variants=${vcf_dir}/${sampleid}_aligned.sorted.strain.variant.vcf.gz
vcf_strain=${vcf_dir}/${sampleid}_aligned.sorted.strain.vcf.gz
vcf_raw=${vcf_dir}/${sampleid}_aligned.sorted.strain.gz

ref=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta
bam=${bam_dir}/${sampleid}_rmdup.bam

samtools mpileup -q30 -x -s -O -d3000 -f ${ref} ${bam} > ${pileup}
samtools mpileup -q30 -t SP -d3000 -vf ${ref} ${bam} > ${vcf_raw}
bcftools call -c -Oz -o ${vcf_strain} ${vcf_raw}
bcftools view -Oz -v snps -q .9 ${vcf_strain} > ${variants}
tabix -p vcf ${variants}
