#!/bin/bash -l

# output and error
#SBATCH -o /ptmp/iclight/pestis_work/standardize_bams_for_bedtools_pangenome/array/job_%A.out
#SBATCH -e /ptmp/iclight/pestis_work/standardize_bams_for_bedtools_pangenome/array/job_%A.err

# Job Name:
#SBATCH -J standardize_bams_for_bedtools_pangenome_gen72
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

dir_of_run=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_GEN72_production_run_Ypestis_pseudotb_enterocolitica_pangenome/
dir_for_output=${dir_of_run}/results/bedtools_unified
anno_file=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_pseudotb_enterocolitica_pangenome/genic_accession_for_bedtools.txt
genome_file=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_pseudotb_enterocolitica_pangenome/genome_file_for_bedtools.txt
tsv_file=/u/iclight/eager/pestis_project/ancient_samples_eager_input_GEN72.tsv

gen72_dir=gen72

mkdir -p ${dir_for_output}
sample_name=GEN72
mkdir -p temp_for_merging/${gen72_dir}
grep $sample_name ${tsv_file} | cut -f2 | while read ERR_acc
do
    ln -s ${dir_of_run}/results/deduplication/${ERR_acc}/${ERR_acc}_rmdup.bam temp_for_merging/${gen72_dir}/${ERR_acc}_rmdup_mapped.bam 
done

# run for all reads (no MQ filter)
bam=temp_for_merging/${gen72_dir}/${sample_name}_rmdup_mapped_merged.bam
samtools merge -@ 72 -o ${bam} temp_for_merging/${gen72_dir}/*.bam 
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} > ${dir_for_output}/${sample_name}.breadth
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} -mean > ${dir_for_output}/${sample_name}.depth
gzip ${dir_for_output}/${sample_name}.breadth
gzip ${dir_for_output}/${sample_name}.depth

# run for only MQ20 reads (remove all double-mapped reads)
bam=temp_for_merging/${gen72_dir}/${sample_name}_rmdup_mapped_qual20_merged.bam
samtools view -@ 72 -F4 -q 20 -b temp_for_merging/${gen72_dir}/${sample_name}_rmdup_mapped_merged.bam > ${bam}
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} > ${dir_for_output}/${sample_name}_qual20.breadth
bedtools coverage -nonamecheck -g ${genome_file} -sorted -a ${anno_file} -b ${bam} -mean > ${dir_for_output}/${sample_name}_qual20.depth
gzip ${dir_for_output}/${sample_name}_qual20.breadth
gzip ${dir_for_output}/${sample_name}_qual20.depth

rm -rf temp_for_merging/${gen72_dir}