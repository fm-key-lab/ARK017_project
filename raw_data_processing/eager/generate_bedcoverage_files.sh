#!/bin/bash -l

# output and error
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
# Job Name:
#SBATCH -J generate_bedcoverage_files
#
# Directory:
#SBATCH -D /u/iclight/pestis_work/pestis_evolution_bedcoverage_files
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
#SBATCH --time=24:00:00

echo "Creating bedtool coverage files for all pestis bams, histograms of overall coverage and 0 coverage regions."
echo "Merging all regions with 0 coverage that are within 500bp of each other, then computing proportions of these regions that have 0 coverage to filter within analysispy."
echo "Also creating a separate file which has all regions of 0 coverage"

bam_directory=/raven/ptmp/iclight/snakemake/pestis_production_run_minoraf_update_30_04_2024/mapping/data
conda activate basic_genomics

outdir=out_files_30_04_2024

mkdir ${outdir}

for file_name in $(ls -1 $bam_directory)
    do
    echo "start ${file_name}"
    touch ${bam_directory}/${file_name}/${file_name}*.bam*
    samtools view -b -F 4 -q 30 ${bam_directory}/${file_name}/${file_name}.bam > ${bam_directory}/${file_name}/${file_name}_filtered.bam

    bedtools genomecov -ibam ${bam_directory}/${file_name}/${file_name}_filtered.bam > ${outdir}/${file_name}_genome_coverage_hist.tsv
    bedtools genomecov -bga -ibam ${bam_directory}/${file_name}/${file_name}_filtered.bam > ${outdir}/${file_name}_genome_coverage_allpos.tsv 
    grep -w '0$' ${outdir}/${file_name}_genome_coverage_allpos.tsv > ${outdir}/${file_name}.zerocoverage
    bedtools merge -d 500 -i ${outdir}/${file_name}.zerocoverage  > tmp
    bedtools coverage -b ${bam_directory}/${file_name}/${file_name}_filtered.bam -a tmp -hist | awk -F'\t' '{if ($4 == 0) print}' > ${outdir}/${file_name}_merged_zero_covg_regions.tsv
    gzip ${file_name}*.tsv
    gzip ${outdir}/${file_name}.zerocoverage
    echo "end ${file_name}"
done



for file_name in $(ls -1 $bam_directory)
    do
    bedtools merge -d 500 -i ${outdir}/${file_name}.zerocoverage  > tmp
    bedtools coverage -b ${bam_directory}/${file_name}/${file_name}_filtered.bam -a tmp -hist | awk -F'\t' '{if ($4 == 0) print}' > ${outdir}/${file_name}_merged_zero_covg_regions.tsv
    gzip ${file_name}*.tsv
    gzip ${outdir}/${file_name}.zerocoverage
    echo "end ${file_name}"
done
