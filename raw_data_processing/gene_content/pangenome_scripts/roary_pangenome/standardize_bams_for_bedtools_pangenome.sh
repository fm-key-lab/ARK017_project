#!/bin/bash -l

# output and error
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
# Job Name:
#SBATCH -J standardize_bedtoools_merging
#
# Directory:
#SBATCH -D /ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run_Ypestis_pseudotb_pangenome/results/bedtools
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

tail -n +2 ~/eager/pestis_project/ancient_pestis.tsv | cut -f1 | sort | uniq | while read sample_name
    do
    find . -wholename "./${sample_name}*.breadth.gz" > tmp.txt
    if [[ $(wc -l tmp.txt) -gt 1 ]] 
        then
        mkdir temp_for_merging
        grep $sample_name ~/eager/pestis_project/ancient_pestis.tsv | cut -f2 | while read ERR_acc
        do
            samtools view -F4 ../deduplication/${ERR_acc}/${ERR_acc}_rmdup.bam > temp_for_merging/${ERR_acc}_rmdup_mapped.bam
        done
        samtools merge -o temp_for_merging/${sample_name}_rmdup_mapped_merged.bam temp_for_merging/*.bam 
        bam=temp_for_merging/${sample_name}_rmdup_mapped_merged.bam
        bedtools coverage -nonamecheck -g genome.txt -sorted -a ${anno_file} -b ${bam} > ${sample_name}_all_libmerged.breadth
        gzip ${sample_name}_all_libmerged.breadth
        rm -rf temp_for_merging
    fi
done


anno_file=genic_accession_for_bedtools.txt
ls -1 ERR* | while read file_name
    do 
    basename=$(echo $file_name | cut -f1 -d'_')
    sample_name=$(grep $basename ~/eager/pestis_project/ancient_pestis.tsv | cut -f1)
    grep ${sample_name} ~/eager/pestis_project/ancient_pestis.tsv > tmp.tsv
    if [[ $(wc -l tmp.tsv) -gt 1 ]] 
        then
        mkdir temp_for_merging
        cut -f2 tmp.tsv | while read ERR_acc
        do
            samtools view -F4 ../deduplication/${ERR_acc}/${ERR_acc}_rmdup.bam > temp_for_merging/${ERR_acc}_rmdup_mapped.bam
        done
        samtools merge -o temp_for_merging/${sample_name}_rmdup_mapped_merged.bam temp_for_merging/*.bam 
        bam=temp_for_merging/${sample_name}_rmdup_mapped_merged.bam
        bedtools coverage -nonamecheck -g genome.txt -sorted -a ${anno_file} -b ${bam} > ${sample_name}_all_libmerged.breadth
        gzip ${sample_name}_all_libmerged.breadth
        rm -rf temp_for_merging
    else
        cp $file_name ${sample_name}_all_libmerged.breadth
        gzip ${sample_name}_all_libmerged.breadth
    fi
done


