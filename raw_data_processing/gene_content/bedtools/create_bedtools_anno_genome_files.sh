#! /bin/bash
# This script creates necessary input files for bedtools to annotate breadth and depth on genic regions

# generated from genome.gff for pestis reference genome (CO92)

working_directory=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1

awk -F'\t' '{print $1, "1", $2}' ${working_directory}/genome.fasta.fai | sed 's/ /\t/g' > ${working_directory}/genome_accession_for_bedtools.txt
awk -F'\t' '{print $1, $2}' ${working_directory}/genome.fasta.fai | sed 's/ /\t/g' > ${working_directory}/genome_file_for_bedtools.txt

# generated from genome.gff for pseudotb reference genome (IP 32953)

working_directory=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365

zcat GCF_000834295.1_ASM83429v1_genomic.gff.gz | grep -E '[[:space:]]RefSeq[[:space:]]gene' | cut -f1,4,5 > ${working_directory}/genic_accession_for_bedtools.txt
awk -F'\t' '{print $1, $2}' ${working_directory}/genome.fasta.fai | sed 's/ /\t/g' > ${working_directory}/genome_file_for_bedtools.txt

# generated from aida Stable4
working_directory=/u/iclight/pestis_work
working_directory_refs=/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1

cut -d'|' -f4,5 ${working_directory}/co92_virulence_genes.tsv | tr '|' ' ' | cut -f1,2,3 | sed 's/ //g' | tail -n+2 > ${working_directory}/co92_virulence_genes_for_bedtools.tsv
cp ${working_directory}/co92_virulence_genes_for_bedtools.tsv ${working_directory_refs}/co92_virulence_genes_for_bedtools.tsv 

# generate mapping from bedtools to gene name
cut -d'|' -f4,5,6 ${working_directory}/co92_virulence_genes.tsv | tr '|' ' ' | cut -f1,2,3,4 | sed 's/ //g' | tail -n+2 > ${working_directory}/co92_virulence_genes_for_bedtools_mapping.tsv
