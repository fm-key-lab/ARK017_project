#! /bin/bash

# pseudotb
fetchngs_dir=/ptmp/iclight/pestis_data_fetchngs_pseudotb
mkdir -p ${fetchngs_dir} && cd ${fetchngs_dir}/  

grep 'Y. pseudotuberculosis' modern_accessions.tsv | cut -f5 > input.csv

nextflow run nf-core/fetchngs --input input.csv --force_sratools_download -r 1.11.0 --outdir . -profile singularity
cd metadata

grep 'Y. pseudotuberculosis' ../modern_accessions.tsv | while read line
    do
    id=$(echo $line | cut -f1 -d' ')
    accession=$(echo $line | cut -f6 -d' '| sed 's/ /_/g')
    fasta_basename=$(tail -n+2 ${accession}.runinfo_ftp.tsv | cut -f1)
    echo "${id} ${id} 1 4 PE   ${fetchngs_dir}/fastq/${fasta_basename}_1.fastq.gz ${fetchngs_dir}/fastq/${fasta_basename}_2.fastq.gz NA" >> /u/iclight/eager/pestis_project/pseudotb.tsv
done


# pestis
fetchngs_dir=/ptmp/iclight/pestis_data_fetchngs
mkdir -p ${fetchngs_dir} && cd ${fetchngs_dir}/  

grep 'Y. pestis' modern_accessions.tsv | cut -f5 > input.csv

nextflow run nf-core/fetchngs --input input.csv --force_sratools_download -r 1.11.0 --outdir . -profile singularity
cd metadata

grep 'Y. pseudotuberculosis' ../modern_accessions.tsv | while read line
    do
    id=$(echo $line | cut -f1 -d' ')
    accession=$(echo $line | cut -f6 -d' '| sed 's/ /_/g')
    fasta_basename=$(tail -n+2 ${accession}.runinfo_ftp.tsv | cut -f1)
    echo "${id} ${id} 1 4 PE unknown double none ${fetchngs_dir}/fastq/${fasta_basename}_1.fastq.gz ${fetchngs_dir}/fastq/${fasta_basename}_2.fastq.gz NA" >> /u/iclight/eager/pestis_project/pestis.tsv
done


cat /u/iclight/eager/other_input_sheets/header.tsv  /u/iclight/eager/pestis_project/pestis.tsv /u/iclight/eager/pestis_project/pseudotb.tsv > /u/iclight/eager/pestis_project/all_modern_samples.tsv

cat /u/iclight/eager/other_input_sheets/header.tsv /u/iclight/eager/pestis_project/pestis.tsv > /u/iclight/eager/pestis_project/modern_pestis.tsv
grep -E 'YAC' /u/iclight/eager/pestis_project/pseudotb.tsv >> /u/iclight/eager/pestis_project/modern_pestis.tsv