#!/bin/bash

# download and format reference genomes for pangenome generation pestis/pseudotuberculosis

cd /ptmp/iclight/pestis_work

cat /u/iclight/pestis_work/pestis_pseudotb_pangenome_genomes.tsv | awk -F '\t' '{ print $20 }' | parallel \
--jobs 16 \
' wget {}/$(basename {})_genomic.fna.gz {}/$(basename {})_genomic.gff.gz '

mkdir -p fastas

mv *.fna.gz fastas

cd fastas && gunzip *