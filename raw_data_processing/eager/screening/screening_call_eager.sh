#! /bin/bash
# execution call for screening of ARK data
# ran within a screen session, in parallel across sets of 2 samples into MALT database


tsv_loc=/u/iclight/eager/screening_data
tsv_file_main=ark_screening.tsv
num_per_loop=2
## common for all runs
name=$(echo ${tsv_file_main} | sed 's/\.tsv//g')
total_file_len=$(wc -l <${tsv_loc}/${tsv_file_main})
num_loops_needed=$(( $(( ${total_file_len} + ${num_per_loop} - 1)) / ${num_per_loop}))

## running the loop
for i in $(seq 1 ${num_loops_needed}) 
    do
    ## pre_nextflow things
    ## nextflow variables
    name_of_run=${name}_step_${i}
    ## create tsv
    start_i=$(($((${i}-1))*${num_per_loop}))
    end_i=$((${i}*${num_per_loop}))

    ## change to new folder in /ptmp and start the run
    mkdir -p /ptmp/iclight/eager/screening/${name}/${i} && cd /ptmp/iclight/eager/screening/${name}/${i} || exit
    
    run_loc=/ptmp/iclight/eager/screening/${name}/${i}
    cat ${tsv_loc}/intermediate_files/header.tsv > ${run_loc}/tsv.tmp.final.tsv
    awk -v start_i="$start_i" -v end_i="$end_i" 'NR > start_i && NR <= end_i {print}' ${tsv_loc}/${tsv_file_main} >> ${run_loc}/tsv.tmp.final.tsv

    ## call nextflow
    nextflow run nf-core/eager -bg -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_max \
    -name ${name_of_run} \
    --input './tsv.tmp.final.tsv' \
    --fasta '/u/iclight/eager/reference_genome/GRCh38.p13_genomic.fna' \
    --bwa_index '/u/iclight/eager/reference_genome/' \
    --fasta_index '/u/iclight/eager/reference_genome/GRCh38.p13_genomic.fna.fai' \
    --seq_dict '/u/iclight/eager/reference_genome/GRCh38.p13_genomic.fna.dict' \
    --outdir './results' \
    -w './work' \
    --skip_adapterremoval \
    --complexity_filter_poly_g \
    --run_bam_filtering \
    --bam_unmapped_type 'fastq' \
    --run_metagenomic_screening \
    --metagenomic_tool 'malt' \
    --database '/ptmp/iclight/malt_build_out_step4_2022_clustering' \
    --malt_min_support_mode 'reads' \
    --run_maltextract \
    --maltextract_taxon_list '/u/iclight/bin/HOPS/Resources/default_list_final.txt' \
    --maltextract_ncbifiles '/u/iclight/bin/HOPS/Resources' \
    --maltextract_destackingoff

    sleep 10m
done

# maltextract_taxon_list can be found in this directory.