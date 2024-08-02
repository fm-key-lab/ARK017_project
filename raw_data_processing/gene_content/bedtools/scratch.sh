

tsv_file=/u/iclight/eager/pestis_project/ancient_pestis.tsv

ls -1 | while read dir_num
    do
    sample_name=$(tail -n +2 ${tsv_file} | cut -f1 | sort | uniq | sed -n "${dir_num}p")
    echo $sample_name
done
    bam_merged=temp_for_merging/${dir_num}/${sample_name}_rmdup_merged.bam
    bam=temp_for_merging/${dir_num}/${sample_name}_rmdup_mapped_merged.bam

    samtools view -@ 72 -F4 -q 20 -b ${bam_merged} > ${bam}

done



dir_of_run=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run
dir_for_output=${dir_of_run}/results/bedtools_unified
tsv_file=/u/iclight/eager/pestis_project/ancient_pestis.tsv
# create output dir (if not present already)
mkdir -p ${dir_for_output}

# process genomes
tail -n +2 ${tsv_file} | cut -f1 | sort | uniq | while read sample_name
    do
    mkdir -p temp_for_merging/${sample_name}
    grep $sample_name ${tsv_file} | cut -f2 | while read ERR_acc
        do
        ln -s ${dir_of_run}/results/deduplication/${ERR_acc}/${ERR_acc}_rmdup.bam temp_for_merging/${sample_name}/${ERR_acc}_rmdup_mapped.bam 
    done

    # run for all reads (no MQ filter)
    bam_merged=temp_for_merging/${sample_name}/${sample_name}_rmdup_merged.bam
    bam=temp_for_merging/${sample_name}/${sample_name}_rmdup_mapped_merged.bam

    samtools merge -@ 8 -o ${bam_merged} temp_for_merging/${sample_name}/*.bam 
    samtools view -@ 8 -F4 -b ${bam_merged} > ${bam}
done

dir_of_run=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run
dir_for_output=${dir_of_run}/results/bedtools_unified

ls -1 | while read dir_num
    do
    bam=${dir_num}/${dir_num}_rmdup_mapped_merged.bam
    cp ${bam} ${dir_of_run}/results/bedtools_unified/${dir_num}.bam
done





outdir=blast_results
queries_dir=blast_queries
queries_info=singletons_to_extract_mappinginfo.txt
dir_of_bams=/ptmp/iclight/pestis_bams_clean_production_run_07_02_2024

 cat ${queries_info} | while read line
    do 
    sample_name=$(echo $line | cut -d',' -f1)
    query_id=$(echo $line | cut -d',' -f2)
    query_file_name=${queries_dir}/${sample_name}_${query_id}.fasta
    query_output_name=${outdir}/${sample_name}_${query_id}_out.tsv
    samtools mpileup -r ${query_id} -q30 -x -s -O --output-QNAME -d3000 -f /nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta ${dir_of_bams}/${sample_name}.bam > mpileup_with_readnames/${sample_name}_${query_id}.pileup
done


# only run 101 lines of the query file (total lenght is just over 1000 singletons so 10*101 is enought)
start_i=$(($((${SLURM_ARRAY_TASK_ID}-1))*101))
end_i=$((${SLURM_ARRAY_TASK_ID}*101))

awk -v start_i="$start_i" -v end_i="$end_i" 'NR > start_i && NR <= end_i {print}' ${queries_info} | while read line
    do 
    # set up singleton specific I/O names 
    sample_name=$(echo $line | cut -d',' -f1)
    query_id=$(echo $line | cut -d',' -f2)
    query_file_name=${queries_dir}/${sample_name}_${query_id}.fasta
    query_output_name=${outdir}/${sample_name}_${query_id}_out.tsv

    # get all (q30) reads supports to singleton variant call, format to a fasta file
    samtools mpileup -r ${query_id} -q30 -x -s -O --output-QNAME -d3000 -f /nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta ${dir_of_bams}/${sample_name}.bam > mpileup_with_readnames/${sample_name}_${query_id}.pileup
    
    # blastn all reads supporting singleton variant call
    blastn -query ${query_file_name} -db nt -evalue 0.1 -num_threads 32 -outfmt "6 qacc sacc staxid bitscore pident" -out "${query_output_name}"
done
