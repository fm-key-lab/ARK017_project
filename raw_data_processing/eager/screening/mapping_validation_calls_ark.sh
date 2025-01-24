#! /bin/bash
# eager calls for mapping of ARK seqs to pestis, pseudotuberculosis


# pestis
tsv_loc=/u/iclight/eager/screening_data
tsv_file_main=ark_screening.tsv

run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/ark_seq_ypestis
mkdir -p ${run_location} && cd ${run_location}
cat ${tsv_loc}/intermediate_files/header.tsv ${tsv_loc}/${tsv_file_main} > ${run_loc}/ark_seqs.tsv
nextflow clean -f
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input './ark_seqs.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --complexity_filter_poly_g \
    --bwaalnn '0.01' \
    --bwaalnl '16' \
    --bwaalnk '2' \
    --run_trim_bam \
    --bamutils_clip_double_stranded_half_udg_left '2' \
    --bamutils_clip_double_stranded_half_udg_right '2' \
    --bamutils_clip_double_stranded_none_udg_left '3' \
    --bamutils_clip_double_stranded_none_udg_right '3' \
    --bamutils_clip_single_stranded_half_udg_left '2' \
    --bamutils_clip_single_stranded_half_udg_right '2' \
    --bamutils_clip_single_stranded_none_udg_left '3' \
    --bamutils_clip_single_stranded_none_udg_right '3' 



# pesudotuberculosis
tsv_loc=/u/iclight/eager/screening_data
tsv_file_main=ark_screening.tsv

run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/ark_seq_ypseudotuberculosis
mkdir -p ${run_location} && cd ${run_location}
cat ${tsv_loc}/intermediate_files/header.tsv ${tsv_loc}/${tsv_file_main} > ${run_location}/ark_seqs.tsv
nextflow clean -f
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input './ark_seqs.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --complexity_filter_poly_g \
    --bwaalnn '0.01' \
    --bwaalnl '16' \
    --bwaalnk '2' \
    --run_trim_bam \
    --bamutils_clip_double_stranded_half_udg_left '2' \
    --bamutils_clip_double_stranded_half_udg_right '2' \
    --bamutils_clip_double_stranded_none_udg_left '3' \
    --bamutils_clip_double_stranded_none_udg_right '3' \
    --bamutils_clip_single_stranded_half_udg_left '2' \
    --bamutils_clip_single_stranded_half_udg_right '2' \
    --bamutils_clip_single_stranded_none_udg_left '3' \
    --bamutils_clip_single_stranded_none_udg_right '3' 

