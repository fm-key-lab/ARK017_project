# These commands are what were used to execute all nf-core/eager runs 

############################
# Yersinia pestis mapping production run eager calls
############################
# gen72 on its own due to internal barcode requiring additional adapter removal tails
run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_GEN72_production_run
mkdir -p ${run_location} && cd ${run_location}
name_of_run=pestis_evolution_ancient_GEN72
nextflow clean -f

NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    -name ${name_of_run} \
    --input '/u/iclight/eager/pestis_project/ancient_samples_eager_input_GEN72.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --bwaalnn '0.01' \
    --bwaalnl '16' \
    --bwaalnk '2' \
    --run_post_ar_trimming \
    --post_ar_trim_front '7' \
    --post_ar_trim_tail '0' \
    --post_ar_trim_front2 '0' \
    --post_ar_trim_tail2 '0' \
    --run_trim_bam \
    --bamutils_clip_double_stranded_half_udg_left '2' \
    --bamutils_clip_double_stranded_half_udg_right '2' \
    --bamutils_clip_double_stranded_none_udg_left '3' \
    --bamutils_clip_double_stranded_none_udg_right '3' \
    --bamutils_clip_single_stranded_half_udg_left '2' \
    --bamutils_clip_single_stranded_half_udg_right '2' \
    --bamutils_clip_single_stranded_none_udg_left '3' \
    --bamutils_clip_single_stranded_none_udg_right '3' 

## all ancient samples
run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run
mkdir -p ${run_location} && cd ${run_location}
name_of_run=pestis_evolution_production_run
nextflow clean -f
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_max_slurm \
    -name ${name_of_run} \
    --input '/u/iclight/eager/pestis_project/ancient_pestis_non_GEN72.tsv' \
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


run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_modern_production_run
mkdir -p ${run_location} && cd ${run_location}
name_of_run=pestis_evolution_production_run_pseudotb_ref_modern
nextflow clean -f
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    -name ${name_of_run} \
    --input '/u/iclight/eager/pestis_project/modern_pestis.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --mapper 'bwamem' \
    --skip_adapterremoval \
    --skip_fastqc \
    --skip_damage_calculation \
    --bwaalnn '0.01' \
    --bwaalnl '16' \
    --bwaalnk '2' \

## # # # # PROBES TO PESTIS
run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_probes_production_run_pestis
mkdir -p ${run_location} && cd ${run_location}
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input '/u/iclight/eager/pestis_project/pestis_probes.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypestis_ASM906v1/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --complexity_filter_poly_g \
    --bwaalnn '0.01' \
    --bwaalnl '16' \
    --bwaalnk '2' 

######

######
######

######
######

######

####################################################################################################
### Separate mapping to YAC Ypseudotb for comparisson of genomic content, and reference comparissons ###
####################################################################################################

# gen72 on its own due to internal barcode requiring additional adapter removal tails
run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_GEN72_production_run_pseudotb
mkdir -p ${run_location} && cd ${run_location}
name_of_run=pestis_evolution_ancient_GEN72_pseudotb
nextflow clean -f

NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    -name ${name_of_run} \
    --input '/u/iclight/eager/pestis_project/ancient_samples_eager_input_GEN72.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --bwaalnn '0.01' \
    --bwaalnl '16' \
    --bwaalnk '2' \
    --run_post_ar_trimming \
    --post_ar_trim_front '7' \
    --post_ar_trim_tail '0' \
    --post_ar_trim_front2 '0' \
    --post_ar_trim_tail2 '0' \
    --run_trim_bam \
    --bamutils_clip_double_stranded_half_udg_left '2' \
    --bamutils_clip_double_stranded_half_udg_right '2' \
    --bamutils_clip_double_stranded_none_udg_left '3' \
    --bamutils_clip_double_stranded_none_udg_right '3' \
    --bamutils_clip_single_stranded_half_udg_left '2' \
    --bamutils_clip_single_stranded_half_udg_right '2' \
    --bamutils_clip_single_stranded_none_udg_left '3' \
    --bamutils_clip_single_stranded_none_udg_right '3' 

## all ancient samples
run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run_pseudotb
mkdir -p ${run_location} && cd ${run_location}
name_of_run=pestis_evolution_production_run_pseudotb
nextflow clean -f
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_max_slurm \
    -name ${name_of_run} \
    --input '/u/iclight/eager/pestis_project/ancient_pestis_non_GEN72.tsv' \
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


run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_non_GEN72_ancient_production_run_pseudotb
mkdir -p ${run_location} && cd ${run_location}
name_of_run=pestis_evolution_production_run_pseudotb
nextflow clean -f
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_max_slurm \
    --input '/u/iclight/eager/pestis_project/ancient_pestis_singlestranded_samples.tsv' \
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

# all modern samples (plus additional ypseudo samples)
run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_evolution_modern_production_run_pseudotb
mkdir -p ${run_location} && cd ${run_location}
name_of_run=pestis_evolution_modern_production_run_pseudotb
nextflow clean -f
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    -name ${name_of_run} \
    --input '/u/iclight/eager/pestis_project/modern_pestis.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --mapper 'bwamem' \
    --skip_adapterremoval \
    --skip_fastqc \
    --skip_damage_calculation 


## # # # # PROBES TO TB
run_location=/ptmp/iclight/eager/mapping/pestis_evolution_production_runs/pestis_probes_production_run_pseudotb
mkdir -p ${run_location} && cd ${run_location}
NXF_VER=22.10.6 nextflow run nf-core/eager -r 2.4.6 \
    -c /u/iclight/eager/mpcdf_raven.config \
    -profile singularity,mpcdf,raven,resources_normal \
    --input '/u/iclight/eager/pestis_project/pestis_probes.tsv' \
    --fasta '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta' \
    --bwa_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/' \
    --fasta_index '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta.fai' \
    --seq_dict '/nexus/posix0/MPIIB-keylab/reference_genomes/Ypseudotuberculosis_GCF000047365/genome.fasta.dict' \
    --outdir './results' \
    -w './work' \
    --complexity_filter_poly_g \
    --bwaalnn '0.01' \
    --bwaalnl '16' \
    --bwaalnk '2' 
