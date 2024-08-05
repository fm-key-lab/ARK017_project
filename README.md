# Online Code/Metadata Repository for ARK017 and reanalysis of *Y. pestis* ancient genomes

## Repository Contents

This repository contains code, input data sheets, metadata, and small datasets used for the manuscript by Light-Maka et al.

Additional questions or comments can be sent to light@mpiib-berlin.mpg.de, or via github issues.

### Directory Overview

- **raw_data_processing**: Contains shell scripts and code for processing raw data inputs for downstream analysis.
  - **eager**: Contains nf-core/eager execution code and scripts for eager setup and conversion from eager outputs to snakemake inputs.
    - **screening**: Contains eager execution code for initial data screening, input TSV, and maltextract species ID list.
  - **gene_content**: Contains code for generating data used in gene-content analysis (*Y. pseudotuberculosis* genes retained in early segments of diversity) after running eager.
- **custom_maltextract_postprocessing**: Extensions of maltextract for screening, used after running eager for screening.
- **snakemake**: Further processing of mapped data to generate candidate mutation matrices for downstream filtering and analysis.
- **local_analysis**: QC and analysis scripts for generating results in the publication, including metadata formatting scripts.
  - **analysis_py**: SNV level analysis.
  - **host_validation**: Generation and plotting of PCA.
- **tree_files**: Additional tree files, including uncollapsed tree, bootstrap support trees, and model parameters.
- **metadata_tables**: Tables from various studies whose samples were reanalyzed, including tables generated using existing metadata (see `local_analysis/formatting*.py`)
  - **pca_tables**: Tables of PCA coordinates, and PCA loadings for plotting
  - **c14_unified_processing**: Outputs (pdf + tables) of unified processing of uncalibrated dates using OxCal 4.4 and IntCal20 atmospheric curve
  - **associated_study_metadata**: Minimally processed metadata tables from studies to identify isolation sources of reanalyzed public modern *Y. pestis* genomes.
- **software_versions**: YAML files for conda environments for (screening post-processing)[software_versions/eager_screening_postprocessing.yml], [local analysis](software_versions/analysis_py_conda_environment.yaml), [other bioinformatic processing of mapped data](software_versions/analysis_py_conda_environment.yaml), (snakemake)[software_versions/snakemake_conda_env.yml] and [nf-core/eager pipeline versions information](software_versions/eager_software_versions.csv).


## Execution dependencies
The analysis done can be roughly organized into 4 sections: Screening, SNV-based analysis on *Y. pestis*, and gene content investigations for *Y. pseudotuberculosis* and *Y. pestis*. The dependencies of analysis is described below, but the 4 sections of analysis are independent (although some sections overalp in their initial steps/can be used across multiple parts).

## PART 1: SCREENING
### Basic overview of execution
1) Process screening data with nf-core/eager: [screening_call_eager.sh](raw_data_processing/eager/screening/screening_call_eager.sh)
2) Reprocess maltextract output (results folder) with call to [post_screening_malt](custom_maltextract_postprocessing/post_screening_malt)

## PART 2: SNV-BASED ANALYSIS
### Basic overview of execution
1) Download data to reanalyze with nf-core/fetchngs: [generate_eager_input.sh](raw_data_processing/eager/generate_eager_input.sh)
2) Run alignment and library-specific aDNA base correction: [pestis_evolution_production_run_to_pestis.sh](raw_data_processing/eager/pestis_evolution_production_run_to_pestis.sh)
3) Collapse all libraries per sample into single bam for snakemake: [merge_eager_samples_for_snakemake.py](raw_data_processing/eager/merge_eager_samples_for_snakemake.py)
4) Generate additional bedtools coverage files for identification of coverage islands and coverage peaks: [raw_data_processing/eager/generate_bedcoverage_files.sh](raw_data_processing/eager/generate_bedcoverage_files.sh)
5) Run snakemake variant calling (run within `snakemake/mapping` directory): [snakemake/mapping/snakemakeslurm.sh](snakemake/mapping/snakemakeslurm.sh)
6) Run snakemake candidate mutation table generation (run within `snakemake/case` directory): [snakemake/case/snakemakeslurm.sh](snakemake/case/snakemakeslurm.sh)
7) Generate formatted metadata tables for analysis with [formatting_lnba_metadata.py](local_analysis/formatting_lnba_metadata.py) and [formatting_isolation_metadata.py](local_analysis/formatting_isolation_metadata.py)
8) Run local analysis initial QC script for filtering mutations and basecalls [pestis_evolution_initial_qc.py](local_analysis/analysis_py/pestis_evolution_initial_qc.py)
   1) Run blast+ on identified ancient singletons [blast_execution_ancient_singletons.sh](local_analysis/analysis_py/blast_execution_ancient_singletons.sh) (line 514 break-point)
   2) Continue filtering
9) Run main analysis script

## PART 3: *Y. pseudotuberculosis* GENE CONTENT
### Basic overview of execution
1) See PART 2, steps 1-2, executing eager with *Y. pseudotuberculosis* reference genome
2) Create bedtools genome and regions file: [create_bedtools_anno_genome_files.sh](raw_data_processing/gene_content/bedtools/create_bedtools_anno_genome_files.sh)
3) Collapse all deduplicated, but unfiltered libraries per sample into single bam and run bedtools: [raw_data_processing/gene_content/bedtools](raw_data_processing/gene_content/bedtools)
4) Generate pangenome information for QC of above generated data:
   1) Download reference genomes from database [download_references_and_format.sh](raw_data_processing/gene_content/pangenome_scripts/roary_pangenome/download_references_and_format.sh)
   2) Run [generate_prokka_annotation.sh](raw_data_processing/gene_content/pangenome_scripts/roary_pangenome/generate_prokka_annotation.sh) on the reference genomes to create Roray compatible .gff files
   3) Run [generate_pestis_pseudotb_pangenome_roary.sh](raw_data_processing/gene_content/pangenome_scripts/roary_pangenome/generate_pestis_pseudotb_pangenome_roary.sh) to generate a pangenome with information on gene presence in each reference genome
5) Analyze with [gene_content_analysis.py](local_analysis/gene_content_analysis.py)

## PART 4: *Y. pestis* GENE CONTENT
### Basic overview of execution
1) See PART 2, steps 1-2, executing eager with *Y. pestis* reference genome
2) Collapse all deduplicated, but unfiltered libraries per sample into single bam and run bedtools: [raw_data_processing/gene_content/bedtools](raw_data_processing/gene_content/bedtools)
3) Analyze with [gene_content_analysis.py](local_analysis/gene_content_analysis.py) (second half of analysis script)
