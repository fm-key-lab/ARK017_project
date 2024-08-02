#!/bin/bash -l

# output and error
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
# Job Name:
#SBATCH -D /home/light/pestis_work/raxml/raxml_final
#SBATCH -J raxml_ng_bootstrapping
#SBATCH --nodes=1
#SBATCH -c 32
#SBATCH --time=48:00:00
#SBATCH --mem=100000


conda activate raxml

echo full tree, 1000 bootstraps

raxml-ng --all --threads 32 --msa strand_1_maf_9_heterozygosityAncient_recombinationAll_aidaexclusionsRepeatRNA_minorafcovper95_full_qc.fa --outgroup YAC --model GTR+G --prefix full_tree_alignment_1000bs_ML_all_production_run_fullqc --bs-trees 1000

raxml-ng --support --tree full_tree_alignment_1000bs_ML_all_production_run_fullqc.raxml.bestTree --bs-trees full_tree_alignment_100bs_ML_all_production_fullqc.raxml.bootstraps --prefix full_tree_alignment_1000bs_ML_all_production_run_fullqc_with_bootstraps