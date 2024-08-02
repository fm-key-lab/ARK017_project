#!/bin/bash -l

# output and error
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
# Job Name:
#SBATCH -J generate_pestis_pseudotb_pangenome_roary
#
# Directory:
#SBATCH -D /ptmp/iclight/pestis_work
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

conda activate roary

roary -p 72 -e -f roary_output_pestis_pseudotb_95 prokka_pestis_pseudotb/*.gff
