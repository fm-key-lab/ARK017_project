#!/bin/bash -l

# output and error
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
# Job Name:
#SBATCH -J generate_pestis_pseudotb_prokka_annotations
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

conda activate prokka

mkdir prokka_output

for fna in $(ls -1 fastas/*.fna):
    do 
    basename=$(echo ${fna} | cut -d'.' -f1 | cut -d'/' -f2)
    if [ ! -d "prokka_output/$basename" ] 
        then
        prokka --cpus 72 --kingdom Bacteria --outdir prokka_output/${basename} --genus Yersinia --locustag ${basename} ${fna}
    fi
done

cd /ptmp/iclight/pestis_work/prokka_output
ls -1 | while read name; do for f in $(ls -1 ${name}); do basename_f=$(echo $f | awk -F'.' '{print $1}'); ext=$(echo $f | awk -F'.' '{print $2}'); mv ${name}/$f ${name}/${name}.${ext}; done; done