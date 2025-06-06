#!/bin/bash

#
# check for clean
#
# https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-do-i-remove-all-files-created-by-snakemake-i-e-like-make-clean
# if [ "$1" == "clean" ]; then
#     echo 'rm $(snakemake --summary | tail -n+2 | cut -f1)'
#     snakemake --summary | tail -n+2 | cut -f1
#     rm -f $(snakemake --summary | tail -n+2 | cut -f1)
#     exit 0
# fi


# launch snakemake to run jobs via SLURM
# ntasks
SM_PARAMS="job-name ntasks time mem mail-user mail-type error output"

#SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task} --mem-per-cpu {cluster.mem-per-cpu-mb}"
SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task}" #new

for P in ${SM_PARAMS}; do SM_ARGS="${SM_ARGS} --$P {cluster.$P}"; done
echo "SM_ARGS: ${SM_ARGS}"

# our SLURM error/output paths expect a logs/ subdir in PWD
mkdir -p logs

conda config --set ssl_verify no

## activate snakemake conda env
activate_conda_path=$(echo $CONDA_EXE | sed 's#bin/conda#bin/activate#g')
source ${activate_conda_path} snakemake

## Check if slurm_status.py is executable and change if not
if [[ ! -x "scripts/slurm_status_script.py" ]]
then
    chmod +x scripts/slurm_status_script.py;
    echo "Changed 'scripts/slurm_status_script.py' to executable";
fi

## Check if email has been updated
email_update=$(grep "UPDATE_EMAIL_ADDRESS" cluster.slurm.json)
[ ! ${#email_update} -eq 0 ] && echo "Please modify cluster.slurm.json to include your email address for error reporting" && exit 1

## create tmp and conda storage directory in /ptmp/ (if not already present)
mkdir -p /ptmp/${USER}/tmp
mkdir -p /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/

### run snakemake
# -j defines total number of jobs executed in parallel on cluster
# -n dryrun
# -p print command lines
# --use-conda allows to activate conda env necessary for rule
# --conda-prefix envs: builds environment in envs folder where it will be 
snakemake -p \
    $* \
     --latency-wait 60 \
    -j 100 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch ${SM_ARGS}" \
    --cluster-status scripts/slurm_status_script.py \
    --default-resources "tmpdir='/ptmp/iclight/tmp'" \
    --group-components var2pos=200 \
    --rerun-incomplete \
    --restart-times 1 \
    --keep-going \
    --use-conda \
    --conda-prefix /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/

    #--dry-run \
    #--group-components \
    # --dag \
    # | dot -Tsvg > dag.svg
