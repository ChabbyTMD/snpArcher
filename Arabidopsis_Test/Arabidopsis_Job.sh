#!/bin/bash
# 
#PBS -V
#PBS -l nodes=ram15tb:ppn=48
#PBS -N Arabidopsis_Test
#PBS -joe
#PBS -q batch
#PBS -o Arabidopsis_Test.stdout
#PBS -e Arabidopsis_Test.stderr
#PBS -m abe
#PBS -M tmugoya@sdsu.edu

NCORES=`wc -w < $PBS_NODEFILE`
DATE=`date`
HOST=`hostname`
cd $PBS_O_WORKDIR 
echo " "
echo "running on host: $HOSTNAME"
echo "$NCORES cores requested"
echo "job submitted: $DATE"
echo "job STDOUT follows:"
echo " "

source /home/tmugoya/bashrc-ana3


export PATH=/home/tmugoya/.conda/envs/snparcher/bin/:$PATH
export PATH=/usr/local/anaconda3/bin/conda:$PATH

which snakemake
which conda
conda --version


# # For use on 2TB node

# Dry run
# snakemake -n -s ../workflow/Snakefile -d . --cores 48 --conda-frontend conda --use-conda
# Dry run incomplete 
# snakemake -n -s ../workflow/Snakefile -d . --cores 48 --conda-frontend conda --use-conda --rerun-incomplete

# Wet Run
# snakemake -s ../workflow/Snakefile -d . --cores 48 --conda-frontend conda --use-conda

# For use when workflow is forcibly stopped.
# snakemake -n -s ../workflow/Snakefile -d . --cores 48 --conda-frontend conda --use-conda --unlock

# snakemake -s ../workflow/Snakefile -d . --cores 48 --conda-frontend conda --use-conda --rerun-incomplete

# NEW!!! Use ammended snparcher snakemake profile that specifies resources

# Dry run
# snakemake -n -s ../workflow/Snakefile -d . --cores 48 --workflow-profile ../profiles/default --conda-frontend conda --use-conda
# Dry run incomplete 
# snakemake -n -s ../workflow/Snakefile -d . --cores 48 --workflow-profile profiles/default --conda-frontend conda --use-conda --rerun-incomplete
# Wet run
snakemake -s ../workflow/Snakefile -d . --cores 48 --workflow-profile ../profiles/default --conda-frontend conda --use-conda
# Wet run incomplete
# snakemake -s ../workflow/Snakefile -d . --cores 48  --workflow-profile ../profiles/default --conda-frontend conda --use-conda --rerun-incomplete

# 512 GB node
snakemake -n -s ../workflow/Snakefile -d . --cores 20 --workflow-profile ../profiles/default --conda-frontend conda --use-conda

# Wet run
snakemake -s ../workflow/Snakefile -d . --cores 20 --workflow-profile ../profiles/default --conda-frontend conda --use-conda

# Wet run incomplete
snakemake -s ../workflow/Snakefile -d . --cores 20  --workflow-profile ../profiles/default --conda-frontend conda --use-conda --rerun-incomplete