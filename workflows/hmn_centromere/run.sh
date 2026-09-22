#!/bin/bash
#SBATCH --job-name=run_pgkmc_cm
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=5G 
#SBATCH --time=200:00:00      
#SBATCH --output run_pgkmc_cm.log
#SBATCH --mail-type=END
#SBATCH --mail-user=gianfranco@ds.dfci.harvard.edu

source /homes2/gianfranco/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --unlock
snakemake --use-conda --cores 64 --rerun-incomplete