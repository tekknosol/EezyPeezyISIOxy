#! /bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --ntasks=98
#SBATCH --array=0-8
#SBATCH --job-name="isioxy_full"
#SBATCH --output=hpc/output/%x-%A_%a.out
#SBATCH --mail-user=philipp.keller@igb-berlin.de
#SBATCH --mail-type=ALL

cd ~/EezyPeezyISIOxy/
  
module load StdEnv/2020 r/4.1.0

Rscript isioxy_run.R $SLURM_ARRAY_TASK_ID $1 $2