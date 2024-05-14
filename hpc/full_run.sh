#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --mail-user=philipp.keller@igb-berlin.de
#SBATCH --mail-type=ALL

cd ~/EezyPeezyISIOxy/

module load StdEnv/2020 r/4.1.0

Rscript run.R