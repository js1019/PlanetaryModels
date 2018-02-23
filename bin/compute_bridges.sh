#!/bin/bash
#SBATCH -J nm_model
#SBATCH -p LM
#SBATCH -C LM
#SBATCH --mem 800GB
#SBATCH --time=12:00:00
##SBATCH -A TG-ASC160064
#SBATCH --mail-user=shijia1019@gmail.com
#SBATCH --mail-type=ALL

echo "SLURM_NODELIST"=$SLURM_JOB_NODELIST
cd /home/js116/Earth/modelbuilder
module load matlab
matlab -nodisplay -r "PREM_mesh; quit"
#matlab -nodisplay -r "Gravity; quit"

