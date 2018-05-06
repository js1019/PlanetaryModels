#!/bin/bash
#SBATCH -J nm_model
#SBATCH -p LM
#SBATCH -C LM
#SBATCH --mem 400GB
#SBATCH --time=4:00:00
##SBATCH -A TG-ASC160064
#SBATCH --mail-user=shijia1019@gmail.com
#SBATCH --mail-type=all

echo "SLURM_NODELIST"=$SLURM_JOB_NODELIST
cd /home/js116/GitHub/PNMsG_models/modelbuilder
module load matlab/R2018a
#matlab -nodisplay -r "CONST_mesh; quit"
matlab -nodisplay -r "Gravity; quit"

