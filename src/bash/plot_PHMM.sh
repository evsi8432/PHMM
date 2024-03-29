#!/bin/bash
#SBATCH --account=def-nheckman

#SBATCH --mail-user=evan.sidrow@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

#SBATCH --time=29:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=0-99

module load StdEnv/2020
module load gcc/9.3.0
module load r/4.3.1
module load proj/9.0.1
module load gdal/3.5.1
module load udunits/2.2.28

Rscript ../HMM/plot_PHMM.R $1 $SLURM_ARRAY_TASK_ID
