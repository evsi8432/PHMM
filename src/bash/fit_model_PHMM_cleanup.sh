#!/bin/bash
#SBATCH --account=def-nheckman

#SBATCH --mail-user=evan.sidrow@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

#SBATCH --time=11:59:00
#SBATCH --mem-per-cpu=8G

module load StdEnv/2020
module load gcc/9.3.0
module load r/4.3.1
module load proj/9.0.1
module load gdal/3.5.1
module load udunits/2.2.28

Rscript ../HMM/fit_model_PHMM.R $1
