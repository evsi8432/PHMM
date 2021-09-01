#!/bin/bash
#SBATCH --account=def-nheckman

#SBATCH --mail-user=evan.sidrow@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

#SBATCH --time=48:00:00
#SBATCH --array=1-4
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
  python Prep_data.py CarHMM
fi
if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then
  python Prep_data.py HHMM
fi
if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then
  python Prep_data.py CarHHMM1
fi
if [ $SLURM_ARRAY_TASK_ID -eq 4 ]
then
  python Prep_data.py CarHHMM2
fi
