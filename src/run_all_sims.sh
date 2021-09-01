#!/bin/bash
#SBATCH --account=def-nheckman

#SBATCH --mail-user=evan.sidrow@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-2000

if [ $SLURM_ARRAY_TASK_ID -lt 501 ]
then
  python Run_Sim.py $SLURM_ARRAY_TASK_ID CarHMM
fi
if [ $SLURM_ARRAY_TASK_ID -gt 500 ] && [ $SLURM_ARRAY_TASK_ID -lt 1001 ]
then
  python Run_Sim.py $SLURM_ARRAY_TASK_ID HHMM
fi
if [ $SLURM_ARRAY_TASK_ID -gt 1000 ] && [ $SLURM_ARRAY_TASK_ID -lt 1501 ]
then
  python Run_Sim.py $SLURM_ARRAY_TASK_ID CarHHMM1
fi
if [ $SLURM_ARRAY_TASK_ID -gt 1500 ]
then
  python Run_Sim.py $SLURM_ARRAY_TASK_ID CarHHMM2
fi
