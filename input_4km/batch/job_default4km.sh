#!/bin/bash

#SBATCH --time=0-1:0:0
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=240GB
#SBATCH --switches=1
#SBATCH --output = output_4km


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd ..

srun --cpus-per-task=$SLURM_CPUS_PER_TASK ./mitgcmuv
