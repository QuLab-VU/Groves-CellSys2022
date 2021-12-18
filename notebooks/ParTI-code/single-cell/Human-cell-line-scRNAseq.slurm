#!/bin/bash
#SBATCH --mail-user=sarah.f.maddox@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=2
#SBATCH --tasks-per-node=8
#SBATCH --mem=50G
#SBATCH --time=20:00:00
#SBATCH --output=matlab_job_slurm.out

module load MATLAB   # load the default version of Matlab

matlab -nodisplay -nosplash < Human-cell-line-scRNAseq.m