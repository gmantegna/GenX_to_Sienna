#!/bin/bash

#SBATCH --job-name=ldes         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=4        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=16GB         # memory per cpu-core (4G is default)
#SBATCH --time=23:59:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=fail          # send email when job ends
#SBATCH --mail-user=gm1710@princeton.edu

module purge
module load julia/1.11.1
module load gurobi/12.0.0

julia --project='/home/gm1710/GenX_to_Sienna' GenX_to_Sienna.jl
