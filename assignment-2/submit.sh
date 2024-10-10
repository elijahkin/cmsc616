#!/bin/bash

# Allocate 128 cores on a single node for 5 minutes
#SBATCH -N 1
#SBATCH --ntasks=128
#SBATCH -t 00:05:00
#SBATCH -A cmsc416-class
#SBATCH --mem-bind=local
#SBATCH --exclusive

# This is to suppress the warning about not finding a GPU resource
export OMPI_MCA_mpi_cuda_support=0

# Load OpenMPI
module load openmpi/gcc

declare -a nums_procs=(4 8 16 32 64 128)

> life_nonblocking.out
for n in "${nums_procs[@]}"
do
  mpirun -np $n ./life-nonblocking data/life.512x512.data 500 512 512 &>> life_nonblocking.out
done
