#!/bin/bash

# Request 64 cores on a single node for 5 minutes
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 00:05:00
#SBATCH -A cmsc416-class
#SBATCH --mem-bind=local
#SBATCH --exclusive

# This is to suppress the warning about not finding a GPU resource
export OMPI_MCA_mpi_cuda_support=0

# Env variable to reduce performance variability
export OMP_PROCESSOR_BIND=true

declare -a nums_threads=(1 2 4 8 16 32 64)

for n in "${nums_threads[@]}"
do
  # Set the number of OpenMP threads to n
  export OMP_NUM_THREADS=$n

  # Run the executables
  ./problem1 16384 &> problem1_$n.out
  ./problem2 8192 &> problem2_$n.out
  ./problem3 67108864 &> problem3_$n.out
  ./problem4 8192 &> problem4_$n.out
done
