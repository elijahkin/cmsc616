#!/bin/bash

# Request 8 cores on a single node for 5 minutes
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 00:05:00
#SBATCH -A cmsc416-class
#SBATCH --mem-bind=local


# This is to suppress the warning about not finding a GPU resource
export OMPI_MCA_mpi_cuda_support=0

# Env variable to reduce performance variability
export OMP_PROCESSOR_BIND=true

# Set the number of OpenMP threads to 8
export OMP_NUM_THREADS=8


# Run the executable
./problem1 &> problem1.out
./problem2 &> problem2.out
./problem3 &> problem3.out
./problem4 &> problem4.out
./problem1_serial &> problem1_serial.out
./problem2_serial &> problem2_serial.out
./problem3_serial &> problem3_serial.out
