#!/bin/bash -x

module load spectrum-mpi

#####################################################################################################
# Launch N tasks per compute node allocated. Per below this launches 32 MPI rank per compute node.
# taskset insures that hyperthreaded cores are skipped.
#####################################################################################################
taskset -c 0-159:4 mpirun -N 2 /gpfs/u/home/PCPC/PCPCchnr/scratch/project/final-result/a.out