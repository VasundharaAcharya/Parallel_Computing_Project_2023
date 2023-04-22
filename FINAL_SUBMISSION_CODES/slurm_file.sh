#!/bin/sh


taskset -c 0-159:4 mpirun -N  2 ./mpi_sil
taskset -c 0-159:4 mpirun -N  4 ./mpi_sil

taskset -c 0-159:4 mpirun -N 8 ./mpi_sil

taskset -c 0-159:4 mpirun -N  32 ./mpi_sil


