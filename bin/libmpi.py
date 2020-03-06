#!/usr/bin/env python

import os

try:
    from mpi4py import MPI

    MPI_COMM = MPI.COMM_WORLD
    MPI_SIZE = MPI_COMM.Get_size()
    MPI_RANK = MPI_COMM.Get_rank()
    MPI_KING = 0
    MPI_ANY_SOURCE = MPI.ANY_SOURCE
except:
    MPI_COMM = None
    MPI_SIZE = 1
    MPI_RANK = 0
    MPIK_ING = 0
    MPI_ANY_SOURCE = None

if 'CUDA_VISIBLE_DEVICES' in os.environ:
    GPU_s = os.environ['CUDA_VISIBLE_DEVICES'].split(",")
else:
    GPU_s = ['0']
n_GPU = len(GPU_s)
MPI_GPU_BINDING = [None] + ['%d'%(i%n_GPU) for i in range(MPI_SIZE-1)]
