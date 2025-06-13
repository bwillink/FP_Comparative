#!/bin/bash

DATA="../data/phylogeny"
OUT="../data/processed/OU"
rb_FP="/home/bw566/Applications/revbayes/projects/cmake/build-mpi"


# run correlated evolution analysis
rb_command="source(\"./OU/OU_Setup.Rev\");"

echo $rb_command | mpirun -np 2 $rb_FP/rb-mpi

