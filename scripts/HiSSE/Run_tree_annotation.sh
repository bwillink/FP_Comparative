#!bin/bash

rb_FP="/home/bwillink/Applications/revbayes_master/projects/cmake/build-mpi"



rb_command="source(\"./HiSSE/Annot_tree.Rev\");"
echo $rb_command | mpirun -np 1 $rb_FP/rb-mpi
