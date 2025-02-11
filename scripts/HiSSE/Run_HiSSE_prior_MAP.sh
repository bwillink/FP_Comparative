#!/bin/bash

DATA="../data/phylogeny"
OUT="../data/processed/HiSSE"
rb_FP="/home/bwillink/Applications/revbayes_master/projects/cmake/build-mpi"

cat $DATA/G1_beta_strong.673464_MAP.tree | sed -E 's/Risiocnemis_erythrua/Risiocnemis_erythrura/' | \
                                 sed -E 's/Teinobasis_aeridis/Teinobasis_aerides/' > $OUT/Willink_MAP.tree

# run HiSSE from prior
rb_command="source(\"./HiSSE/HiSSE_setup_MAP_prior_pruned.Rev\");"
echo $rb_command | mpirun -np 1 $rb_FP/rb-mpi
#echo $rb_command | rb
