#!/bin/bash

DATA="../data/phylogeny"
OUT="../data/processed/CorrEvol"
rb_FP="/home/bw566/Applications/revbayes/projects/cmake/build-mpi"

cr1="Hab_full"
cat ../data/processed/CorrEvol/"$cr1".nex | sed 's/23456789//' > ../data/processed/CorrEvol/"$cr1"_states.nex

cr2="FP_hab"
cat ../data/processed/CorrEvol/"$cr2".nex | sed 's/23456789//' > ../data/processed/CorrEvol/"$cr2"_states.nex

# run correlated evolution analysis
rb_command="source(\"./CorrEvol/Corr_Evol_Setup.Rev\");"

echo $rb_command | mpirun -np 2 $rb_FP/rb-mpi

