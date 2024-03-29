#!/bin/sh
#PBS -N sim-0329
#PBS -j oe
#PBS -l ncpus=6
#PBS -q low
#export OMP_NUM_THREADS=6

cd $HOME/sim20240318_us2d_tauV-xi_Vlinit
./All.sh

#$HOME/work/All.sh
