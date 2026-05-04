#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --job-name=structure
#SBATCH --partition=mpcs_hifmb.p

source ~/.bashrc
conda activate ipyrad

INPATH=$1
INFILE=$2
KVAL=$3
REP=$4
SEED=$5

cd $INPATH

structure -K ${KVAL} -i ${INFILE} -o output${KVAL}_${REP} -D ${SEED}
