#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --job-name=admixture
#SBATCH --partition=mpcs_hifmb.p

source ~/.bashrc
conda activate ipyrad

INPREFIX=$1
OUTDIR=$2
KVAL=$3
REP=$4
SEED=$5

cd ${OUTDIR}

cp ${INPREFIX}.bed ./tmp_K${KVAL}_${REP}.bed
cp ${INPREFIX}.bim ./tmp_K${KVAL}_${REP}.bim
cp ${INPREFIX}.fam ./tmp_K${KVAL}_${REP}.fam
admixture -j${SLURM_NTASKS} \
    --cv=10 --seed=${SEED} \
    ./tmp_K${KVAL}_${REP}.bed ${KVAL} \
    | tee tmp_K${KVAL}_${REP}.${KVAL}.log
rm ./tmp_K${KVAL}_${REP}.bed
rm ./tmp_K${KVAL}_${REP}.bim
rm ./tmp_K${KVAL}_${REP}.fam
