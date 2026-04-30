#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=256GB
#SBATCH --time=08:00:00
#SBATCH --job-name=lifto
#SBATCH --partition=mpcs_hifmb.p

source ~/.bashrc
conda activate gatk
cd /fs/dss/work/noge4093/apal_popgen/jupyter_notebooks
gatk --java-options '-Xmx64G' LiftoverVcf \
    I=../filtered_calls/allsamples.vcf.gz \
    O=../filtered_calls/allsamples_lifted_to_jaAcrPala1.3.vcf.gz \
    CHAIN=/fs/dss/groups/agmarinecons/apal_reference_panel/apal_reference_panel/liftover_process/chains/adig_vs_apal.chain \
    RECOVER_SWAPPED_REF_ALT=true \
    REJECT=../filtered_calls/allsamples_lifted_to_jaAcrPala1.3_rejected.vcf.gz \
    R=/fs/dss/groups/agmarinecons/apal_reference_panel/apal_reference_panel/references/GCF_964030605.1_jaAcrPala1.3_genomic.fna
