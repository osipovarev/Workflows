#!/bin/bash
#SBATCH -J WGR_snpEff
#SBATCH -o ./logs/snpEff/WGR_snpEff_%a_%A.log
#SBATCH -e ./logs/snpEff/WGR_snpEff_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 10:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=32000  # Requested Memory
#SBATCH --mail-type=END

module load snpeff/2017-11-24

WDIR=/nese/meclab/Katya/snpEff_turtles/snpEff_combined
CONFIG=$WDIR/snpEff.config
OUTDIR=$WDIR/snpEff_out
VCF=/nese/meclab/Lisa/Dc_WGR_analyses_docsMay2024/Final_filtered_combvcf/DerCor_combined_filtered.vcf
GENINPUT=rDerCor1.pri.cur.20210524

snpEff eff -nodownload -c $CONFIG $GENINPUT $VCF > $OUTDIR//DerCor_combined_filtered.snpEff.vcf

