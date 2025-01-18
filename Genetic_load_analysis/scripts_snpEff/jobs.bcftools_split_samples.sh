#!/bin/bash
#SBATCH -J bcftools
#SBATCH -o ./logs/bcftools/bcftools_%a_%A.log
#SBATCH -e ./logs/bcftools/bcftools_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 2:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=2000  # Requested Memory
#SBATCH --array=1-10
#SBATCH --mail-type=END

module load bcftools/1.19

WDIR=/nese/meclab/Katya/snpEff_turtles/snpEff_combined/snpEff_out/
SCRATCH=/scratch/workspace/lkomoroske_umass_edu-Leatherback_WGRv2/

IMPACT=$1
VCF=$WDIR/$IMPACT.DerCor_combined_filtered.snpEff.vcf.gz
#SAMPLES=$SCRATCH/Dc_SNPs_May24_samplelist.txt
SAMPLES=$WDIR/failed.smaples.lst

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLES)

bcftools view -s $SAMPLE $VCF | bcftools view -i 'GT="het"' -O z -o $WDIR/split_by_sample/$IMPACT.het.$SAMPLE.snpEff.vcf
bcftools view -s $SAMPLE $VCF | bcftools view -i 'GT="hom"' -O z -o $WDIR/split_by_sample/$IMPACT.hom.$SAMPLE.snpEff.vcf

