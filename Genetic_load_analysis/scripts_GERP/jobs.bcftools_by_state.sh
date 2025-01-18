#!/bin/bash
#SBATCH -J bcftools
#SBATCH -o ./logs/bcftools/bcftools_%a_%A.log
#SBATCH -e ./logs/bcftools/bcftools_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 2:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=2000  # Requested Memory
#SBATCH --array=1-153

module load bcftools/1.19

VCF=$1

SCRATCH=/scratch/workspace/lkomoroske_umass_edu-Leatherback_WGRv2/
SAMPLES=$SCRATCH/Dc_SNPs_May24_samplelist.txt

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLES)

for i in snps indels; \
do \
	hom=$(bcftools view -s $SAMPLE $VCF |  bcftools view -v ${i} -i 'GT!="0/0" && GT="hom"' | wc -l)
	het=$(bcftools view -s $SAMPLE $VCF |  bcftools view -v ${i} -i 'GT="het"' | wc -l )
	echo -e "$SAMPLE\t$i\t$hom\thom"; \
	echo -e "$SAMPLE\t$i\t$het\thet"; \
done 

