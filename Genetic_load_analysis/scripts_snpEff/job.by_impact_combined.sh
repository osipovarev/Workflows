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

module load uri/main
module load bcftools/1.19
module load tabixpp/1.1.2-GCC-12.3.0

WDIR=/nese/meclab/Katya/snpEff_turtles/snpEff_combined/snpEff_out
VCF=DerCor_combined_filtered.snpEff.vcf.gz

bcftools view -i 'ANN~"HIGH"' $WDIR/$VCF -O z -o $WDIR/HIGH.$VCF;
tabix $WDIR/HIGH.$VCF;

bcftools view -i 'ANN~"MODERATE"' $WDIR/$VCF -O z -o $WDIR/MODERATE.$VCF;
tabix $WDIR/MODERATE.$VCF;

bcftools view -i 'ANN~"LOW"' $WDIR/$VCF -O z -o $WDIR/LOW.$VCF;
tabix $WDIR/LOW.$VCF;

bcftools view -i 'INFO/ANN~"MODIFIER"' $WDIR/$VCF | grep -v -e "HIGH" -e "MODERATE" -e "LOW" > $WDIR/MODIFIER.$VCF;
bgzip $WDIR/MODIFIER.$VCF;
tabix $WDIR/MODIFIER.$VCF.gz;

