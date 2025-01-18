#!/bin/bash
#SBATCH -J impact_annovar
#SBATCH -o ./logs/impact_annovar/impact_annovar_%a_%A.log
#SBATCH -e ./logs/impact_annovar/impact_annovar_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 6:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=32000  # Requested Memory
#SBATCH --mail-type=END

module load uri/main
module load bcftools/1.19
module load tabixpp/1.1.2-GCC-12.3.0

WDIR=/nese/meclab/Katya/ANNOVAR_turtles/
VCF=derCor_annovar.derCor_multianno.vcf.gz

## HIGH
#bcftools filter -i 'INFO/ExonicFunc.genes ~ "frameshift_insertion" || INFO/ExonicFunc.genes ~ "frameshift_deletion" || INFO/ExonicFunc.genes ~ "frameshift_block_substitution" || INFO/ExonicFunc.genes ~ "stopgain" || INFO/ExonicFunc.genes ~ "stoploss"' $VCF -O z -o $WDIR/HIGH.derCor_annovar.vcf.gz;
#tabix $WDIR/HIGH.derCor_annovar.vcf.gz;

## MODERATE
bcftools filter -i 'INFO/ExonicFunc.genes ~ "nonframeshift_insertion" || INFO/ExonicFunc.genes ~ "nonframeshift_deletion" || INFO/ExonicFunc.genes ~ "nonframeshift_block_substitution" || INFO/ExonicFunc.genes ~ "nonsynonymous_SNV"' $VCF -O z -o $WDIR/MODERATE.derCor_annovar.vcf.gz;
tabix $WDIR/MODERATE.derCor_annovar.vcf.gz;

## LOW
bcftools filter -i 'INFO/ExonicFunc.genes ~ "synonymous_SNV"' $VCF -O z -o $WDIR/LOW.derCor_annovar.vcf.gz;
tabix $WDIR/LOW.derCor_annovar.vcf.gz;

## MODIFIER
bcftools filter -i 'INFO/ExonicFunc.genes ~ "unknown"' $VCF -O z -o $WDIR/MODIFIER.derCor_annovar.vcf.gz;
tabix $WDIR/MODIFIER.derCor_annovar.vcf.gz;
