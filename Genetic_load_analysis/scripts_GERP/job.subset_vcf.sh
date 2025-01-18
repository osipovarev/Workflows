#!/bin/bash
#SBATCH -J randomLines
#SBATCH -o ./logs/randomLines/randomLines_%a_%A.log
#SBATCH -e ./logs/randomLines/randomLines_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 2:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=64000  # Requested Memory

load module bcftools/1.19

VCF=/nese/meclab/Katya/snpEff_turtles/snpEff_combined/snpEff_out/DerCor_combined_filtered.snpEff.vcf.gz

bcftools sort <(cat <(zgrep "^#" $VCF) <(zgrep -v "^#" $VCF | randomLines stdin 100000 stdout)) -O z -o  gerp_less_1.derCor_153.vcf.gz

